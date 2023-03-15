#!/usr/bin/env nextflow

// NXF ver 19.08+ needed because of the use of tuple instead of set
if( !nextflow.version.matches('>=19.08') ) {
    println "This workflow requires Nextflow version 19.08 or greater and you are running version $nextflow.version"
    exit 1
}

// ANSI escape codes to color output messages
ANSI_GREEN = "\033[1;32m"
ANSI_RED   = "\033[1;31m"
ANSI_RESET = "\033[0m"

// pipeline input parameters 
params.readsdir = null
params.sample_sheet = null
params.fasta = null
params.outdir = "${workflow.launchDir}/results"
params.no_update_freyja = false
params.help = ""

if (params.help || !params.sample_sheet || !params.fasta || !params.readsdir) {
    helpMessage()
    exit(0)
}

log.info """\
 Used parameters:
 ----------------
 --readsdir               : ${params.readsdir}
 --sample_sheet           : ${params.sample_sheet}
 --outdir                 : ${params.outdir}
 --fqpattern              : ${params.fqpattern}
 --fasta                  : ${params.fasta}
 --primerfile             : ${params.primerfile}
 --amplicon_primerfile_v3 : ${params.amplicon_primerfile_v3}
 --amplicon_primerfile_v4 : ${params.amplicon_primerfile_v4}
 --no_update_freyja       : ${params.no_update_freyja}

 Runtime data:
 -------------
 Running with profile     : ${ANSI_GREEN}${workflow.profile}${ANSI_RESET}
 Running as user          : ${ANSI_GREEN}${workflow.userName}${ANSI_RESET}
 Launch dir               : ${ANSI_GREEN}${workflow.launchDir}${ANSI_RESET}
 Base dir                 : ${ANSI_GREEN}${baseDir}${ANSI_RESET}
"""

def helpMessage() {
log.info """\
 Usage:
 ------
 --readsdir               : Directory with input fastq files, required
 --sample_sheet           : The sample sheet of samples to analyse, required
 --outdir                 : Where results will be saved, default is "results" in the working directory
 --fqpattern              : Regex pattern to match fastq files, default is "*R{1,2}*.fastq.gz"
 --fasta                  : Reference genome to use for alignment, required
 --primerfile             : Paired End primer file (bedpe)
 --amplicon_primerfile_v3 : Amplicon primer scheme v3-v4.01
 --amplicon_primerfile_v4 : Amplicon primer scheme v4.02
 --no_update_freyja       : Do not update freyja lineages
 """
}

// Read the sample sheet CSV file cols:
//  - sample_id
//  - sample_site_code
//  - timestamp_sample_collected
//  - sequencing_lab_code
//  - sequencing_sample_id
//  - sequencing_run_id
// This will generate the input folder locations and metadata:
//  - fq_input: ${sequencing_run_id}/${sequencing_sample_id}
//  - metadata: data from other columns as a map (?)

// Add trailing slash
readsdir_repaired = "${params.readsdir}".endsWith("/") ? "${params.readsdir}" : "${params.readsdir}/"
outdir = "${params.outdir}".endsWith("/") ? "${params.outdir}" : "${params.outdir}/"

// Read the input CSV, find the files and make the metadata map
Channel
    .fromPath(params.sample_sheet)
    .splitCsv(header:true, sep:",")
    .map{ create_fastq_channels(it, readsdir_repaired) }
    .set{ ch_input_from_csv }

// Create a channel from the input genome
Channel
    .from(params.fasta)
    .map{ f -> file(f) }
    .set{ ch_input_genome }
// Duplicate 
ch_input_genome.into { ch_input_genome_0; ch_input_genome_1 }

// Create a channel from the PE bed file
Channel
    .from(params.primerfile)
    .map{ f -> file(f) }
    .set{ ch_input_bedpe }

// Create a channel from the amplicon bed file (V3.*, 4.01)
Channel
    .from(params.amplicon_primerfile_v3)
    .map{ f -> file(f) }
    .set{ ch_amplicon_bed_v3 }
// Amplicon bed is used twice, so duplicate channel
ch_amplicon_bed_v3.into { ch_input_amplicon_bed_v3; ch_input_amplicon_bed_v3_1 }

// V4.02 bed channel
Channel
    .from(params.amplicon_primerfile_v4)
    .map{ f -> file(f) }
    .set{ ch_amplicon_bed_v4 }
// Amplicon bed is used twice, so duplicate channel
ch_amplicon_bed_v4.into { ch_input_amplicon_bed_v4; ch_input_amplicon_bed_v4_1 }

// Create a channel from the variant table
Channel
    .from(params.variant_table)
    .map{ f -> file(f) }
    .set{ ch_variant_table }

// Create a channel from the SNP list
Channel
    .from(params.snp_list)
    .map{ f -> file(f) }
    .set{ ch_snps }

// Create a channel for ANNOVAR inputs
Channel
    .from(params.annovar_covdb)
    .map{ f -> file(f) }
    .set{ ch_annovar_covdb }

// Create a channels lineage maps
Channel
    .from(params.who_map)
    .map{ f -> file(f) }
    .set{ ch_who_map }
Channel
    .from(params.sublineage_map)
    .map{ f -> file(f) }
    .set{ ch_sub_map }

// Workflow steps:
//  - fastp;    trimming, mergeing, qcing data
//  - bwa_mem;  align merged fastq
//  - samtools; merge, sort, index
//  - samtools ampliconstats; get coverage
//  - mosdepth; genome coverage
//  - freyja 

// Outputs:
//  - merged, trimmed, qc'd fastq
//  - merged, trimmed, sorted, indexed bams
//  - variants from various pipelines
//  - genome/amplicon coverage
// ========================================================================= //

process collapse_lanes {
    cpus 8
    tag "${meta.lab}:${meta.run_id}:${meta.sequencing_sample_id}"

    input:
        tuple val(meta), file(fq1), file(fq2) from ch_input_from_csv
     output:
        tuple val(meta), file("*_1.collapsed.fastq.gz"), file("*_2.collapsed.fastq.gz") into collapsed_reads

    script:
        def split_cpus = Math.round(Math.floor(task.cpus/2))
        """
        pigz -${split_cpus} -cd ${fq1} | pigz -${split_cpus} > ${meta.sequencing_sample_id}_1.collapsed.fastq.gz
        pigz -${split_cpus} -cd ${fq2} | pigz -${split_cpus} > ${meta.sequencing_sample_id}_2.collapsed.fastq.gz
        """
}

process fastp {
    cpus 2
    tag "${meta.lab}:${meta.run_id}:${meta.sequencing_sample_id}"

    publishDir "${params.outdir}", 
        mode: 'copy', 
        pattern: '*.trim.fastq.gz',
        saveAs: { filename -> saveFile(filename, meta, "trimmed") }

    input:
        tuple val(meta), file(fq1), file(fq2) from collapsed_reads
    output:
        file("*.json")                             into fastp_json    // multiqc
        tuple val(meta), file("*.trim.fastq.gz")   into trimmed_reads // alignment

    script:
        """
        fastp \\
            --in1 ${fq1} \\
            --in2 ${fq2} \\
            --out1 ${meta.sequencing_sample_id}_1.trim.fastq.gz \\
            --out2 ${meta.sequencing_sample_id}_2.trim.fastq.gz \\
            --thread ${task.cpus} \\
            --json ${meta.sequencing_sample_id}.fastp.json \\
            ${params.modules.fastp.args}
        """
}

process bwa_index {
    tag { genome.baseName }

    input:
        file genome from ch_input_genome_0

    output:
        tuple file("ref.fa"), file("ref.fa.*") into ch_index

    script:
        """
        ln -s ${genome} ref.fa
        bwa index ref.fa
        """
}

process bwa_mem {
    cpus 4
    tag "${meta.lab}:${meta.run_id}:${meta.sequencing_sample_id}"

    input:
        tuple val(meta), file(reads), file(fasta), file(ref_star) from trimmed_reads.combine(ch_index)

    output:
        tuple val(meta), file("*.bam") into unsorted_bams

    script:
    def split_cpus = Math.round(Math.floor(task.cpus/2))
        """
        bwa mem \\
            -t ${split_cpus} \\
            $fasta \\
            $reads \\
            | samtools view -@ ${split_cpus} -bh -o ${meta.sequencing_sample_id}.bam -
        """
}

process samtools_sort_idx {
    cpus 4
    tag "${meta.lab}:${meta.run_id}:${meta.sequencing_sample_id}"

    publishDir "${params.outdir}", 
        mode: "copy", 
        pattern: "*.sorted.bam*",
        saveAs: { filename -> saveFile(filename, meta, "alignments") }

    input:
        tuple val(meta), file(bam) from unsorted_bams

     output:
        tuple val(meta), file("*.sorted.bam"), file("*.sorted.bam.bai") into sorted_bams

    script:
        """
        samtools sort \\
            -@ ${task.cpus} \\
            -o ${meta.sequencing_sample_id}.sorted.bam \\
            ${meta.sequencing_sample_id}.bam
        samtools index \\
            ${meta.sequencing_sample_id}.sorted.bam \\
            ${meta.sequencing_sample_id}.sorted.bam.bai
        """
}

process iVar_trim {
    cpus 2
    tag "${meta.lab}:${meta.run_id}:${meta.sequencing_sample_id}"

    publishDir "${params.outdir}", 
        mode: "copy", 
        pattern: "*.iVar_trimmed.sorted.bam*",
        saveAs: { filename -> saveFile(filename, meta, "ivar") }

    input:
        tuple val(meta), file(bam), file(bai), file(bed_v3), file(bed_v4) from sorted_bams.combine(ch_input_amplicon_bed_v3_1).combine(ch_input_amplicon_bed_v4_1)

    output:
        tuple val(meta), file("*.iVar_trimmed.sorted.bam"), file("*.iVar_trimmed.sorted.bam.bai") into trimmed_bams, trimmed_bams_1, trimmed_bams_2, trimmed_bams_3

    script:
        // Maybe use samtools view -F4 to filter bams before iVar (remove unmapped reads)
        if (meta.primer_version == "4.02")
            """
            ivar trim -e \\
                -i ${bam} \\
                -b ${bed_v4} \\
                -p ${meta.sequencing_sample_id}.iVar_trimmed \\
                ${params.modules.ivar.args}
            samtools sort \\
                -@ ${task.cpus} \\
                -o ${meta.sequencing_sample_id}.iVar_trimmed.sorted.bam \\
                ${meta.sequencing_sample_id}.iVar_trimmed.bam
            samtools index \\
                ${meta.sequencing_sample_id}.iVar_trimmed.sorted.bam \\
                ${meta.sequencing_sample_id}.iVar_trimmed.sorted.bam.bai
            """
        else

            """
            ivar trim -e \\
                -i ${bam} \\
                -b ${bed_v3} \\
                -p ${meta.sequencing_sample_id}.iVar_trimmed \\
                ${params.modules.ivar.args}
            samtools sort \\
                -@ ${task.cpus} \\
                -o ${meta.sequencing_sample_id}.iVar_trimmed.sorted.bam \\
                ${meta.sequencing_sample_id}.iVar_trimmed.bam
            samtools index \\
                ${meta.sequencing_sample_id}.iVar_trimmed.sorted.bam \\
                ${meta.sequencing_sample_id}.iVar_trimmed.sorted.bam.bai
            """
}

process freyja_update {
    cpus 1
    tag "${update_flag}"
    cache false

    output: 
        val updated_ok into freyja_updated

    script:
        updated_ok = true
        update_flag = params.no_update_freyja ? "skipped" : "updated"
        if (params.no_update_freyja)
            """
            """
        else
            """
            freyja update
            """
}

process freyja_variants {
    cpus 4
    tag "${meta.lab}:${meta.run_id}:${meta.sequencing_sample_id}"

    publishDir "${params.outdir}", 
        mode: "copy", 
        pattern: "${meta.sequencing_sample_id}*",
        saveAs: { filename -> saveFile(filename, meta, "freyja") }

    input:
        tuple val(meta), file(bam), file(bai), val(freyja_update_flag) from trimmed_bams.combine(freyja_updated)

    output:
        tuple val(meta), file("*variants.tsv"), file("*depths") into freyja_outputs

    script:
        """
        freyja variants \\
            $bam \\
            --variants ${meta.sequencing_sample_id}.variants \\
            --depths ${meta.sequencing_sample_id}.depths 
        """
}

process freyja_demix {
    cpus 1
    tag "${meta.lab}:${meta.run_id}:${meta.sequencing_sample_id}"

    publishDir "${params.outdir}", 
        mode: "copy", 
        pattern: "*.demix.txt",
        saveAs: { filename -> saveFile(filename, meta, "freyja") }

    input:
        tuple val(meta), file(variants), file(depths) from freyja_outputs

    output:
        tuple val(meta), file("*demix.txt") into demix, demix_1

    // NOTE: `freyja demix' may fail, this is most likely due `freyja variants' 
    //       producing empty outputs. We can test the size of the `depths' file 
    //       and length of the `variants' file. If either fail checks we should 
    //       use a template file instead.
    script:
        """
        freyja demix \\
            ${variants} \\
            ${depths} \\
            --output ${meta.sequencing_sample_id}.demix.txt \\
            || sed "s/VARIANT_FILE/${variants}/" ${baseDir}/static/demix_template.txt > ${meta.sequencing_sample_id}.demix.txt
        """

        // if (depths.size() > 0 && variants.countLines() > 1) 
        //     """
        //     freyja demix \\
        //         ${variants} \\
        //         ${depths} \\
        //         --output ${meta.sequencing_sample_id}.demix.txt
        //     """
        // else
        //     """
        //     sed "s/VARIANT_FILE/${variants}/" ${baseDir}/static/demix_template.txt > ${meta.sequencing_sample_id}.demix.txt
        //     """
}

process samtools_amplicon_stats {
    cpus 1
    tag "${meta.lab}:${meta.run_id}:${meta.sequencing_sample_id}"

    publishDir "${params.outdir}", 
        mode: "copy", 
        pattern: "${meta.sequencing_sample_id}*",
        saveAs: { filename -> saveFile(filename, meta, "ampliconstats") }

    input:
        tuple val(meta), file(bam), file(bai), file(amplicon_bed_v3), file(amplicon_bed_v4) from trimmed_bams_1.combine(ch_input_amplicon_bed_v3).combine(ch_input_amplicon_bed_v4)

    output:
        tuple val(meta), file("*.ampliconstats.txt") into ampliconstats

    script:
        if (meta.primer_version == "4.02")
            """
            samtools ampliconstats -o ${meta.sequencing_sample_id}.ampliconstats.txt $amplicon_bed_v4 $bam 
            """
        else
            """
            samtools ampliconstats -o ${meta.sequencing_sample_id}.ampliconstats.txt $amplicon_bed_v3 $bam 
            """
}

// TODO: Would be good to (somehow) group by the `run_id' parameter in the 
//       meta hashmap. This would allow us to output per-run plots.
// process plot_ampliconstats {
//     cpus 1
//     tag "${meta.lab}:${meta.run_id}:${meta.sequencing_sample_id}"
// 
//     publishDir "${params.outdir}", 
//         mode: "copy", 
//         pattern: "${meta.sequencing_sample_id}*",
//         saveAs: { filename -> saveFile(filename, meta, "plot_ampliconstats") }
// 
//     input:
//         tuple val(meta), file(bam), file(bai) from trimmed_bams_3
//     output:
//         tuple val(meta), 
// }

process run_varscan {
    errorStrategy "ignore"
    cpus 2
    tag "${meta.lab}:${meta.run_id}:${meta.sequencing_sample_id}"

    publishDir "${params.outdir}", 
        mode: "copy", 
        pattern: "${meta.sequencing_sample_id}*",
        saveAs: { filename -> saveFile(filename, meta, "varscan") }

    input:
        tuple val(meta), file(bam), file(bai), file(ref) from trimmed_bams_3.combine(ch_input_genome_1)
    output:
        // tuple val(meta), file("*.pileup"), file("*_varscan.tsv"), file("*.acgt.tsv"), file("*_varscan.indel.tsv") into ampliconstats
        tuple val(meta), file("*.pileup"), file("*_varscan.tsv"), file("*_varscan.indel.tsv") into varscan_outputs, varscan_outputs_1, varscan_outputs_2

    // TODO: abstract samtools mpileup to another process
    script:
        """
        echo ${meta.sequencing_sample_id}
        samtools mpileup -f ${ref} -q 10 -d 1000000 ${bam} > ${meta.sequencing_sample_id}.pileup
        varscan pileup2snp ${meta.sequencing_sample_id}.pileup -p-value 0.05 > ${meta.sequencing_sample_id}_varscan.tsv;
        varscan pileup2indel ${meta.sequencing_sample_id}.pileup -p-value 0.05 > ${meta.sequencing_sample_id}_varscan.indel.tsv;
        """
}

process run_varscan_cov {
    errorStrategy "ignore"
    cpus 1
    tag "${meta.lab}:${meta.run_id}:${meta.sequencing_sample_id}"

    conda "bioconda::sequenza-utils==3.0.0"

    publishDir "${params.outdir}", 
        mode: "copy", 
        pattern: "${meta.sequencing_sample_id}*",
        saveAs: { filename -> saveFile(filename, meta, "varscan") }

    input:
        tuple val(meta), file(pileup), file(varscan_tsv), file(varscan_indel_tsv) from varscan_outputs_1
    output:
        tuple val(meta), file("*_cov.tsv"), file("*_acgt.tsv") into varscan_cov

    script:
        """
        awk 'BEGIN {FS="\t";OFS="\t";print "#chr","position","refbase","coverage"} {print \$1,\$2,\$3,\$4}' ${pileup} > ${meta.sequencing_sample_id}_cov.tsv
        sequenza-utils pileup2acgt -p ${pileup} > ${meta.sequencing_sample_id}_acgt.tsv
        """
}

process run_varscan_mpileup {
    errorStrategy "ignore"
    cpus 2
    tag "${meta.lab}:${meta.run_id}:${meta.sequencing_sample_id}"

    publishDir "${params.outdir}", 
        mode: "copy", 
        pattern: "${meta.sequencing_sample_id}*",
        saveAs: { filename -> saveFile(filename, meta, "varscan_mpileup") }

    input:
        tuple val(meta), file(pileup), file(varscan_tsv), file(varscan_indel_tsv) from varscan_outputs
    output:
        tuple val(meta), file("*_varscan.vcf"), file("*_varscan.indel.vcf") into mpileup_varscan_outputs, mpileup_varscan_outputs_1

    script:
        """
        varscan mpileup2snp ${pileup} --min-var-freq 0.01 -p-value 0.05 --output-vcf > ${meta.sequencing_sample_id}_varscan.vcf;
        varscan mpileup2indel ${pileup} --min-var-freq 0.01 -p-value 0.05 --output-vcf > ${meta.sequencing_sample_id}_varscan.indel.vcf;
        """
}


process mosdepth {
    cpus 4
    tag "${meta.lab}:${meta.run_id}:${meta.sequencing_sample_id}"

    publishDir "${params.outdir}", 
        mode: "copy", 
        pattern: "${meta.sequencing_sample_id}*",
        saveAs: { filename -> saveFile(filename, meta, "mosdepth") }

    input:
        tuple val(meta), file(bam), file(bai) from trimmed_bams_2

    output:
        tuple val(meta), file("*.mosdepth.summary.txt") into md_summary

    script:
        """
        mosdepth --no-per-base --threads $task.cpus ${meta.sequencing_sample_id} $bam
        """
}

process join_ampstats_freyja_mosdepth {
    cpus 1
    tag "${meta.lab}:${meta.run_id}:${meta.sequencing_sample_id}"
    errorStrategy 'ignore'

    conda "python=3.9"

    publishDir "${params.outdir}", 
        mode: "copy", 
        pattern: "${meta.run_id}_${meta.sequencing_sample_id}*",
        saveAs: { filename -> saveFile(filename, meta, "stats_csv") }

    input:
        tuple val(meta), file(ampstats), file(mdsummary), file(demix) from ampliconstats.join(md_summary).join(demix)

    output:
        file("*.joined.csv") into joined_stats_vars

    script:
    def meta_sample_id = meta.sample_id ? "--meta-sample-id ${meta.sample_id}" : "" 
    def meta_sample_site_code = meta.sample_site_code ? "--meta-sample-site-code ${meta.sample_site_code}" : "" 
    def meta_timestamp_sample_collected = meta.timestamp ? "--meta-timestamp ${meta.timestamp}" : "" 
    def meta_sequencing_lab_code = meta.sequencing_lab_code ? "--meta-sequencing-lab-code ${meta.sequencing_lab_code}" : "" 
    def meta_original_sample_id = meta.original_sample_id ? "--meta-original-sample-id ${meta.original_sample_id}" : "" 
    def meta_sequencing_sample_id = meta.sequencing_sample_id ? "--meta-sequencing-sample-id ${meta.sequencing_sample_id}" : "" 
    def meta_sequencing_run_id = meta.run_id ? "--meta-sequencing-run-id ${meta.run_id}" : "" 
        """
        collate_results.py ${demix} ${ampstats} ${mdsummary} \\
            --output ${meta.run_id}_${meta.sequencing_sample_id}.joined.csv \\
            ${meta_sample_id} \\
            ${meta_sample_site_code} \\
            ${meta_timestamp_sample_collected} \\
            ${meta_sequencing_lab_code} \\
            ${meta_original_sample_id} \\
            ${meta_sequencing_sample_id} \\
            ${meta_sequencing_run_id} \\
            ${params.modules.agg_outputs.args}
        """
}

process annovar {
    cpus 1
    tag "${meta.lab}:${meta.run_id}:${meta.sequencing_sample_id}"

    conda "conda-forge::perl=5.32.1 conda-forge::sed=4.7"

    publishDir "${params.outdir}", 
        mode: "copy", 
        pattern: "${meta.sequencing_sample_id}*",
        saveAs: { filename -> saveFile(filename, meta, "variant_annotation") }

    input:
        tuple val(meta), file(snp_vcf), file(indel_vcf), file(annovar_covdb) from mpileup_varscan_outputs_1.combine(ch_annovar_covdb)

    output:
        file("${meta.sequencing_sample_id}.annovar*") into annovar_output

    script:
        """
        convert2annovar.pl -format vcf4 ${snp_vcf} > ${meta.sequencing_sample_id}.annovar;
        sed -i s/NC_045512.2/NC_045512v2/ ${meta.sequencing_sample_id}.annovar;
        table_annovar.pl -buildver NC_045512v2 ${meta.sequencing_sample_id}.annovar ${annovar_covdb}/ -protocol avGene -operation g;
        """
}

process join_varscan {
    cpus 1
    tag "${meta.lab}:${meta.run_id}:${meta.sequencing_sample_id}"

    conda "python=3.9 pandas=1.4.1"

    publishDir "${params.outdir}", 
        mode: "copy", 
        pattern: "${meta.run_id}_${meta.sequencing_sample_id}*",
        saveAs: { filename -> saveFile(filename, meta, "varscan") }

    input:
        tuple val(meta), file(coverage), file(acgt), file(pileup), file(varscan_tsv), file(indel), file(vartab), file(snps) from varscan_cov.join(varscan_outputs_2).combine(ch_variant_table).combine(ch_snps)

    output:
        file("*.varsites.tsv") into joined_varscan

    script:
        """
        collate_varscan.py ${vartab} ${coverage} ${acgt} ${indel} ${snps} \\
            --output ${meta.run_id}_${meta.sequencing_sample_id}.varsites.tsv \\
            --sample ${meta.run_id}_${meta.sequencing_sample_id} \\
            -vvv 
        """
}

process all_lineages {
    cpus 1
    tag "${meta.lab}:${meta.run_id}:${meta.sequencing_sample_id}"
    errorStrategy 'ignore'

    conda "python=3.9"

    publishDir "${params.outdir}", 
        mode: "copy", 
        pattern: "${meta.run_id}_${meta.sequencing_sample_id}*",
        saveAs: { filename -> saveFile(filename, meta, "stats_csv") }

    input:
        tuple val(meta), file(demix) from demix_1

    output:
        file("*.all_lineages.csv") into all_lineages

    script:
    def meta_sample_id = meta.sample_id ? "--meta-sample-id ${meta.sample_id}" : "" 
    def meta_sample_site_code = meta.sample_site_code ? "--meta-sample-site-code ${meta.sample_site_code}" : "" 
    def meta_timestamp_sample_collected = meta.timestamp ? "--meta-timestamp ${meta.timestamp}" : "" 
    def meta_sequencing_lab_code = meta.sequencing_lab_code ? "--meta-sequencing-lab-code ${meta.sequencing_lab_code}" : "" 
    def meta_original_sample_id = meta.original_sample_id ? "--meta-original-sample-id ${meta.original_sample_id}" : "" 
    def meta_sequencing_sample_id = meta.sequencing_sample_id ? "--meta-sequencing-sample-id ${meta.sequencing_sample_id}" : "" 
    def meta_sequencing_run_id = meta.run_id ? "--meta-sequencing-run-id ${meta.run_id}" : "" 
        """
        collate_lineages.py ${demix} \\
            --output ${meta.run_id}_${meta.sequencing_sample_id}.all_lineages.csv \\
            ${meta_sample_id} \\
            ${meta_sample_site_code} \\
            ${meta_timestamp_sample_collected} \\
            ${meta_sequencing_lab_code} \\
            ${meta_original_sample_id} \\
            ${meta_sequencing_sample_id} \\
            ${meta_sequencing_run_id} \\
            -vvv
        """
}

process aggregate_joined {
    cpus 1
    // https://github.com/nextflow-io/nextflow/issues/1385
    beforeScript 'ulimit -Ss unlimited'

    publishDir "${params.outdir}", mode: "copy"

    input:
        file joined_csvs from joined_stats_vars.collect()

    output:
        file("aggregated.csv")

    script:
        """
        xsv cat rows ${joined_csvs} > aggregated.csv
        """
}

process aggregate_lineages {
    cpus 1
    // https://github.com/nextflow-io/nextflow/issues/1385
    beforeScript 'ulimit -Ss unlimited'

    publishDir "${params.outdir}", mode: "copy"

    input:
        file lineages from all_lineages.collect()

    output:
        file("all_lineages.csv") into ch_all_lineages

    script:
        """
        xsv cat rows ${lineages} > all_lineages.csv
        """
}

process apply_WHO_labels {
    cpus 1
    publishDir "${params.outdir}", mode: "copy"

    conda "conda-forge::python=3.9 conda-forge::pandas=1.4.2"

    input:
        file lineages from ch_all_lineages
        file who_map from ch_who_map
        file sub_map from ch_sub_map

    output:
        file("all_lineages_with_WHO_labels.csv")

    script:
        """
        apply_who_labels.py ${lineages} ${who_map} ${sub_map} --output all_lineages_with_WHO_labels.csv -vvv
        """
}

process aggregate_varscan {
    cpus 1
    // https://github.com/nextflow-io/nextflow/issues/1385
    beforeScript 'ulimit -Ss unlimited'

    publishDir "${params.outdir}", mode: "copy"

    input:
        file varsites from joined_varscan.collect()

    output:
        file("all_varsites.csv")

    script:
        """
        xsv cat rows ${varsites} > all_varsites.csv
        """
}

workflow.onComplete {
    if (workflow.success) {
        log.info """
            ${ANSI_GREEN}Finished successfully${ANSI_RESET}
            """
            .stripIndent()
    } else {
        log.info """
            ${ANSI_RED}Finished with errors!${ANSI_RESET}
            """
            .stripIndent()
    }
}

// ========================================================================= //
// FUNCTIONS 
// ========================================================================= //

/**
 * Convert CSV sample sheet to an array [ meta, fastq_1, fastq_2 ]
 *
 * row fields (sample sheet cols):
 *  - sample_id
 *  - sample_site_code (optional)
 *  - timestamp_sample_collected (optional)
 *  - sequencing_lab_code
 *  - sequencing_sample_id
 *  - sequencing_run_id
 * 
 * readsdir:
 *   Top level directory path that contains FASTQ under the structure
 *     readsdir/lab/run_id/sequencing_sample_id
 *
 * Returns: array
 *   [ meta (map), fastq_1 (file), fastq_2 (file) ]
 **/
def create_fastq_channels(LinkedHashMap row, readsdir) {
    def NAME_MAP = [:]
    NAME_MAP.nottingham_uni = "nottingham"
    NAME_MAP.ea = "nottingham"
    NAME_MAP.liverpool_uni = "liverpool"
    NAME_MAP.exeter_uni = "exeter"

    def meta = [:]
    meta.sample_id                     = row.sample_id
    meta.run_id                        = row.sequencing_run_id
    meta.sample_site_code              = row.sample_site_code ?: ""
    meta.timestamp                     = row.timestamp_sample_collected ?: "" // TODO: Can't be permissive if selecting date_collected
    meta.sequencing_lab_code           = row.sequencing_lab_code
    meta.lab                           = NAME_MAP[row.sequencing_lab_code.toLowerCase()] ?: row.sequencing_lab_code.toLowerCase()
    meta.original_sample_id            = row.sequencing_sample_id
    meta.sequencing_sample_id          = fixSampleID(meta.lab, row.sequencing_sample_id)
    meta.sequencing_primer_set_version = row.sequencing_primer_set_version
    meta.primer_version                = fixPrimerVersion(row.sequencing_primer_set_version) 
    meta.input_file_pattern            = getInputFilePattern(params.fqpattern, meta, readsdir)

    def array = []
    array = splitReads([ meta, file(meta.input_file_pattern) ])
    return array
}

/**
 * Get files from sample sheet row
 *
 * fq_pattern:
 *   Glob style pattern for finding files, default "*R{1,2}*.fastq.gz"
 * meta:
 *   Metadata map from the sample sheet row
 * readsdir:
 *   Path to reads directory
 *
 * Returns: str
 *   Reads glob pattern
 **/
def getInputFilePattern(fq_pattern, meta, readsdir) {
    return "${readsdir}${meta['lab']}/${meta['run_id']}/${meta['sequencing_sample_id']}${fq_pattern}"
}

/**
 * Fix sequencing_sample_id for each center, usually replacing single characters
 * 
 * lab:
 *   Lab name
 * sample_id:
 *   sequencing_sample_id to be corrected
 *
 * Returns: str
 *   Corrected sample_id if a rule is defined for the lab, otherwise orignal sample_id
 **/
def fixSampleID(lab, sample_id) {
    if (lab == "nottingham") {
        return sample_id.replace("_", "-")
    } else if (lab == "liverpool") {
        return sample_id.replace(".", "_")
    } else {
        return sample_id
    }
}

/**
 * Fix primer set version, removing "v" and "V" and parsing as a Float
 * 
 * primer_str:
 *   Sequencing primer set version as a string
 *
 * Returns: Float
 *   Corrected primer version 
 **/
def fixPrimerVersion(primer_str) {
    primer_str = primer_str - "v"  // Remove lowercase `v'
    primer_str = primer_str - "V"  // Remove uppercase `V'
    // return Float.parseFloat(primer_str)
    return primer_str
}

/**
 * Function to group FASTQs based on read1 or read2
 *
 * entry: 
 *   array of [ meta, fastq ]
 *
 * Returns:
 *   array of [ meta, fastq_1, fastq_2 ]
 **/
def splitReads(entry) {
    def meta = entry[0]
    def fqs = entry[1].flatten()
    def r1 = fqs.findAll { it =~ /_R1_/ }.sort()
    def r2 = fqs.findAll { it =~ /_R2_/ }.sort()
    return [ meta, r1, r2 ]
}

/**
 * Function to create the output directory structure from meta and a final
 * folder name
 * 
 * filename:
 *   output filename
 * meta:
 *   Metadata map from the sample sheet row, must contain `lab` and `run_id`
 * final_dir:
 *   directory where output will be copied/moved to
 **/
def saveFile(filename, meta, final_dir) {
    return "${meta.lab}/${meta.run_id}/${final_dir}/${filename}"
}
