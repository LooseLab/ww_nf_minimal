Wastewater analysis
===================

Introduction
============

**ww_minimal** is a bioinformatics analysis pipeline used to perform the initial quality control and variant analysis on wastewater sequencing samples.
This pipeline supports Illumina short-reads prepared using the [Nimagen primer scheme][nimagen-primers] on various platforms (NovaSeq, NextSeq, MiSeq).


Pipeline summary
----------------

1. Merge sequencing FASTQ files ([`pigz`][pigz])
1. Adapter trimming ([`fastp`][fastp])
1. Variant calling
    1. Read alignment ([`bwa mem`][bwa-mem])
    1. Sort and index alignments ([`Samtools`][samtools])
    1. Primer sequence removal ([`BAMClipper`][bamclipper])
    1. Genome-wide and amplicon coverage ([`mosdepth`][mosdepth], [`Samtools ampliconstats`][ampliconstats])
    1. Variant calling ([`freyja variants/demix`][freyja]; samples may fail on this step due to low coverage, these are omitted from further analysis in the pipeline, they are not excluded overall)
    1. Extract WHO and pango lineages ([`collate_results.py`][collate_results], [`collate_lineages.py`][collate_lineages])
    1. Aggregate all sample outputs ([`xsv`][xsv])

Quickstart
==========

This pipeline uses [`conda`][conda] for environment and package management (recommended to use [`miniconda`][miniconda]).

Initialise environment
----------------------

With [`[mini]conda`][miniconda] installed:

```bash
git clone https://github.com/LooseLab/ww_nf_minimal
cd ww_nf_minimal
conda env create -f environment.yml
```

Run test profile
----------------

```bash
conda activate ww_minimal
nextflow run main.nf -profile test
```

Running an actual run
---------------------

After successfully running the test subset you can attempt to run on other samples.
Read the [input](#input) section for how to setup the FASTQ directory and sample sheet.
Once these are done the pipeline can be run like so:

```bash
nextflow run /path/to/main.nf --readsdir <FASTQ INPUT DIRECTORY> --sample_sheet <SAMPLE SHEET CSV> -with-report report.html
```

If nextflow crashes while running, you can add the flag `-resume` to the previous command to check for cached jobs so the entire pipeline does not need to be re-run.


Input and Output
================

### Input

There are two required user supplied inputs, the sample sheet and the FASTQ reads directory.
These can be supplied by editing the [`nextflow.config`][nf-config] file adding the `sample_sheet` and `readsdir` attributes to the `params` or supplied on the command line using `--sample_sheet` and `--readsdir`.
In addition there are three static inputs that are provided with the workflow (these may change in the future as the primer scheme changes).
These are the reference genome, paired-end primer file, and amplicon primer file.

#### FASTQ

This pipeline expects FASTQ files to be structured inside an input directory with subfolders for each sequencing lab and then further subfolders for each run ID.
For most labs the share directory can be used directly, however samples from [Exeter require symlinking][symlink-help].
An example input directory structure can be seen below:

```none
input
├── <LAB1>
│  ├── <RUN1>
│  │  ├── SAMPLE_R1_L002_001.fastq.gz
│  │  └── SAMPLE_R2_L002_001.fastq.gz
│  └── <RUN2>
│     ├── SAMPLE_R1_L002_001.fastq.gz
│     └── SAMPLE_R2_L002_001.fastq.gz
└── <LAB2>
   └── <RUN1>
      ├── SAMPLE_L001_R1_001.fastq.gz
      ├── SAMPLE_L001_R2_001.fastq.gz
      ├── ...
      ├── SAMPLE_L004_R1_001.fastq.gz
      └── SAMPLE_L004_R2_001.fastq.gz
```

### Sample sheet CSV

The CSV sample sheet is required as this informs the pipeline which samples should be analysed.
It currently requires 6 fields:
1. `sample_id`
1. `sample_site_code`
1. `timestamp_sample_collected`
1. `sequencing_lab_code`
1. `sequencing_sample_id`
1. `sequencing_run_id`

These are used to find the input FASTQ files in the `readsdir`.
All fields are passed through to the aggregation steps at the end of the pipeline.

### Output

Output files are written, by default to a `results` directory where the pipeline is called from.
This folder is organised for each step that emits files and results like so:

```none
results
├── aggregated.csv
├── all_lineages.csv
├── <LAB1>
│  ├── <RUN1>
│  │  ├── alignments
│  │  ├── ampliconstats
│  │  ├── bamclipper
│  │  ├── freyja
│  │  ├── mosdepth
│  │  ├── stats_csv
│  │  └── trimmed
│  └── <RUN2>
│     ├── alignments
│     ├── ampliconstats
│     ├── bamclipper
│     ├── freyja
│     ├── mosdepth
│     ├── stats_csv
│     └── trimmed
└── <LAB2>
   └── <RUN1>
      ├── alignments
      ├── ampliconstats
      ├── bamclipper
      ├── freyja
      ├── mosdepth
      ├── stats_csv
      └── trimmed
```

Outputs that are organised in directories under a `<RUN ID>` are the raw outputs from the steps in [pipeline summary](#pipeline-summary).
The aggregated outputs are placed at the top level as these combine data from all of the sequencing labs and runs.

#### _`aggregated.csv`_

This CSV file aggregates the WHO lineages, their frequencies, and sequencing depths for all the samples that are able to complete analysis.
As multiple lineages maybe present multiple rows can be returned for a single sample.

| Column                     | Description                                                       |
| :------------------------- | :---------------------------------------------------------------- |
| amplicon_mean              | Mean coverage over all amplicons including zeros                  |
| non_zero_amplicon_mean     | Mean coverage over amplicons excluding zeros                      |
| amplicon_median            | Median coverage over all amplicons including zeros                |
| non_zero_amplicon_median   | Median coverage over amplicons excluding zeros                    |
| count_gte_20               | Count of amplicons with at least (≥) 20× coverage                 |
| count_lt_20                | Count of amplicons with less than (<) 20× coverage                |
| stdev                      | Standard deviation of coverage over all amplicons                 |
| non_zero_stdev             | Standard deviation of coverage over all amplicons excluding zeros |
| lineage                    | WHO lineage assigned by Freyja                                    |
| abundance                  | Abundance of this WHO lineage                                     |
| mean_genome_coverage       | Mean coverage over whole genome from mosdepth                     |
| sample_id                  | Sample ID used in the pipeline                                    |
| sample_site_code           | Sample site location code                                         |
| timestamp_sample_collected | Timestamp sample collected                                        |
| sequencing_lab_code        | Sequencing lab                                                    |
| original_sample_id         | Original metadata sample id                                       |
| sequencing_sample_id       | Sample ID used in the pipeline                                    |
| sequencing_run_id          | Run ID for the sample                                             |
| amplicon_001_mean_depth    | Coverage over this individual amplicon                            |
| ...                        | ...                                                               |
| amplicon_154_mean_depth    | repeated for all amplicons                                        |


#### _`all_lineages.csv`_

This CSV file aggregates Pango lineages that are assigned by Freyja.
It is a more fine-grained breakdown of the sample composition than the WHO lineages.

| Column                     | Description                      |
| :------------------------- | :------------------------------- |
| lineage                    | Pango lineage assigned by Freyja |
| abundance                  | Abundance of this lineage        |
| sample_id                  | Sample ID used in the pipeline   |
| sample_site_code           | Sample site location code        |
| timestamp_sample_collected | Timestamp sample collected       |
| sequencing_lab_code        | Sequencing lab                   |
| original_sample_id         | Original metadata sample id      |
| sequencing_sample_id       | Sample ID used in the pipeline   |
| sequencing_run_id          | Run ID for the sample            |

 [nimagen-primers]: https://www.nimagen.com/covid19
 [pigz]: https://zlib.net/pigz/
 [fastp]: https://github.com/OpenGene/fastp
 [bwa-mem]: https://github.com/lh3/bwa
 [samtools]: https://www.htslib.org/
 [bamclipper]: https://github.com/tommyau/bamclipper
 [mosdepth]: https://github.com/brentp/mosdepth/
 [ampliconstats]: https://www.htslib.org/doc/samtools-ampliconstats.html
 [freyja]: https://github.com/andersen-lab/freyja
 [conda]: https://docs.conda.io/en/latest/
 [miniconda]: https://docs.conda.io/en/latest/miniconda.html
 [collate_results]: bin/collate_results.py
 [collate_lineages]: bin/collate_lineages.py
 [xsv]: https://github.com/BurntSushi/xsv
 [nf-config]: nextflow.config
 [symlink-help]: scripts
