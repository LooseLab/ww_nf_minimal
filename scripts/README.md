Helper scripts
==============

Scripts here are to help with creating input for the nextflow pipeline, they are not directly run during the workflow.


_`symlink_exeter.py`_

Input FASTQ files are expected to be in the folder structure `{lab}/{run_id}/{sample_id}*.fastq.gz` but samples from exeter are not in this format.
To get around this constraint, [symlinks](https://en.wikipedia.org/wiki/Symbolic_link) to the input are created using `symlink_exeter.py` script.
This should be run with the `ww_minimal` [conda](https://conda.io) activated.

```bash
python3 symlink_exeter.py <INPUT SAMPLE CSV> <FASTQ SOURCE DIRECTORY> <DESTINATION DIRECTORY>
```
