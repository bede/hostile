[![Tests](https://github.com/bede/hostile/actions/workflows/test.yml/badge.svg)](https://github.com/bede/hostile/actions/workflows/test.yml)

# Hostile

Rapid FASTQ decontamination by host depletion. Accepts paired fastq.gz files as arguments and outputs paired fastq.gz files. Downloads and caches a custom human reference genome to `$XDG_DATA_DIR`. Replaces read headers with incrementing integers for speed and privacy. Python package with CLI and Python API. Installs with conda/mamba.



## Install

```bash
git clone https://github.com/bede/hostile.git
cd hostile
conda env create -f environment.yml  # Use mamba if impatient
conda activate hostile
pip install .
```



## Command line usage

```bash
% hostile clean --help
usage: hostile clean [-h] --fastq1 FASTQ1 --fastq2 FASTQ2 [--aligner {bowtie2,minimap2}] [--out-dir OUT_DIR] [--threads THREADS] [--debug]

Remove human reads from paired fastq.gz files

options:
  -h, --help            show this help message and exit
  --fastq1 FASTQ1       path to forward fastq.gz file
  --fastq2 FASTQ2       path to reverse fastq.gz file
  --aligner {bowtie2,minimap2}
                        alignment algorithm
                        (default: bowtie2)
  --out-dir OUT_DIR     output directory for decontaminated fastq.gz files
                        (default: /Users/bede/Research/Git/hostile)
  --threads THREADS     number of CPU threads to use
                        (default: 10)
  --debug               show debug messages
                        (default: False)

```


```bash
% hostile clean --fastq1 reads.r1.fastq.gz --fastq2 reads.r2.fastq.gz
INFO: Using Bowtie2
INFO: Using cached human index (~/Library/Application Support/hostile/human-bowtie2)
Cleaning: 100%|█████████████████████████████████████████████| 1/1 [00:00<00:00,  2.40it/s]
[
    {
        "fastq1_in_name": "reads.r1.fastq.gz",
        "fastq2_in_name": "reads.r2.fastq.gz",
        "fastq1_in_path": "tests/data/reads.r1.fastq.gz",
        "fastq2_in_path": "tests/data/reads.r2.fastq.gz",
        "fastq1_out_name": "reads.r1.dehosted_1.fastq.gz",
        "fastq2_out_name": "reads.r2.dehosted_2.fastq.gz",
        "fastq1_out_path": "~/Research/Git/hostile/reads.r1.dehosted_1.fastq.gz",
        "fastq2_out_path": "~/Research/Git/hostile/reads.r2.dehosted_2.fastq.gz",
        "reads_in": 20,
        "reads_out": 20,
        "reads_removed": 0,
        "reads_removed_proportion": 0.0
    }
]

```



## Python usage

```python
from hostile.lib import clean_paired_fastqs

decontamination_statistics = clean_paired_fastqs(fastqs=[("h37rv_10.r1.fastq.gz", fastq2="h37rv_10.r1.fastq.gz")])
```



## Development

```bash
git clone https://github.com/bede/hostile.git
cd hostile
conda env create -f environment.yml  # Use mamba if impatient
conda activate hostile
pip install --editable '.[dev]'
pytest
```
