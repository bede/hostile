# Hostile

Host decontamination using minimap2 & samtools. Accepts paired fastq.gz files as arguments and outputs paired fastq.gz files. Downloads and caches a custom human reference genome to `$XDG_DATA_DIR`. Replaces read headers with incrementing integers for speed and privacy. Python package with CLI and Python API. Installs with conda/mamba.



## Install

```
conda create -c bioconda -c conda-forge -n hostile python=3 minimap2 samtools seqkit pigz
conda activate hostile
git clone https://github.com/bede/hostile.git
cd hostile
pip install .
```



## Command line usage

```
% hostile -h    
usage: hostile [-h] --fastq1 FASTQ1 [--fastq2 FASTQ2] [-o OUT_DIR] [--version]

Dehost fastqs using minimap2

options:
  -h, --help            show this help message and exit
  --fastq1 FASTQ1       path to fastq.gz file
  --fastq2 FASTQ2       path to optional second fastq.gz file (for paired reads)
                        (default: None)
  -o OUT_DIR, --out-dir OUT_DIR
                        output directory for decontaminated fastq.gz files
                        (default: /Users/bede/Research/Git/hostile)
  --version             show program's version number and exit
```

```bash
hostile --fastq1 reads.r1.fastq.gz --fastq2 reads.r2.fastq.gz
{
    "reads.r1.fastq.gz": "46443dbf65d372ac7a6857d867cbc8b4763aa0c2fce8778fb0e051eda30cc4f6",
    "reads.r1.dehosted_1.fastq.gz": "ead10ee41f4bb0945f1792d963e6f02cacd1d589a8bc1b941fb72a60958eebed",
    "reads.r2.fastq.gz": "a254ef8341493056ab0b08b315bbb1e4d77020b47fbbd658e57991507d3e08a0",
    "reads.r2.dehosted_2.fastq.gz": "43192f7e2227e7e1c6f973b8712f9790612861929219d24ec004678851c96e9c"
}
```



## Python usage

```python
from hostile.lib import dehost_fastqs

checksums = dehost_fastqs(fastq1="h37rv_10.r1.fastq.gz", fastq2="h37rv_10.r1.fastq.gz")
```



## Development

```
conda create -c bioconda -c conda-forge -n hostile python=3 minimap2 samtools seqkit pigz
conda activate hostile
git clone https://github.com/bede/hostile.git
cd hostile
pip install --editable '.[dev]'
pytest
```

