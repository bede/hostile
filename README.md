# Hostile

Host decontamination using minimap2.



## Install

```
git clone https://github.com/bede/hostile.git
cd hostile
pip install .
```



## Usage

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



```
hostile --fastq1 reads_1.fastq.gz --fastq2 reads_2.fastq.gz
# generates reads_1.dehosted_1.fastq.gz reads_2.dehosted_1.fastq.gz
```



## Development

```
git clone https://github.com/bede/hostile.git
cd hostile
pip install --editable '.[dev]'
pytest
```

