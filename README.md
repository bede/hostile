[![Tests](https://github.com/bede/hostile/actions/workflows/test.yml/badge.svg)](https://github.com/bede/hostile/actions/workflows/test.yml)

# Hostile

Hostile removes host sequences from short and long reads, consuming paired or unpaired `fastq[.gz]` input and producing `fastq.gz` output. Batteries are included – Hostile downloads and saves a human T2T-CHM13v2.0 + HLA reference genome to `$XDG_DATA_DIR` when run for the first time. Read headers are replaced with incrementing integers for privacy and more compressible FASTQs. Hostile is implemented as a Python package with a CLI and Python API, but all of the heavy lifting is done by fast compiled code (Minimap2/Bowtie2 and Samtools). When used with a masked reference genome, Hostile achieves near-perfect retention of microbial reads while removing >99.5% of human reads. Please read the [BioRxiv preprint](https://www.biorxiv.org/content/10.1101/2023.07.04.547735) for further information and open a GitHub issue, [tweet](https://twitter.com/beconsta) or [toot](https://mstdn.science/@bede) me to report bugs or suggest improvements.



## Reference genomes

The default `human-t2t-hla` reference is downloaded when running Hostile for the first time. This can be overriden by specifying a custom `--index`. Bowtie2 indexes need to be untarred before use. The databases `human-t2t-hla` and `human-t2t-hla-argos985-mycob140`  were compared in the [paper](https://www.biorxiv.org/content/10.1101/2023.07.04.547735).

|               Name                |                         Composition                          |                      Genome (Minimap2)                       |                        Bowtie2 index                         |
| :-------------------------------: | :----------------------------------------------------------: | :----------------------------------------------------------: | :----------------------------------------------------------: |
|   `human-t2t-hla` **(default)**   | [T2T-CHM13v2.0](https://www.ncbi.nlm.nih.gov/assembly/11828891) + [IPD-IMGT/HLA](https://www.ebi.ac.uk/ipd/imgt/hla/) v3.51 | [human-t2t-hla.fa.gz](https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o/human-t2t-hla.fa.gz) | [human-t2t-hla.tar](https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o/human-t2t-hla.tar) |
|     `human-t2t-hla-argos985`      | [T2T-CHM13v2.0](https://www.ncbi.nlm.nih.gov/assembly/11828891) & [IPD-IMGT/HLA](https://www.ebi.ac.uk/ipd/imgt/hla/) v3.51; masked with [985](https://github.com/bede/hostile/blob/main/paper/supplementary-table-2.tsv) [FDA-ARGOS](https://www.ncbi.nlm.nih.gov/bioproject/231221) 150mers | [human-t2t-hla-argos985.fa.gz](https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o/human-t2t-hla-argos985.fa.gz) | [human-t2t-hla-argos985.tar](https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o/human-t2t-hla-argos985.tar) |
| `human-t2t-hla-argos985-mycob140` | [T2T-CHM13v2.0](https://www.ncbi.nlm.nih.gov/assembly/11828891) & [IPD-IMGT/HLA](https://www.ebi.ac.uk/ipd/imgt/hla/) v3.51; masked with [985](https://github.com/bede/hostile/blob/main/paper/supplementary-table-2.tsv) [FDA-ARGOS](https://www.ncbi.nlm.nih.gov/bioproject/231221) & [140](https://github.com/bede/hostile/blob/main/paper/supplementary-table-2.tsv) mycobacterial 150mers | [human-t2t-hla-argos985-mycob140.fa.gz](https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o/human-t2t-hla-argos985-mycob140.fa.gz) | [human-t2t-hla-argos985-mycob140.tar](https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o/human-t2t-hla-argos985-mycob140.tar) |



## Install

Hostile is tested with Ubuntu Linux 22.04 and MacOS 12, and WSL2. Since it has non-Python depedencies (Minimap2, Bowtie2, Samtools & Bedtools), [Mini]conda and Docker are the recommended ways to use it.

### Conda

```bash
curl -OJ https://raw.githubusercontent.com/bede/hostile/main/environment.yml
conda env create -f environment.yml  # Use Mamba if impatient
conda activate hostile
pip install hostile
```



### Docker

*Coming soon*



### Development install

```bash
git clone https://github.com/bede/hostile.git
cd hostile
conda env create -f environment.yml  # Use Mamba if impatient
conda activate hostile
pip install --editable '.[dev]'
pytest
```




## Command line usage

```bash
hostile clean --help
usage: hostile clean [-h] --fastq1 FASTQ1 [--fastq2 FASTQ2] [--aligner {bowtie2,minimap2}] [--index INDEX] [--out-dir OUT_DIR] [--threads THREADS] [--debug]

Remove human reads from paired fastq(.gz) files

options:
  -h, --help            show this help message and exit
  --fastq1 FASTQ1       path to forward fastq(.gz) file
  --fastq2 FASTQ2       optional path to reverse fastq(.gz) file
                        (default: None)
  --aligner {bowtie2,minimap2}
                        alignment algorithm
                        (default: bowtie2)
  --index INDEX         path to custom genome/index. Bowtie2 requires an index given without the .bt2 extension
                        (default: None)
  --out-dir OUT_DIR     output directory for decontaminated fastq.gz files
                        (default: /Users/bede/Research/Git/hostile)
  --threads THREADS     number of CPU threads to use
                        (default: 10)
  --debug               show debug messages
                        (default: False)
```


```bash
% hostile clean --fastq1 reads.r1.fastq.gz --fastq2 reads.r2.fastq.gz
INFO: Paired read input
INFO: Found Bowtie2
INFO: Found cached index (/Users/bede/Library/Application Support/hostile/human-t2t-hla)
INFO: Cleaning…
[
    {
        "fastq1_in_name": "reads.r1.fastq.gz",
        "fastq2_in_name": "reads.r2.fastq.gz",
        "fastq1_in_path": "/path/to/hostile/reads.r1.fastq.gz",
        "fastq2_in_path": "/path/to/hostile/reads.r2.fastq.gz",
        "fastq1_out_name": "reads.r1.dehosted_1.fastq.gz",
        "fastq2_out_name": "reads.r2.dehosted_2.fastq.gz",
        "fastq1_out_path": "/path/to/hostile/reads.r1.dehosted_1.fastq.gz",
        "fastq2_out_path": "/path/to/hostile/reads.r2.dehosted_2.fastq.gz",
        "reads_in": 20,
        "reads_out": 20,
        "reads_removed": 0,
        "reads_removed_proportion": 0.0
    }
]
```



## Python usage

```python
from pathlib import Path
from hostile.lib import clean_paired_fastqs

stats = clean_paired_fastqs(
    fastqs=[(Path("h37rv_10.r1.fastq.gz"), Path("h37rv_10.r1.fastq.gz"))]
)

print(stats)
```



## Masking reference genomes

The `mask` subcommand makes it easy to create custom-masked reference genomes and achieve maximum retention of specific target organisms:
```bash
hostile mask human.fasta lots-of-bacterial-genomes.fasta --threads 8
```
You may wish to use one of the existing [reference genomes](#reference-genomes) as a starting point. Masking uses Minimap2's `asm10` preset to align the supplied target genomes with the reference genome, and bedtools to mask out all aligned regions. This feature requires a [development install](#development-install) until release in version 0.0.3. For Bowtie2—the default aligner for decontaminating short reads—you will also need to build an index before you can use your masked genome with Hostile.
```bash
bowtie2-build masked.fasta masked-index
hostile clean --index masked-index --fastq1 reads_1.fastq.gz --fastq2 reads_2.fastq.gz
```
