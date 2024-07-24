[![Tests](https://github.com/EIT-Pathogena/hostile-eit/actions/workflows/test.yml/badge.svg)](https://github.com/EIT-Pathogena/hostile-eit
/actions/workflows/test.yml) [![PyPI version](https://img.shields.io/pypi/v/hostile)](https://pypi.org/project/hostile/) [![Bioconda version](https://anaconda.org/bioconda/hostile/badges/version.svg)](https://anaconda.org/bioconda/hostile/) [![Downloads](https://img.shields.io/conda/dn/bioconda/hostile.svg)](https://anaconda.org/bioconda/hostile) [![DOI:10.1101/2023.07.04.547735](https://img.shields.io/badge/citation-10.1093/bioinformatics/btad728-blue)](https://doi.org/10.1093/bioinformatics/btad728)

<p align="center">
    <img src="logo-250px.png" alt="logo">
</p>

# Hostile EIT

This repository is a fork of the original Hostile project, created to ensure the long-term support and maintainability
of our products which make use of it. By forking the original repository, we will manage updates, fix bugs, and
implement features more effectively while staying aligned with the core functionality of the original codebase.
This approach allows for dedicated attention to our product's specific needs, ensuring it remains robust and up-to-date. 
Through active maintenance and regular updates, we aim to enhance the product's longevity and performance, 
committing general fixes and support upstream.

## Hostile

Hostile accurately removes host sequences from short and long read (meta)genomes, consuming paired or
unpaired `fastq[.gz]` input. Batteries are included – a human reference genome is downloaded when run for the first
time. Hostile is precise by default, removing
an [order of magnitude fewer microbial reads](https://log.bede.im/2023/08/29/precise-host-read-removal.html#evaluating-accuracy)
than existing approaches while removing >99.5% of real human reads from 1000 Genomes Project samples. For the best
possible retention of microbial reads, use an existing index masked against bacterial and/or viral genomes, or make your
own using the built-in masking utility. Read headers can be replaced with integers (using `--rename`) for privacy and
smaller FASTQs. Heavy lifting is done with fast existing tools (Minimap2/Bowtie2 and Samtools). Bowtie2 is the default
aligner for short (paired) reads while Minimap2 is default aligner for long reads. In benchmarks, bacterial Illumina
reads were decontaminated at 32Mbp/s (210k reads/sec) and bacterial ONT reads at 22Mbp/s, using 8 alignment threads. By
default, Hostile requires 4GB of RAM for decontaminating short reads and 13GB for long reads (Minimap2). Further
information and benchmarks can be found in the [paper](https://doi.org/10.1093/bioinformatics/btad728)
and [blog post](https://log.bede.im/2023/08/29/precise-host-read-removal.html). Please open an issue to report
problems [or](mailto:b@bede.im) [otherwise](https://twitter.com/beconsta) [reach](https://bsky.app/profile/bedec.bsky.social) [out](https://mstdn.science/@bede)
for help, advice etc.

## Indexes

The default index `human-t2t-hla` comprises [T2T-CHM13v2.0](https://www.ncbi.nlm.nih.gov/assembly/11828891)
and [IPD-IMGT/HLA](https://www.ebi.ac.uk/ipd/imgt/hla/) v3.51, and is downloaded automatically when running Hostile
unless another index is specified. Slightly higher microbial sequence retention is may be possible using masked indexes,
listed below. The index `human-t2t-hla-argos985` is masked
against [985 reference grade bacterial genomes](https://www.ncbi.nlm.nih.gov/bioproject/231221) including common human
pathogens, while `human-t2t-hla.argos-bacteria-985_rs-viral-202401_ml-phage-202401` is further masked comoprehensively
against all known virus and phage genomes. The latter should be used when retention of viral sequences is a priority. To
use a standard index, simply pass its name as the value of the `--index` argument which takes care of downloading and
caching the relevant index. Automatic download can be disabled using the `--offline` flag, and `--index` can accept a
path to a custom reference genome or Bowtie2
index. [Object storage](https://objectstorage.uk-london-1.oraclecloud.com/n/lr3yhdniv6gu/b/human-genome-indices/o) is
provided by the [ModMedMicro research unit](https://www.expmedndm.ox.ac.uk/modernising-medical-microbiology) at the
University of Oxford.

|                                                                                                                    Name                                                                                                                    |                                                                        Composition                                                                         | Date                   | Masked positions       |
|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|:----------------------------------------------------------------------------------------------------------------------------------------------------------:|------------------------|------------------------|
|                                                                                                       `human-t2t-hla` **(default)**                                                                                                        |                [T2T-CHM13v2.0](https://www.ncbi.nlm.nih.gov/assembly/11828891) + [IPD-IMGT/HLA](https://www.ebi.ac.uk/ipd/imgt/hla/) v3.51                 | 2023-07                | 0 (**0%**)             |
|                                                                                                          `human-t2t-hla-argos985`                                                                                                          |                                 `human-t2t-hla` masked with 150mers for [985](https://github.com/EIT-Pathogena/hostile-eit                                 
                                                       /blob/main/paper/supplementary-table-2.tsv) [FDA-ARGOS](https://www.ncbi.nlm.nih.gov/bioproject/231221) **bacterial** genomes                                                        |                                                                          2023-07                                                                           | 317,973 (**0.010%**)   |
|                                                                                              `human-t2t-hla.rs-viral-202401_ml-phage-202401`                                                                                               | `human-t2t-hla` masked with 150mers for 18,719 RefSeq **viral** and 26,928 [Millard Lab **phage**](https://millardlab.org/phage-genomes-jan-2024/) genomes | 2024-01                | 1,172,993 (**0.037%**) |
|                                                                                     `human-t2t-hla.argos-bacteria-985_rs-viral-202401_ml-phage-202401`                                                                                     |                                 `human-t2t-hla` masked with 150mers for [985](https://github.com/EIT-Pathogena/hostile-eit                                 
 /blob/main/paper/supplementary-table-2.tsv) [FDA-ARGOS](https://www.ncbi.nlm.nih.gov/bioproject/231221) **bacterial**, 18,719 RefSeq **viral**, and 26,928 [Millard Lab **phage**](https://millardlab.org/phage-genomes-jan-2024/) genomes |                                                                          2024-01                                                                           | 1,473,260 (**0.046%**) |
|                                                                                                     `human-t2t-hla-argos985-mycob140`                                                                                                      |                                 `human-t2t-hla` masked with 150mers for [985](https://github.com/EIT-Pathogena/hostile-eit                                 

/blob/main/paper/supplementary-table-2.tsv) [FDA-ARGOS](https://www.ncbi.nlm.nih.gov/bioproject/231221) **bacterial
** & [140](https://github.com/EIT-Pathogena/hostile-eit
/blob/main/paper/supplementary-table-2.tsv) **mycobacterial** genomes | 2023-07 | 319,752 (**0.010%**)   |

*Performance of `human-t2t-hla` and `human-t2t-hla-argos985-mycob140` was evaluated in
the [paper](https://doi.org/10.1093/bioinformatics/btad728)*

## Install  [![Install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square&logo=anaconda)](https://biocontainers.pro/tools/hostile) [![Install with Docker](https://img.shields.io/badge/install%20with-docker-important.svg?style=flat-square&logo=docker)](https://biocontainers.pro/tools/hostile)

Installation with conda/mamba or Docker is recommended due to non-Python dependencies (Bowtie2, Minimap2, Samtools and
Bedtools). Hostile is tested with Ubuntu Linux 22.04, MacOS 12, and under WSL for Windows.

**Conda/mamba**

```bash
conda create -y -n hostile_eit -c conda-forge -c bioconda hostile_eit
conda activate hostile_eit
```

**Docker**

```bash
wget https://raw.githubusercontent.com/EIT-Pathogena/hostile-eit
/main/Dockerfile
docker build . --platform linux/amd64
```

A [Biocontainer image](https://biocontainers.pro/tools/hostile) is also available, but beware that this often lags
behind the latest released version

## Index installation (optional)

Hostile automatically downloads and caches the default index `human-t2t-hla` when run for the first time, meaning that
there is no need to download an index in advance. Neverthless:

- To download and cache the default index (`human-t2t-hla`), run `hostile fetch`
- To list available indexes, run `hostile fetch --list`
- To download and cache another standard index, run e.g. `hostile fetch --name human-t2t-hla-argos985`. This will
  download and cache both short read (Bowtie2) and long read (Minimap2) indexes, unless restricted to one or the other
  using e.g. `--aligner minimap2`.
- To use a custom genome (made with `hostile mask` or otherwise), run `hostile clean`
  with  `--index path/to/genome.fa` (for minimap2) or `--index path/to/bowtie2-index-name` (for Bowtie2). Note that
  Minimap2 mode accepts a path to a genome in fasta format, whereas Bowtie2 mode accepts a path to a precomputed index,
  minus the `.x.bt2` suffix. A Bowtie2 index can be built for use with Hostile using
  e.g. `bowtie2-build genome.fa index-name`.

- To change where indexes are stored, set the environment variable `HOSTILE_CACHE_DIR` to a directory of your choice.
  Run `hostile fetch --list` to verify.

## Command line usage

```bash
$ hostile clean -h
usage: hostile clean [-h] --fastq1 FASTQ1 [--fastq2 FASTQ2] [--aligner {bowtie2,minimap2,auto}] [--index INDEX]
                     [--invert] [--rename] [--reorder] [--out-dir OUT_DIR] [--threads THREADS]
                     [--aligner-args ALIGNER_ARGS] [--force] [--offline] [--debug]

Remove reads aligning to an index from fastq[.gz] input files

options:
  -h, --help            show this help message and exit
  --fastq1 FASTQ1       path to forward fastq[.gz] file
  --fastq2 FASTQ2       optional path to reverse fastq[.gz] file
                        (default: None)
  --aligner {bowtie2,minimap2,auto}
                        alignment algorithm. Default is Bowtie2 (paired reads) & Minimap2 (unpaired reads)
                        (default: auto)
  --index INDEX         name of standard index or path to custom genome/index
                        (default: human-t2t-hla)
  --invert              keep only reads aligning to the target genome (and their mates if applicable)
                        (default: False)
  --rename              replace read names with incrementing integers
                        (default: False)
  --reorder             ensure deterministic output order
                        (default: False)
  --out-dir OUT_DIR     path to output directory
                        (default: /Users/bede/Research/Git/hostile)
  --threads THREADS     number of alignment threads. A sensible default is chosen automatically
                        (default: 5)
  --aligner-args ALIGNER_ARGS
                        additional arguments for alignment
                        (default: )
  --force               overwrite existing output files
                        (default: False)
  --offline             disable automatic index download
                        (default: False)
  --debug               show debug messages
                        (default: False)
```

**Short reads, default index**

```bash
$ hostile clean --fastq1 human_1_1.fastq.gz --fastq2 human_1_2.fastq.gz
INFO: Hostile version 1.0.0. Mode: paired short read (Bowtie2)
INFO: Found cached standard index human-t2t-hla
INFO: Cleaning…
INFO: Cleaning complete
[
    {
        "version": "1.0.0",
        "aligner": "bowtie2",
        "index": "human-t2t-hla",
        "options": [],
        "fastq1_in_name": "human_1_1.fastq.gz",
        "fastq1_in_path": "/Users/bede/human_1_1.fastq.gz",
        "fastq1_out_name": "human_1_1.clean_1.fastq.gz",
        "fastq1_out_path": "/Users/bede/human_1_1.clean_1.fastq.gz",
        "reads_in": 2,
        "reads_out": 0,
        "reads_removed": 2,
        "reads_removed_proportion": 1.0,
        "fastq2_in_name": "human_1_2.fastq.gz",
        "fastq2_in_path": "/Users/bede/human_1_2.fastq.gz",
        "fastq2_out_name": "human_1_2.clean_2.fastq.gz",
        "fastq2_out_path": "/Users/bede/human_1_2.clean_2.fastq.gz"
    }
]
```

**Short reads, masked index, save log**

```bash
$ hostile clean --fastq1 human_1_1.fastq.gz --fastq2 human_1_2.fastq.gz --index human-t2t-hla-argos985 > log.json
INFO: Hostile version 1.0.0. Mode: paired short read (Bowtie2)
INFO: Found cached standard index human-t2t-hla
INFO: Cleaning…
INFO: Cleaning complete
```

**Short unpaired reads, save log**

By default, single fastqs are assumed to be long reads. Override this by specifying `--aligner bowtie2` when
decontaminating unpaired short reads.

```bash
$ hostile clean --aligner bowtie2 --fastq1 tests/data/human_1_1.fastq.gz > log.json
INFO: Hostile version 1.0.0. Mode: short read (Bowtie2)
INFO: Found cached standard index human-t2t-hla
INFO: Cleaning…
INFO: Cleaning complete
```

**Long reads**

```bash
$ hostile clean --fastq1 tests/data/tuberculosis_1_1.fastq.gz
INFO: Hostile version 1.0.0. Mode: long read (Minimap2)
INFO: Found cached standard index human-t2t-hla
INFO: Cleaning…
INFO: Cleaning complete
[
    {
        "version": "1.0.0",
        "aligner": "minimap2",
        "index": "human-t2t-hla",
        "options": [],
        "fastq1_in_name": "tuberculosis_1_1.fastq.gz",
        "fastq1_in_path": "/Users/bede/Research/Git/hostile/tests/data/tuberculosis_1_1.fastq.gz",
        "fastq1_out_name": "tuberculosis_1_1.clean.fastq.gz",
        "fastq1_out_path": "/Users/bede/Research/Git/hostile/tuberculosis_1_1.clean.fastq.gz",
        "reads_in": 1,
        "reads_out": 1,
        "reads_removed": 0,
        "reads_removed_proportion": 0.0
    }
]

```

## Python usage

```python
from pathlib import Path
from hostile.lib import clean_fastqs, clean_paired_fastqs

# Long reads, defaults
clean_fastqs(
    fastqs=[Path("reads.fastq.gz")],
)

# Paired short reads, various options, capture log
log = clean_paired_fastqs(
    fastqs=[(Path("reads_1.fastq.gz"), Path("reads_2.fastq.gz"))],
    index="human-t2t-hla-argos985",
    out_dir=Path("decontaminated-reads"),
    rename=True,
    force=True,
    threads=4
)

print(log)
```

## Masking reference genomes

The `mask` subcommand makes it easy to create custom-masked indexes in order to achieve maximum retention of specific
target organisms:

```bash
hostile mask human.fasta lots-of-bacterial-genomes.fasta --threads 8
```

You may wish to use one of the existing [reference genomes](#reference-genomes--indexes) as a starting point. Masking
uses Minimap2 to align 150mers of the supplied target genomes with the reference genome, and bedtools to mask all
aligned regions with N. Both a masked genome (for Minimap2) and a masked Bowtie2 index is created.

## Limitations

- Hostile prioritises retaining microbial sequences above discarding host sequences. If you strive to remove every last
  human sequence, other approaches may serve you better.
- Performance is not always improved by using all available CPU cores. A sensible default is therefore chosen
  automatically at runtime based on the number of available CPU cores.
- Minimap2 has an overhead of 30-90s for human genome indexing prior to starting decontamination. Surprisingly, loading
  a prebuilt index is not significantly faster. I hope to mitigate this in a future release.

## Citation

Bede Constantinides, Martin Hunt, Derrick W Crook, Hostile: accurate decontamination of microbial host sequences,
*Bioinformatics*, 2023; btad728, https://doi.org/10.1093/bioinformatics/btad728

```latex
@article{10.1093/bioinformatics/btad728,
    author = {Constantinides, Bede and Hunt, Martin and Crook, Derrick W},
    title = {Hostile: accurate decontamination of microbial host sequences},
    journal = {Bioinformatics},
    volume = {39},
    number = {12},
    pages = {btad728},
    year = {2023},
    month = {12},
    issn = {1367-4811},
    doi = {10.1093/bioinformatics/btad728},
    url = {https://doi.org/10.1093/bioinformatics/btad728},
    eprint = {https://academic.oup.com/bioinformatics/article-pdf/39/12/btad728/54850422/btad728.pdf},
}
```

## Development install

```bash
git clone https://github.com/EIT-Pathogena/hostile-eit
.git
cd hostile-eit
conda env create -y -f environment.yml
conda activate hostile-eit
pip install --editable '.[dev]'
pytest
```
