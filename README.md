[![DOI:10.1101/2023.07.04.547735](https://img.shields.io/badge/citation-10.1093/bioinformatics/btad728-blue)](https://doi.org/10.1093/bioinformatics/btad728) [![PyPI version](https://img.shields.io/pypi/v/hostile)](https://pypi.org/project/hostile/) [![Bioconda version](https://anaconda.org/bioconda/hostile/badges/version.svg)](https://anaconda.org/bioconda/hostile/) [![Downloads](https://img.shields.io/conda/dn/bioconda/hostile.svg)](https://anaconda.org/bioconda/hostile) [![Tests](https://github.com/bede/hostile/actions/workflows/test.yml/badge.svg)](https://github.com/bede/hostile/actions/workflows/test.yml)

<p align="center">
    <img width="250" src="logo.png">
</p>

# Hostile

Hostile accurately removes host sequences from short and long read (meta)genomes, consuming single-read or paired `fastq[.gz]` input. Batteries are included – a human reference genome is downloaded when run for the first time. Hostile is precise by default, removing an [order of magnitude fewer microbial reads](https://log.bede.im/2023/08/29/precise-host-read-removal.html#evaluating-accuracy) than existing approaches while removing >99.5% of real human reads from 1000 Genomes Project samples. For the best possible retention of microbial reads, use an existing index masked against bacterial and/or viral genomes, or make your own using the built-in masking utility. Read headers can be replaced with integers (using `--rename`) for privacy and smaller FASTQs. Heavy lifting is done with fast existing tools (Minimap2/Bowtie2 and Samtools). In benchmarks, bacterial Illumina reads were decontaminated at 32Mbp/s (210k reads/sec) and bacterial ONT reads at 22Mbp/s, using 8 alignment threads. In typical use, Hostile requires 4GB of RAM for decontaminating short reads (Bowtie2) and 13GB for long reads (Minimap2). Further information and benchmarks can be found in the [paper](https://doi.org/10.1093/bioinformatics/btad728) and [blog post](https://log.bede.im/2023/08/29/precise-host-read-removal.html). Please open an issue to report problems or otherwise [reach](https://bsky.app/profile/bedec.bsky.social) [out](mailto:b@bede.im) for help and advice, and please cite the paper if you use Hostile in your work.



## Indexes

The default index `human-t2t-hla` comprises [T2T-CHM13v2.0](https://www.ncbi.nlm.nih.gov/assembly/11828891) and [IPD-IMGT/HLA](https://www.ebi.ac.uk/ipd/imgt/hla/) v3.51, and is downloaded automatically when running Hostile unless another index is specified. Higher microbial sequence retention may be possible using masked indexes, which are very easy to use. The index `human-t2t-hla-argos985` is masked against [985 reference grade bacterial genomes](https://www.ncbi.nlm.nih.gov/bioproject/231221) including common human pathogens, while `human-t2t-hla.argos-bacteria-985_rs-viral-202401_ml-phage-202401` is further masked against all known virus and phage genomes. The latter should be used when retention of viral sequences is a priority. To use a standard index, simply pass its name as the value of the `--index` argument, which takes care of downloading and caching the relevant index. [Object storage](https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o) is provided by the [ModMedMicro research unit](https://www.expmedndm.ox.ac.uk/modernising-medical-microbiology) at the University of Oxford. Custom indexes are also supported (see below).

|                             Name                             |                         Composition                          | Date    | Masked positions       |
| :----------------------------------------------------------: | :----------------------------------------------------------: | ------- | ---------------------- |
|                `human-t2t-hla` **(default)**                 | [T2T-CHM13v2.0](https://www.ncbi.nlm.nih.gov/assembly/11828891) + [IPD-IMGT/HLA](https://www.ebi.ac.uk/ipd/imgt/hla/) v3.51 | 2023-07 | 0 (**0%**)             |
|                   `human-t2t-hla-argos985`                   | `human-t2t-hla` masked with 150mers for [985](https://github.com/bede/hostile/blob/main/paper/supplementary-table-2.tsv) [FDA-ARGOS](https://www.ncbi.nlm.nih.gov/bioproject/231221) **bacterial** genomes | 2023-07 | 317,973 (**0.010%**)   |
|       `human-t2t-hla.rs-viral-202401_ml-phage-202401`        | `human-t2t-hla` masked with 150mers for 18,719 RefSeq **viral** and 26,928 [Millard Lab **phage**](https://millardlab.org/phage-genomes-jan-2024/) genomes | 2024-01 | 1,172,993 (**0.037%**) |
| `human-t2t-hla.argos-bacteria-985_rs-viral-202401_ml-phage-202401` | `human-t2t-hla` masked with 150mers for [985](https://github.com/bede/hostile/blob/main/paper/supplementary-table-2.tsv) [FDA-ARGOS](https://www.ncbi.nlm.nih.gov/bioproject/231221) **bacterial**, 18,719 RefSeq **viral**, and 26,928 [Millard Lab **phage**](https://millardlab.org/phage-genomes-jan-2024/) genomes | 2024-01 | 1,473,260 (**0.046%**) |
|              `human-t2t-hla-argos985-mycob140`               | `human-t2t-hla` masked with 150mers for [985](https://github.com/bede/hostile/blob/main/paper/supplementary-table-2.tsv) [FDA-ARGOS](https://www.ncbi.nlm.nih.gov/bioproject/231221) **bacterial** & [140](https://github.com/bede/hostile/blob/main/paper/supplementary-table-2.tsv) **mycobacterial** genomes | 2023-07 | 319,752 (**0.010%**)   |
|                         `mouse-mm39`                         | `GRCm39` ([`GCF_000001635.27`](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001635.27)) | 2024-11 | 0 (**0%**)             |

*Performance of `human-t2t-hla` and `human-t2t-hla-argos985-mycob140` was evaluated in the [paper](https://doi.org/10.1093/bioinformatics/btad728)*



## Install  [![Install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square&logo=anaconda)](https://biocontainers.pro/tools/hostile) [![Install with Docker](https://img.shields.io/badge/install%20with-docker-important.svg?style=flat-square&logo=docker)](https://biocontainers.pro/tools/hostile)

Installation with conda/mamba or Docker is recommended due to non-Python dependencies (Bowtie2, Minimap2, Samtools and Bedtools). Hostile is tested with Ubuntu Linux 22.04, MacOS 12, and under WSL for Windows.

**Conda/mamba**

```bash
conda create -y -n hostile -c conda-forge -c bioconda hostile
conda activate hostile
```

**Docker**

```bash
git clone https://github.com/bede/hostile.git
cd hostile
docker build . --platform linux/amd64
```

A [Biocontainer image](https://biocontainers.pro/tools/hostile) is also available, but beware that this often lags behind the latest released version



## Getting started

```bash
# Long reads
hostile clean --fastq1 long.fastq.gz  # Creates long.clean.fastq.gz
hostile clean --fastq1 --index mouse-mm39  # Use mouse index
cat reads.fastq | hostile clean --fastq1 -  # Read from stdin
hostile clean --fastq1 long.fastq.gz -o - > long.clean.fastq  # Write to stdout
hostile clean --fastq1 long.fastq.gz --invert  # Keep only host reads

# Short reads
hostile clean --fastq1 short.r1.fq.gz --aligner bowtie2  # Single/unpaired
hostile clean --fastq1 short.r1.fq.gz --fastq2 short.r2.fq.gz  # Paired
cat interleaved.fastq | hostile clean --fastq1 - --fastq2 -  # Read interleaved reads from stdin
hostile clean --fastq1 short.r1.fq.gz --fastq2 short.r2.fq.gz -o - > clean.fq  # Write interleaved reads to stdout
```



## Custom indexes

- To list available standard indexes, run `hostile index list`.
- To optionally download and cache the default index (`human-t2t-hla`) ahead of time, run `hostile index fetch`. Include `--minimap2` or `--bowtie2` to download only the respective long or short read index rather than both. To download and cache another standard index, provide its name with e.g. `hostile index fetch --name human-t2t-hla-argos985 --minimap2`.
- To use a custom genome/index (made with `hostile mask` or otherwise), run `hostile clean` with  `--index path/to/genome.fa` (for minimap2) or `--index path/to/bowtie2-index-name` (for Bowtie2). Note that Minimap2 mode accepts a path to a genome in fasta format or .mmi, whereas Bowtie2 mode accepts a path to a precomputed index, minus the `.x.bt2` suffix. A Bowtie2 index can be built for use with Hostile using e.g. `bowtie2-build genome.fa index-name`.
- To change where indexes are stored, set the environment variable `HOSTILE_CACHE_DIR` to a directory of your choice. Run `hostile index list` to verify.
- If you wish to use your own remote repository of indexes, set the environment variable `HOSTILE_REPOSITORY_URL`. Hostile will then look for indexes inside `{HOSTILE_REPOSITORY_URL}/manifest.json`.
- From version 2.0.0 onwards, Hostile automatically builds and reuses .mmi files to speed up long read decontamination with Minimap2. If building an MMI is interrupted, you may receive an error about index corruption. If this happens, run `hostile index delete --mmi`, or if using a custom index, delete the .mmi created in the same directory.



## Command line usage

```bash
$ hostile clean -h
usage: hostile clean [-h] --fastq1 FASTQ1 [--fastq2 FASTQ2] [--aligner {bowtie2,minimap2,auto}] [--index INDEX] [--invert] [--rename] [--reorder] [-c] [-o OUTPUT]
                     [--aligner-args ALIGNER_ARGS] [-t THREADS] [--force] [--airplane] [-d]

Remove reads aligning to an index from fastq[.gz] input files or stdin.

options:
  -h, --help            show this help message and exit
  --fastq1 FASTQ1       path to forward fastq[.gz] file
  --fastq2 FASTQ2       optional path to reverse fastq[.gz] file (short reads only)
                        (default: )
  --aligner {bowtie2,minimap2,auto}
                        alignment algorithm. Defaults to minimap2 (long read) given fastq1 only or bowtie2 (short read)
                        given fastq1 and fastq2. Override with bowtie2 for single/unpaired short reads
                        (default: auto)
  --index INDEX         name of standard index or path to custom genome (Minimap2) or Bowtie2 index
                        (default: human-t2t-hla)
  --invert              keep only reads aligning to the index (and their mates if applicable)
                        (default: False)
  --rename              replace read names with incrementing integers
                        (default: False)
  --reorder             ensure deterministic output order
                        (default: False)
  -c, --casava          use Casava 1.8+ read header format
                        (default: False)
  -o, --output OUTPUT   path to output directory or - for stdout
                        (default: /Users/bede/Research/git/hostile)
  --aligner-args ALIGNER_ARGS
                        additional arguments for alignment
                        (default: )
  -t, --threads THREADS
                        number of alignment threads. A sensible default is chosen automatically
                        (default: 10)
  --force               overwrite existing output files
                        (default: False)
  --airplane            disable automatic index download (offline mode)
                        (default: False)
  -d, --debug           show debug messages
                        (default: False)
```



**Long reads**

Writes compressed fastq.gz files to working directory, sends log to stdout
```bash
$ hostile clean --fastq1 tests/data/tuberculosis_1_1.fastq.gz
INFO: Hostile v2.0.0. Mode: long read (Minimap2)
INFO: Found cached standard index human-t2t-hla (MMI available)
INFO: Cleaning…
INFO: Cleaning complete
[
    {
        "version": "2.0.0",
        "aligner": "minimap2",
        "index": "human-t2t-hla",
        "options": [],
        "fastq1_in_name": "tuberculosis_1_1.fastq.gz",
        "fastq1_in_path": "/Users/bede/Research/git/hostile/tests/data/tuberculosis_1_1.fastq.gz",
        "reads_in": 1,
        "reads_out": 1,
        "reads_removed": 0,
        "reads_removed_proportion": 0.0,
        "fastq1_out_name": "tuberculosis_1_1.clean.fastq.gz",
        "fastq1_out_path": "/Users/bede/Research/git/hostile/tuberculosis_1_1.clean.fastq.gz"
    }
]
```

**Long reads (non-default index, save log)**

```bash
$ hostile clean --fastq1 tests/data/tuberculosis_1_1.fastq.gz --index human-t2t-hla-argos985-mycob140 > log.json
INFO: Hostile v2.0.0. Mode: long read (Minimap2)
INFO: Found cached standard index human-t2t-hla (MMI available)
INFO: Cleaning…
INFO: Cleaning complete
```

**Long reads (`--stdout`)**

Reads sent to stdout, log sent to stderr

```bash
$ hostile clean --fastq1 tests/data/tuberculosis_1_1.fastq.gz -o - > out.fastq
INFO: Hostile v2.0.0. Mode: long read (Minimap2)
INFO: Found cached standard index human-t2t-hla (MMI available)
INFO: Cleaning…
INFO: Cleaning complete
[
    {
        "version": "2.0.0",
        "aligner": "minimap2",
        "index": "human-t2t-hla",
        "options": [
            "stdout"
        ],
        "fastq1_in_name": "tuberculosis_1_1.fastq.gz",
        "fastq1_in_path": "/Users/bede/Research/git/hostile/tests/data/tuberculosis_1_1.fastq.gz",
        "reads_in": 1,
        "reads_out": 1,
        "reads_removed": 0,
        "reads_removed_proportion": 0.0
    }
]
```

**Short paired reads**

When providing both `--fastq1` and `--fastq2`, Hostile asssumes you are providing short reads and uses Bowtie2 automatically.

```bash
$ hostile clean --fastq1 human_1_1.fastq.gz --fastq2 human_1_2.fastq.gz
INFO: Hostile v2.0.0. Mode: paired short read (Bowtie2)
INFO: Found cached standard index human-t2t-hla
INFO: Cleaning…
INFO: Cleaning complete
[
    {
        "version": "2.0.0",
        "aligner": "bowtie2",
        "index": "human-t2t-hla",
        "options": [],
        "fastq1_in_name": "human_1_1.fastq.gz",
        "fastq1_in_path": "/Users/bede/Research/git/hostile/tests/data/human_1_1.fastq.gz",
        "reads_in": 2,
        "reads_out": 0,
        "reads_removed": 2,
        "reads_removed_proportion": 1.0,
        "fastq2_in_name": "human_1_2.fastq.gz",
        "fastq2_in_path": "/Users/bede/Research/git/hostile/tests/data/human_1_2.fastq.gz",
        "fastq1_out_name": "human_1_1.clean_1.fastq.gz",
        "fastq1_out_path": "/Users/bede/Research/git/hostile/human_1_1.clean_1.fastq.gz",
        "fastq2_out_name": "human_1_2.clean_2.fastq.gz",
        "fastq2_out_path": "/Users/bede/Research/git/hostile/human_1_2.clean_2.fastq.gz"
    }
]
```

**Short single/unpaired reads (save log)**

When decontaminating single/unpaired short reads, you must specify `--aligner bowtie2` to override the default long read setting for single/unpaired input. Interleaved input is not supported.

```bash
$ hostile clean --fastq1 human_1_1.fastq.gz --aligner bowtie2 > log.json
INFO: Hostile v2.0.0. Mode: paired short read (Bowtie2)
INFO: Found cached standard index human-t2t-hla-argos985
INFO: Cleaning…
INFO: Cleaning complete
```

**Short paired reads (`--stdout`)**

When using stdout mode with paired input, Hostile sends interleaved paired reads to stdout.

```bash
$ hostile clean --fastq1 human_1_1.fastq.gz --fastq2 human_1_2.fastq.gz -o - > interleaved.fastq
INFO: Hostile v2.0.0. Mode: paired short read (Bowtie2)
INFO: Found cached standard index human-t2t-hla
INFO: Cleaning…
INFO: Cleaning complete
[
    {
        "version": "2.0.0",
        "aligner": "bowtie2",
        "index": "human-t2t-hla",
        "options": [],
        "fastq1_in_name": "human_1_1.fastq.gz",
        "fastq1_in_path": "/Users/bede/Research/git/hostile/tests/data/human_1_1.fastq.gz",
        "reads_in": 2,
        "reads_out": 0,
        "reads_removed": 2,
        "reads_removed_proportion": 1.0,
        "fastq2_in_name": "human_1_2.fastq.gz",
        "fastq2_in_path": "/Users/bede/Research/git/hostile/tests/data/human_1_2.fastq.gz",
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
    threads=12
)

print(log)
```



## Masking reference genomes

The `mask` subcommand makes it easy to create custom-masked indexes in order to achieve maximum retention of specific target organisms:
```bash
hostile mask human.fasta lots-of-bacterial-genomes.fasta --threads 8
```
You may wish to use one of the existing [reference genomes](#reference-genomes--indexes) as a starting point. Masking uses Minimap2 to align 150mers of the supplied target genomes with the reference genome, and bedtools to mask all aligned regions with N. Both a masked genome (for Minimap2) and a masked Bowtie2 index is created.



## Limitations

- Hostile prioritises retaining microbial sequences above discarding host sequences. If you strive to remove every last human sequence, other approaches may serve you better ([blog post](https://log.bede.im/2023/08/29/precise-host-read-removal.html)).
- Performance is not always improved by using all available CPU cores. A sensible default is therefore chosen automatically at runtime based on the number of available CPU cores. For maximum performance you may wish to use `--stdout` mode and compress the fastq stream with zstandard, a faster gzip alternative.



## Citation

Please cite Hostile if you find it useful.

Bede Constantinides, Martin Hunt, Derrick W Crook,  Hostile: accurate decontamination of microbial host sequences, *Bioinformatics*, 2023; btad728, https://doi.org/10.1093/bioinformatics/btad728

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
git clone https://github.com/bede/hostile.git
cd hostile
conda env create -y -f environment.yml
conda activate hostile
pip install --editable '.[dev]'
pytest
pre-commit install
```
