import json
import logging
import sys

from enum import Enum
from pathlib import Path

import defopt

from hostile import lib, util


class ALIGNER(Enum):
    """Provides auto enum for CLI, not to be confused with lib.ALIGNER"""

    bowtie2 = "bowtie2"
    minimap2 = "minimap2"
    auto = "auto"


def clean(
    *,
    fastq1: str,
    fastq2: str = "",
    aligner: ALIGNER = ALIGNER.auto,
    index: str = util.DEFAULT_INDEX_NAME,
    invert: bool = False,
    rename: bool = False,
    reorder: bool = False,
    casava: bool = False,
    output: Path = util.CWD,
    aligner_args: str = "",
    threads: int = util.CPU_COUNT,
    force: bool = False,
    airplane: bool = False,
    debug: bool = False,
) -> None:
    """
    Remove reads aligning to an index from fastq[.gz] input files or stdin.

    :arg fastq1: path to forward fastq[.gz] file
    :arg fastq2: optional path to reverse fastq[.gz] file (short reads only)
    :arg aligner: alignment algorithm. Defaults to minimap2 (long read) given fastq1 only or bowtie2 (short read)
        given fastq1 and fastq2. Override with bowtie2 for single/unpaired short reads
    :arg index: name of standard index or path to custom genome (Minimap2) or Bowtie2 index
    :arg invert: keep only reads aligning to the index (and their mates if applicable)
    :arg rename: replace read names with incrementing integers
    :arg reorder: ensure deterministic output order
    :arg casava: use Casava 1.8+ read header format
    :arg output: path to output directory or - for stdout
    :arg aligner_args: additional arguments for alignment
    :arg threads: number of alignment threads. A sensible default is chosen automatically
    :arg force: overwrite existing output files
    :arg airplane: disable automatic index download (offline mode)
    :arg debug: show debug messages
    """

    if debug:
        logging.getLogger().setLevel(logging.DEBUG)
    aligner_paired = (
        lib.ALIGNER.bowtie2
        if aligner == ALIGNER.auto or aligner == ALIGNER.bowtie2
        else lib.ALIGNER.minimap2
    )
    aligner_single = (
        lib.ALIGNER.minimap2
        if aligner == ALIGNER.auto or aligner == ALIGNER.minimap2
        else lib.ALIGNER.bowtie2
    )
    if fastq2:
        stats = lib.clean_paired_fastqs(
            [(fastq1, fastq2)],
            aligner=aligner_paired,
            index=index,
            invert=invert,
            rename=rename,
            reorder=reorder,
            casava=casava,
            output=output,
            aligner_args=aligner_args,
            threads=threads,
            force=force,
            airplane=airplane,
        )
    else:
        stats = lib.clean_fastqs(
            [fastq1],
            aligner=aligner_single,
            index=index,
            invert=invert,
            rename=rename,
            reorder=reorder,
            casava=casava,
            output=output,
            aligner_args=aligner_args,
            threads=threads,
            force=force,
            airplane=airplane,
        )
    print(
        json.dumps(stats, indent=4),
        file=sys.stderr if str(output) == "-" else sys.stdout,
    )


def mask(
    reference: Path,
    target: Path,
    output: Path = Path("masked"),
    kmer_length: int = 150,
    kmer_step: int = 10,
    threads: int = util.CPU_COUNT,
) -> None:
    """
    Mask reference genome against target genome(s)

    :arg reference: path to reference genome in fasta(.gz) format
    :arg target: path to target genome(s) in fasta(.gz) format
    :arg kmer_length: length of target genome k-mers
    :arg kmer_step: interval between target genome k-mer start positions
    :arg output: path to output directory
    :arg threads: number of threads to use
    """
    lib.mask(
        reference=reference,
        target=target,
        output=output,
        kmer_length=kmer_length,
        kmer_step=kmer_step,
        threads=threads,
    )


def fetch_index(
    name: str = util.DEFAULT_INDEX_NAME,
    minimap2: bool = False,
    bowtie2: bool = False,
) -> None:
    """
    Download and cache indexes from object storage for use with hostile clean

    :arg name: name of index to download
    :arg minimap2: fetch Minimap2 index
    :arg bowtie2: fetch Bowtie2 index
    """
    lib.fetch_index(name=name, minimap2=minimap2, bowtie2=bowtie2)


def list_indexes(airplane: bool = False):
    """
    List available remote and local cached indexes

    :arg airplane: list only local cached indexes (offline mode)
    """
    lib.list_indexes(airplane=airplane)


def delete_index(name: str = "", all: bool = False, mmi: bool = False) -> None:
    """
    Delete cached indexes

    :arg name: name of cached index to delete
    :arg all: delete all cached indexes
    :arg mmi: delete all cached Minimap2 indexes
    """
    lib.delete_index(name=name, all=all, mmi=mmi)


def main():
    defopt.run(
        {
            "clean": clean,
            "mask": mask,
            "index": {
                "delete": delete_index,
                "list": list_indexes,
                "fetch": fetch_index,
            },
        },
        no_negated_flags=True,
        strict_kwonly=False,
    )


if __name__ == "__main__":
    main()
