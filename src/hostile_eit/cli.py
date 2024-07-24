import json
import logging

from enum import Enum
from pathlib import Path
from typing import Literal

import defopt

from hostile_eit import lib, util


class ALIGNER(Enum):
    """Provides auto enum for CLI, not to be confused with lib.ALIGNER"""

    bowtie2 = "bowtie2"
    minimap2 = "minimap2"
    auto = "auto"


def clean(
    *,
    fastq1: Path,
    fastq2: Path | None = None,
    aligner: ALIGNER = ALIGNER.auto,
    index: str = util.DEFAULT_INDEX_NAME,
    invert: bool = False,
    rename: bool = False,
    reorder: bool = False,
    out_dir: Path = util.CWD,
    threads: int = util.THREADS,
    aligner_args: str = "",
    force: bool = False,
    offline: bool = False,
    debug: bool = False,
) -> None:
    """
    Remove reads aligning to an index from fastq[.gz] input files

    :arg fastq1: path to forward fastq[.gz] file
    :arg fastq2: optional path to reverse fastq[.gz] file
    :arg aligner: alignment algorithm. Default is Bowtie2 (paired reads) & Minimap2 (unpaired reads)
    :arg index: name of standard index or path to custom genome/index
    :arg invert: keep only reads aligning to the target genome (and their mates if applicable)
    :arg rename: replace read names with incrementing integers
    :arg reorder: ensure deterministic output order
    :arg out_dir: path to output directory
    :arg threads: number of alignment threads. A sensible default is chosen automatically
    :arg aligner_args: additional arguments for alignment
    :arg force: overwrite existing output files
    :arg offline: disable automatic index download
    :arg debug: show debug messages
    """

    if debug:
        logging.getLogger().setLevel(logging.DEBUG)
    aligner_paired = (
        lib.ALIGNER.bowtie2
        if aligner == ALIGNER.auto or aligner == ALIGNER.bowtie2
        else lib.ALIGNER.minimap2
    )
    aligner_unpaired = (
        lib.ALIGNER.minimap2
        if aligner == ALIGNER.auto or aligner == ALIGNER.minimap2
        else lib.ALIGNER.bowtie2
    )
    if fastq2:
        stats = lib.clean_paired_fastqs(
            [(fastq1, fastq2)],
            index=index,
            invert=invert,
            rename=rename,
            reorder=reorder,
            out_dir=out_dir,
            aligner=aligner_paired,
            aligner_args=aligner_args,
            threads=threads,
            force=force,
            offline=offline,
        )
    else:
        stats = lib.clean_fastqs(
            [fastq1],
            index=index,
            invert=invert,
            rename=rename,
            reorder=reorder,
            out_dir=out_dir,
            aligner=aligner_unpaired,
            aligner_args=aligner_args,
            threads=threads,
            force=force,
            offline=offline,
        )
    print(json.dumps(stats, indent=4))


def mask(
    reference: Path,
    target: Path,
    out_dir: Path = Path("masked"),
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
    :arg out_dir: path to output directory
    :arg threads: number of threads to use
    """
    lib.mask(
        reference=reference,
        target=target,
        out_dir=out_dir,
        kmer_length=kmer_length,
        kmer_step=kmer_step,
        threads=threads,
    )


def fetch(
    name: str = util.DEFAULT_INDEX_NAME,
    aligner: Literal["minimap2", "bowtie2", "both"] = "both",
    list: bool = False,
) -> None:
    """
    Download and cache indexes from object storage for use with hostile_eit clean

    :arg name: name of index to download
    :arg aligner: aligner(s) for which to download an index
    :arg list: list available indexes
    """
    logging.info(f"Cache directory: {util.CACHE_DIR}")
    logging.info(f"Manifest URL: {util.BUCKET_URL}/manifest.json")
    if list:
        manifest = util.fetch_manifest()
        for name in manifest.keys():
            print(name)
    else:
        if aligner == "minimap2" or aligner == "both":
            logging.info(f"Looking for Minimap2 index {name}")
            lib.ALIGNER.minimap2.value.check_index(name)
        if aligner == "bowtie2" or aligner == "both":
            logging.info(f"Looking for Bowtie2 index {name}")
            lib.ALIGNER.bowtie2.value.check_index(name)


def main():
    defopt.run(
        {"clean": clean, "mask": mask, "fetch": fetch},
        no_negated_flags=True,
        strict_kwonly=False,
        short={},
    )


# def clean_many(
#     *fastqs: str,
#     aligner: lib.ALIGNER = lib.ALIGNER.bowtie2,
#     index: Path | None = None,
#     out_dir: Path = util.CWD,
#     threads: int = lib.THREADS,
#     debug: bool = False,
# ) -> None:
#     """
#     Remove human reads from comma-separated pairs of fastq(.gz) files

#     :arg fastqs: path to fastq(.gz) or bam file(s). Paired fastq paths should be comma-separated, e.g. reads_1.fastq.gz,reads_2.fastq.gz
#     :arg aligner: alignment algorithm
#     :arg index: path to custom genome or index. For Bowtie2, provide an index path without the .bt2 extension
#     :arg out_dir: path to output directory
#     :arg threads: number of threads to use
#     :arg debug: show debug messages
#     """
#     if "," in fastqs[0]:  # Paired fastq
#         paired_fastqs = [tuple(pair.split(",")) for pair in fastqs]
#         paired_fastqs = [tuple([Path(fq1), Path(fq2)]) for fq1, fq2 in paired_fastqs]
#         stats = lib.clean_paired_fastqs(
#             paired_fastqs,
#             out_dir=out_dir,
#             threads=threads,
#             aligner=aligner,
#             index=index,
#         )
#         print(json.dumps(stats, indent=4))
#     else:
#         raise NotImplementedError(
#             "Forward and reverse fastq(.gz) paths should be separated with a comma"
#         )
