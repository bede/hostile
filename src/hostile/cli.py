import json
import logging
from enum import Enum
from pathlib import Path

import defopt

from hostile import lib


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
    index: Path | None = None,
    rename: bool = False,
    out_dir: Path = lib.CWD,
    threads: int = lib.THREADS,
    force: bool = False,
    debug: bool = False,
) -> None:
    """
    Remove host reads from paired fastq(.gz) files

    :arg fastq1: path to forward fastq(.gz) file
    :arg fastq2: optional path to reverse fastq(.gz) file
    :arg aligner: alignment algorithm
    :arg index: path to custom genome or index. For Bowtie2, provide an index path without the .bt2 extension
    :arg rename: replace read names with incrementing integers
    :arg out_dir: path to output directory
    :arg threads: number of CPU threads to use
    :arg force: overwrite existing output files
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
            rename=rename,
            out_dir=out_dir,
            aligner=aligner_paired,
            threads=threads,
            force=force,
        )
    else:
        stats = lib.clean_fastqs(
            [fastq1],
            index=index,
            rename=rename,
            out_dir=out_dir,
            aligner=aligner_unpaired,
            threads=threads,
            force=force,
        )
    print(json.dumps(stats, indent=4))


def mask(
    reference: Path, target: Path, out_dir: Path = Path("masked"), threads: int = 1
) -> None:
    """
    Mask reference genome against target genome[s]

    :arg reference: path to reference genome in fasta[.gz] format
    :arg target: path to target genome(s) in fasta[.gz] format
    :arg out_dir: path to output directory
    :arg threads: number of threads to use
    """
    lib.mask(reference=reference, target=target, out_dir=out_dir, threads=threads)


def main():
    defopt.run(
        {"clean": clean, "mask": mask},
        no_negated_flags=True,
        strict_kwonly=False,
        short={},
    )


# def clean_many(
#     *fastqs: str,
#     aligner: lib.ALIGNER = lib.ALIGNER.bowtie2,
#     index: Path | None = None,
#     out_dir: Path = lib.CWD,
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
