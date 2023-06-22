import json
from pathlib import Path

import defopt

from hostile import lib


def clean(
    *,
    fastq1: Path,
    fastq2: Path,
    aligner: lib.ALIGNERS = lib.ALIGNERS.bowtie2,
    out_dir: Path = lib.CWD,
    threads: int = lib.THREADS,
    debug: bool = False,
) -> None:
    """
    Remove human reads from paired fastq(.gz) files

    :arg fastq1: path to forward fastq(.gz) file
    :arg fastq2: path to reverse fastq(.gz) file
    :arg aligner: alignment algorithm
    :arg out_dir: output directory for decontaminated fastq.gz files
    :arg threads: number of CPU threads to use
    :arg debug: show debug messages
    """
    stats = lib.clean_paired_fastqs(
        [(fastq1, fastq2)], out_dir=out_dir, threads=threads, aligner=aligner
    )
    print(json.dumps(stats, indent=4))


def clean_many(
    *reads: str,
    aligner: lib.ALIGNERS = lib.ALIGNERS.bowtie2,
    out_dir: Path = lib.CWD,
    threads: int = lib.THREADS,
    debug: bool = False,
) -> None:
    """
    Remove human reads from comma-separated pairs of fastq(.gz) files

    :arg reads: path to fastq(.gz) or bam file(s). Paired fastq paths should be comma-separated, e.g. reads_1.fastq.gz,reads_2.fastq.gz
    :arg aligner: alignment algorithm
    :arg out_dir: output directory for decontaminated fastq.gz files
    :arg threads: number of threads to use
    :arg debug: show debug messages
    """
    if "," in reads[0]:  # Paired fastq
        paired_fastqs = [tuple(pair.split(",")) for pair in reads]
        paired_fastqs = [tuple([Path(fq1), Path(fq2)]) for fq1, fq2 in paired_fastqs]
        stats = lib.clean_paired_fastqs(
            paired_fastqs, out_dir=out_dir, threads=threads, aligner=aligner
        )
        print(json.dumps(stats, indent=4))
    else:
        raise NotImplementedError(
            "Forward and reverse fastq(.gz) paths should be separated with a comma"
        )


def main():
    defopt.run(
        {
            "clean": clean,
            "clean-many": clean_many,
        },
        no_negated_flags=True,
        strict_kwonly=False,
        short={},
    )
