import json
from pathlib import Path

import defopt

from hostile import lib


def clean(
    *,
    fastq1: Path,
    fastq2: Path | None = None,
    aligner: lib.ALIGNERS = lib.ALIGNERS.bowtie2,
    custom_index: Path | None = None,
    out_dir: Path = lib.CWD,
    threads: int = lib.THREADS,
    debug: bool = False,
) -> None:
    """
    Remove human reads from paired fastq(.gz) files

    :arg fastq1: path to forward fastq(.gz) file
    :arg fastq2: optional path to reverse fastq(.gz) file
    :arg aligner: alignment algorithm
    :arg custom_index: path to custom index
    :arg out_dir: output directory for decontaminated fastq.gz files
    :arg threads: number of CPU threads to use
    :arg debug: show debug messages
    """
    if fastq2:
        stats = lib.clean_paired_fastqs(
            [(fastq1, fastq2)],
            out_dir=out_dir,
            threads=threads,
            aligner=aligner,
            custom_index=custom_index,
        )
    else:
        stats = lib.clean_fastqs(
            [fastq1],
            out_dir=out_dir,
            threads=threads,
            aligner=aligner,
            custom_index=custom_index,
        )
    print(json.dumps(stats, indent=4))


def clean_many(
    *fastqs: str,
    aligner: lib.ALIGNERS = lib.ALIGNERS.bowtie2,
    custom_index: Path | None = None,
    out_dir: Path = lib.CWD,
    threads: int = lib.THREADS,
    debug: bool = False,
) -> None:
    """
    Remove human reads from comma-separated pairs of fastq(.gz) files

    :arg fastqs: path to fastq(.gz) or bam file(s). Paired fastq paths should be comma-separated, e.g. reads_1.fastq.gz,reads_2.fastq.gz
    :arg aligner: alignment algorithm
    :arg custom_index: path to custom index
    :arg out_dir: output directory for decontaminated fastq.gz files
    :arg threads: number of threads to use
    :arg debug: show debug messages
    """
    if "," in fastqs[0]:  # Paired fastq
        paired_fastqs = [tuple(pair.split(",")) for pair in fastqs]
        paired_fastqs = [tuple([Path(fq1), Path(fq2)]) for fq1, fq2 in paired_fastqs]
        stats = lib.clean_paired_fastqs(
            paired_fastqs,
            out_dir=out_dir,
            threads=threads,
            aligner=aligner,
            custom_index=custom_index,
        )
        print(json.dumps(stats, indent=4))
    else:
        raise NotImplementedError(
            "Forward and reverse fastq(.gz) paths should be separated with a comma"
        )


def mask(
    reference: Path, target: Path, out_dir: Path = Path("masked"), threads: int = 1
) -> None:
    """
    Mask a reference genome against fasta files in a specified directory

    :arg reference: path to reference genome in fasta[.gz] format
    :arg target: path to target genome(s) in fasta[.gz] format
    :arg out_dir: path of output directory
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
