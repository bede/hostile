from pathlib import Path

import defopt

from hostile import lib


def dehost(*, fastq1: Path, fastq2: Path | None = None, out_dir: Path = lib.CWD) -> None:
    """
    Dehost fastqs using minimap2

    :arg fastq1: path to fastq.gz file
    :arg fastq2: path to optional second fastq.gz file (for paired reads)
    :arg out_dir: output directory for decontaminated fastq.gz files
    :arg debug: show debug messages
    """
    lib.dehost_fastqs(fastq1, fastq2, out_dir=out_dir)


def main():
    defopt.run(dehost, no_negated_flags=True)
