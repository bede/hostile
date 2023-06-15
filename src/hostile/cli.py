import json
from pathlib import Path

import defopt

from hostile import lib


def dehost(
    *reads: str,
    aligner: lib.ALIGNERS = lib.ALIGNERS.bowtie2,
    out_dir: Path = lib.CWD,
    threads: int = lib.THREADS,
    debug: bool = False,
) -> None:
    """
    Dehost fastqs using minimap2

    :arg reads: path to fastq.gz or bam file(s). Paired fastq paths should be comma-separated, e.g. reads_1.fastq.gz,reads_2.fastq.gz
    :arg aligner: alignment algorithm
    :arg out_dir: output directory for decontaminated fastq.gz files
    :arg threads: number of CPU threads to use
    :arg debug: show debug messages
    """
    if "," in reads[0]:  # Paired fastq
        paired_fastqs = [tuple(pair.split(",")) for pair in reads]
        paired_fastqs = [tuple([Path(fq1), Path(fq2)]) for fq1, fq2 in paired_fastqs]
        checksums = lib.dehost_paired_fastqs(
            paired_fastqs, out_dir=out_dir, threads=threads, aligner=aligner
        )
        print(json.dumps(checksums, indent=4))


def main():
    defopt.run(dehost, no_negated_flags=True)
