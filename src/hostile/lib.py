import logging
import json
import multiprocessing

from enum import Enum
from dataclasses import dataclass
from pathlib import Path

from platformdirs import user_data_dir

from hostile import util
from hostile.aligner import Aligner


logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)


CWD = Path.cwd()
XDG_DATA_DIR = Path(user_data_dir("hostile", "Bede Constantinides"))
THREADS = multiprocessing.cpu_count()


ALIGNERS = Enum(
    "Aligner",
    {
        "bowtie2": Aligner(
            name="Bowtie2",
            short_name="bt2",
            bin_path=Path("bowtie2"),
            # bin_path=Path("/Users/bede/Downloads/bowtie2-2.5.1-macos-arm64/bowtie2"),
            cdn_base_url=f"http://178.79.139.243/hostile",
            working_dir=XDG_DATA_DIR,
            cmd=(
                "{BIN_PATH} -x '{INDEX_PATH}' -1 '{FASTQ1}' -2 '{FASTQ2}'"
                " -k 1 --mm -p {THREADS}"
            ),
            idx_archive_fn="human-bowtie2.tar",
            idx_name="human-bowtie2",
            idx_paths=(
                XDG_DATA_DIR / "human-bowtie2.1.bt2",
                XDG_DATA_DIR / "human-bowtie2.2.bt2",
                XDG_DATA_DIR / "human-bowtie2.3.bt2",
                XDG_DATA_DIR / "human-bowtie2.4.bt2",
                XDG_DATA_DIR / "human-bowtie2.rev.1.bt2",
                XDG_DATA_DIR / "human-bowtie2.rev.2.bt2",
            ),
        ),
        "minimap2": Aligner(
            name="Minimap2",
            short_name="mm2",
            bin_path=Path("minimap2"),
            cdn_base_url=f"http://178.79.139.243/hostile",
            working_dir=XDG_DATA_DIR,
            cmd="{BIN_PATH} -ax sr -m 40 -t {THREADS} '{REF_ARCHIVE_PATH}' '{FASTQ1}' '{FASTQ2}'",
            ref_archive_fn="human.fa.gz",
            idx_name="human.fa.gz",
        ),
    },
)


@dataclass
class SampleReport:
    reads_in: int
    reads_out: int
    delta: int
    proportion_removed: float

    def to_json(self):
        return json.dumps(self.__dict__)


def gather_stats(
    fastqs: list[tuple[Path, Path]], out_dir: Path
) -> dict[str, SampleReport]:
    stats = {}
    for fastq1, fastq2 in fastqs:
        fastq1_stem = util.fastq_path_to_stem(fastq1)
        n_reads_in_path = out_dir / (fastq1_stem + ".reads_in.txt")
        n_reads_out_path = out_dir / (fastq1_stem + ".reads_out.txt")
        n_reads_in = util.parse_count_file(n_reads_in_path)
        n_reads_out = util.parse_count_file(n_reads_out_path)
        n_reads_delta = n_reads_in - n_reads_out
        try:
            proportion_removed = round(n_reads_delta / n_reads_in, 4)
        except ArithmeticError:  # ZeroDivisionError
            proportion_removed = 0
        stats[(fastq1.name, fastq2.name)] = SampleReport(
            reads_in=n_reads_in,
            reads_out=n_reads_out,
            delta=n_reads_delta,
            proportion_removed=proportion_removed,
        )
    return stats


def dehost_paired_fastqs(
    fastqs: list[tuple[Path, Path]],
    out_dir: Path = CWD,
    threads: int = THREADS,
    aligner: ALIGNERS = ALIGNERS.bowtie2,
):
    aligner.value.check()
    backend_cmds = {
        p: aligner.value.gen_paired_dehost_cmd(
            p[0], p[1], out_dir=out_dir, threads=threads
        )
        for p in fastqs
    }
    util.run_bash_parallel(backend_cmds, cwd=out_dir)
    stats = gather_stats(fastqs, out_dir=out_dir)
