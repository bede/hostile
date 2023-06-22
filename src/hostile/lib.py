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


CWD = Path.cwd().resolve()
XDG_DATA_DIR = Path(user_data_dir("hostile", "Bede Constantinides"))
THREADS = multiprocessing.cpu_count()


ALIGNERS = Enum(
    "Aligner",
    {
        "bowtie2": Aligner(
            name="Bowtie2",
            short_name="bt2",
            bin_path=Path("bowtie2"),
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
    fastq1_in_name: str
    fastq2_in_name: str
    fastq1_in_path: str
    fastq2_in_path: str
    fastq1_out_name: str
    fastq2_out_name: str
    fastq1_out_path: str
    fastq2_out_path: str
    reads_in: int
    reads_out: int
    reads_removed: int
    reads_removed_proportion: float


def gather_stats(
    fastqs: list[tuple[Path, Path]], out_dir: Path
) -> dict[str, dict[str : str | int | float]]:
    stats = []
    for fastq1, fastq2 in fastqs:
        fastq1_stem = util.fastq_path_to_stem(fastq1)
        fastq2_stem = util.fastq_path_to_stem(fastq2)
        fastq1_out_path = out_dir / f"{fastq1_stem}.clean_1.fastq.gz"
        fastq2_out_path = out_dir / f"{fastq2_stem}.clean_2.fastq.gz"
        n_reads_in_path = out_dir / (fastq1_stem + ".reads_in.txt")
        n_reads_out_path = out_dir / (fastq1_stem + ".reads_out.txt")
        n_reads_in = util.parse_count_file(n_reads_in_path)
        n_reads_out = util.parse_count_file(n_reads_out_path)
        n_reads_removed = n_reads_in - n_reads_out
        n_reads_in_path.unlink()
        n_reads_out_path.unlink()
        try:
            proportion_removed = round(n_reads_removed / n_reads_in, 5)
        except ArithmeticError:  # ZeroDivisionError
            proportion_removed = float(0)
        stats.append(
            SampleReport(
                fastq1_in_name=fastq1.name,
                fastq2_in_name=fastq2.name,
                fastq1_in_path=str(fastq1),
                fastq2_in_path=str(fastq2),
                fastq1_out_name=fastq1_out_path.name,
                fastq2_out_name=fastq2_out_path.name,
                fastq1_out_path=str(fastq1_out_path),
                fastq2_out_path=str(fastq2_out_path),
                reads_in=n_reads_in,
                reads_out=n_reads_out,
                reads_removed=n_reads_removed,
                reads_removed_proportion=proportion_removed,
            ).__dict__
        )
    return stats


def clean_paired_fastqs(
    fastqs: list[tuple[Path, Path]],
    out_dir: Path = CWD,
    threads: int = THREADS,
    aligner: ALIGNERS = ALIGNERS.bowtie2,
):
    Path(out_dir).mkdir(exist_ok=True, parents=True)
    try:
        aligner.value.check()
    except Exception as e:
        previous_aligner = aligner.name
        if aligner == ALIGNERS.bowtie2:
            aligner = ALIGNERS.minimap2
        elif aligner == ALIGNERS.minimap2:
            aligner = ALIGNERS.bowtie2
        logging.warning(f"Using {aligner.name} instead of {previous_aligner.name}")
        aligner.value.check()

    backend_cmds = {
        p: aligner.value.gen_paired_clean_cmd(
            Path(p[0]), Path(p[1]), out_dir=out_dir, threads=threads
        )
        for p in fastqs
    }
    util.run_bash_parallel(backend_cmds, description="Cleaning")
    stats = gather_stats(fastqs, out_dir=out_dir)
    return stats
