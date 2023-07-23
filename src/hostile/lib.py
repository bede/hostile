import logging
import gzip
import multiprocessing
import shutil
import subprocess

from enum import Enum
from dataclasses import dataclass
from pathlib import Path

from platformdirs import user_data_dir

from hostile import util
from hostile.aligner import Aligner


logging.basicConfig(
    format="%(asctime)s %(levelname)s: %(message)s",
    datefmt="%H:%M:%S",
    level=logging.INFO,
)


CWD = Path.cwd().resolve()
XDG_DATA_DIR = Path(user_data_dir("hostile", "Bede Constantinides"))
THREADS = multiprocessing.cpu_count()


ALIGNER = Enum(
    "Aligner",
    {
        "bowtie2": Aligner(
            name="Bowtie2",
            short_name="bt2",
            bin_path=Path("bowtie2"),
            # cdn_base_url="http://localhost:8000",  # python -m http.server
            cdn_base_url=f"https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o",
            data_dir=XDG_DATA_DIR,
            cmd=("{BIN_PATH} -x '{INDEX_PATH}' -U '{FASTQ}'" " -k 1 --mm -p {THREADS}"),
            paired_cmd=(
                "{BIN_PATH} -x '{INDEX_PATH}' -1 '{FASTQ1}' -2 '{FASTQ2}'"
                " -k 1 --mm -p {THREADS}"
            ),
            idx_archive_fn="human-t2t-hla.tar",
            idx_name="human-t2t-hla",
            idx_paths=(
                XDG_DATA_DIR / "human-t2t-hla.1.bt2",
                XDG_DATA_DIR / "human-t2t-hla.2.bt2",
                XDG_DATA_DIR / "human-t2t-hla.3.bt2",
                XDG_DATA_DIR / "human-t2t-hla.4.bt2",
                XDG_DATA_DIR / "human-t2t-hla.rev.1.bt2",
                XDG_DATA_DIR / "human-t2t-hla.rev.2.bt2",
            ),
        ),
        "minimap2": Aligner(
            name="Minimap2",
            short_name="mm2",
            bin_path=Path("minimap2"),
            # cdn_base_url="http://localhost:8000",  # python -m http.server
            cdn_base_url=f"https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o",
            data_dir=XDG_DATA_DIR,
            cmd="{BIN_PATH} -ax map-ont --secondary no -t {THREADS} '{REF_ARCHIVE_PATH}' '{FASTQ}'",
            paired_cmd="{BIN_PATH} -ax sr -m 40 --secondary no -t {THREADS} '{REF_ARCHIVE_PATH}' '{FASTQ1}' '{FASTQ2}'",
            ref_archive_fn="human-t2t-hla.fa.gz",
            idx_name="human-t2t-hla.fa.gz",
        ),
    },
)


@dataclass
class SampleReport:
    aligner: str
    index: str
    fastq1_in_name: str
    fastq1_in_path: str
    fastq1_out_name: str
    fastq1_out_path: str
    reads_in: int
    reads_out: int
    reads_removed: int
    reads_removed_proportion: float
    fastq2_in_name: str | None = None
    fastq2_in_path: str | None = None
    fastq2_out_name: str | None = None
    fastq2_out_path: str | None = None


def gather_stats(
    fastqs: list[Path], out_dir: Path, aligner: str, index: Path | None
) -> list[dict[str, str | int | float]]:
    stats = []
    for fastq1 in fastqs:
        fastq1_stem = util.fastq_path_to_stem(fastq1)
        fastq1_out_path = out_dir / f"{fastq1_stem}.clean.fastq.gz"
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
        index_fmt = (
            index
            if index
            else Path(ALIGNER[aligner].value.data_dir)
            / Path(ALIGNER[aligner].value.idx_name)
        )
        report = SampleReport(
            aligner=aligner,
            index=str(index_fmt),
            fastq1_in_name=fastq1.name,
            fastq1_in_path=str(fastq1),
            fastq1_out_name=fastq1_out_path.name,
            fastq1_out_path=str(fastq1_out_path),
            reads_in=n_reads_in,
            reads_out=n_reads_out,
            reads_removed=n_reads_removed,
            reads_removed_proportion=proportion_removed,
        ).__dict__
        stats.append({k: v for k, v in report.items() if v is not None})
    return stats


def gather_stats_paired(
    fastqs: list[tuple[Path, Path]], out_dir: Path, aligner: str, index: Path | None
) -> list[dict[str, str | int | float]]:
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
        index_fmt = (
            index
            if index
            else Path(ALIGNER[aligner].value.data_dir)
            / Path(ALIGNER[aligner].value.idx_name)
        )
        stats.append(
            SampleReport(
                aligner=aligner,
                index=str(index_fmt),
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


def choose_aligner(preferred_aligner: ALIGNER, using_custom_index: bool) -> ALIGNER:
    """Fallback to Minimap2 from Bowtie2 if Bowtie2 isn't installed etc"""
    aligner = preferred_aligner
    try:
        aligner.value.check(using_custom_index=using_custom_index)
    except Exception as e:
        if aligner == ALIGNER.bowtie2:
            aligner = ALIGNER.minimap2
            logging.warning(f"Using Minimap2 instead of Bowtie2")
            aligner.value.check(using_custom_index=using_custom_index)
        else:
            raise e
    return aligner


def clean_fastqs(
    fastqs: list[Path],
    index: Path | None = None,
    rename: bool = False,
    out_dir: Path = CWD,
    aligner: ALIGNER = ALIGNER.minimap2,
    threads: int = THREADS,
    force: bool = False,
):
    if aligner == ALIGNER.bowtie2:
        logging.info("Using Bowtie2")
    elif aligner == ALIGNER.minimap2:
        logging.info("Using Minimap2's long read preset")
    fastqs = [Path(path).resolve() for path in fastqs]
    if not all(fastq.is_file() for fastq in fastqs):
        raise FileNotFoundError("One or more fastq files do not exist")
    Path(out_dir).mkdir(exist_ok=True, parents=True)
    aligner = choose_aligner(aligner, using_custom_index=bool(index))
    backend_cmds = [
        aligner.value.gen_clean_cmd(
            fastq=fastq,
            out_dir=out_dir,
            index=index,
            rename=rename,
            threads=threads,
            force=force,
        )
        for fastq in fastqs
    ]
    logging.info("Cleaning…")
    util.run_bash_parallel(backend_cmds, description="Cleaning")
    stats = gather_stats(fastqs, out_dir=out_dir, aligner=aligner.name, index=index)
    logging.info("Complete")
    return stats


def clean_paired_fastqs(
    fastqs: list[tuple[Path, Path]],
    index: Path | None = None,
    rename: bool = False,
    out_dir: Path = CWD,
    aligner: ALIGNER = ALIGNER.bowtie2,
    threads: int = THREADS,
    force: bool = False,
):
    if aligner == ALIGNER.bowtie2:
        logging.info("Using Bowtie2 (paired reads)")
    elif aligner == ALIGNER.minimap2:
        logging.info("Using Minimap2's short read preset (paired reads)")
    fastqs = [(Path(path1).resolve(), Path(path2).resolve()) for path1, path2 in fastqs]
    if not all(path.is_file() for fastq_pair in fastqs for path in fastq_pair):
        raise FileNotFoundError("One or more fastq files do not exist")
    Path(out_dir).mkdir(exist_ok=True, parents=True)
    aligner = choose_aligner(aligner, using_custom_index=bool(index))
    backend_cmds = [
        aligner.value.gen_paired_clean_cmd(
            fastq1=pair[0],
            fastq2=pair[1],
            out_dir=out_dir,
            index=index,
            rename=rename,
            threads=threads,
            force=force,
        )
        for pair in fastqs
    ]
    logging.debug(f"{backend_cmds=}")
    logging.info("Cleaning…")
    util.run_bash_parallel(backend_cmds, description="Cleaning")
    stats = gather_stats_paired(
        fastqs, out_dir=out_dir, aligner=aligner.name, index=index
    )
    logging.info("Complete")
    return stats


def mask(
    reference: Path, target: Path, out_dir=Path("masked"), threads: int = 1
) -> Path:
    """Mask a fasta[.gz] reference genome against fasta.[gz] target genomes"""
    reference_path, target_path = Path(reference), Path(target)
    out_dir.mkdir(exist_ok=True, parents=True)
    bed_path = out_dir / "mask.bed"
    masked_reference_path = out_dir / "masked.fa"

    if reference_path.suffix == ".gz":  # Decompress reference if necessary
        new_reference_path = out_dir / reference_path.stem
        logging.info(f"Decompressing reference into {new_reference_path}")
        with gzip.open(reference_path, "rb") as in_fh:
            with open(new_reference_path, "wb") as out_fh:
                shutil.copyfileobj(in_fh, out_fh)
        reference_path = new_reference_path

    make_cmd = (
        f"minimap2 -x asm10 -t {threads} '{reference_path}' '{target}'"
        f" | awk -v OFS='\t' '{{print $6, $8, $9}}'"
        f" | sort -k1,1 -k2,2n"
        f" | bedtools merge -i stdin > '{bed_path}'"
    )
    logging.info(f"Making mask ({make_cmd=})")
    make_cmd_run = util.run(make_cmd)
    logging.info(make_cmd_run.stderr)

    apply_cmd = (
        f"bedtools maskfasta"
        f" -fi '{reference_path}' -bed '{bed_path}' -fo '{masked_reference_path}'"
    )
    logging.info(f"Applying mask ({apply_cmd=})")
    apply_cmd_run = util.run(apply_cmd)
    logging.info(apply_cmd_run.stderr)

    return masked_reference_path
