import logging
import gzip
import shutil

from dataclasses import dataclass
from pathlib import Path

from hostile_eit import util, __version__
from hostile_eit.aligner import ALIGNER


logging.basicConfig(
    format="%(asctime)s %(levelname)s: %(message)s",
    datefmt="%H:%M:%S",
    level=logging.INFO,
)

logging.getLogger("httpx").setLevel(logging.WARNING)


@dataclass
class SampleReport:
    version: str
    aligner: str
    index: str
    options: list[str]
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
    rename: bool,
    reorder: bool,
    fastqs: list[Path],
    out_dir: Path,
    aligner: str,
    invert: bool,
    index: str,
) -> list[dict[str, str | int | float | list[str]]]:
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
        options = [
            k
            for k, v in {"rename": rename, "reorder": reorder, "invert": invert}.items()
            if v
        ]
        report = SampleReport(
            version=__version__,
            aligner=aligner,
            index=index,
            options=options,
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
    rename: bool,
    reorder: bool,
    fastqs: list[tuple[Path, Path]],
    out_dir: Path,
    aligner: str,
    index: str,
    invert: bool,
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
        options = [
            k
            for k, v in {"rename": rename, "reorder": reorder, "invert": invert}.items()
            if v
        ]
        stats.append(
            SampleReport(
                version=__version__,
                aligner=aligner,
                index=index,
                options=options,
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


def clean_fastqs(
    fastqs: list[Path],
    index: str = util.DEFAULT_INDEX_NAME,
    invert: bool = False,
    rename: bool = False,
    reorder: bool = False,
    out_dir: Path = util.CWD,
    aligner: ALIGNER = ALIGNER.minimap2,
    aligner_args: str = "",
    threads: int = util.THREADS,
    force: bool = False,
    offline: bool = False,
):
    logging.debug(f"clean_fastqs() {threads=}")
    if aligner == ALIGNER.bowtie2:
        logging.info(f"Hostile version {__version__}. Mode: short read (Bowtie2)")
    elif aligner == ALIGNER.minimap2:
        logging.info(f"Hostile version {__version__}. Mode: long read (Minimap2)")
    fastqs = [Path(path).absolute() for path in fastqs]
    if not all(fastq.is_file() for fastq in fastqs):
        raise FileNotFoundError("One or more fastq files do not exist")
    Path(out_dir).mkdir(exist_ok=True, parents=True)
    index_path = aligner.value.check_index(index, offline=offline)
    backend_cmds = [
        aligner.value.gen_clean_cmd(
            fastq=fastq,
            out_dir=out_dir,
            index_path=index_path,
            invert=invert,
            rename=rename,
            reorder=reorder,
            aligner_args=aligner_args,
            threads=threads,
            force=force,
            offline=offline,
        )
        for fastq in fastqs
    ]
    logging.debug(f"{backend_cmds=}")
    logging.info("Cleaning…")
    util.run_bash_parallel(backend_cmds, description="Cleaning")
    stats = gather_stats(
        rename=rename,
        reorder=reorder,
        fastqs=fastqs,
        out_dir=out_dir,
        aligner=aligner.name,
        index=index,
        invert=invert,
    )
    util.fix_empty_fastqs(stats)
    logging.info("Cleaning complete")
    return stats


def clean_paired_fastqs(
    fastqs: list[tuple[Path, Path]],
    index: str = util.DEFAULT_INDEX_NAME,
    invert: bool = False,
    rename: bool = False,
    reorder: bool = False,
    out_dir: Path = util.CWD,
    aligner: ALIGNER = ALIGNER.bowtie2,
    aligner_args: str = "",
    threads: int = util.THREADS,
    force: bool = False,
    offline: bool = False,
):
    logging.debug(f"clean_paired_fastqs() {threads=}")
    if aligner == ALIGNER.bowtie2:
        logging.info(
            f"Hostile version {__version__}. Mode: paired short read (Bowtie2)"
        )
    elif aligner == ALIGNER.minimap2:
        logging.info(
            f"Hostile version {__version__}. Mode: paired short read (Minimap2)"
        )
    fastqs = [
        (Path(path1).absolute(), Path(path2).absolute()) for path1, path2 in fastqs
    ]
    if not all(path.is_file() for fastq_pair in fastqs for path in fastq_pair):
        raise FileNotFoundError("One or more fastq files do not exist")
    Path(out_dir).mkdir(exist_ok=True, parents=True)
    index_path = aligner.value.check_index(index, offline=offline)
    backend_cmds = [
        aligner.value.gen_paired_clean_cmd(
            fastq1=pair[0],
            fastq2=pair[1],
            out_dir=out_dir,
            index_path=index_path,
            invert=invert,
            rename=rename,
            reorder=reorder,
            aligner_args=aligner_args,
            threads=threads,
            force=force,
            offline=offline,
        )
        for pair in fastqs
    ]
    logging.debug(f"{backend_cmds=}")
    logging.info("Cleaning…")
    util.run_bash_parallel(backend_cmds, description="Cleaning")
    stats = gather_stats_paired(
        rename=rename,
        reorder=reorder,
        fastqs=fastqs,
        out_dir=out_dir,
        aligner=aligner.name,
        index=index,
        invert=invert,
    )
    util.fix_empty_fastqs(stats)
    logging.info("Cleaning complete")
    return stats


def mask(
    reference: Path,
    target: Path,
    out_dir=Path("masked"),
    kmer_length: int = 150,
    kmer_step: int = 10,
    threads: int = util.CPU_COUNT,
) -> tuple[Path, int, int]:
    """Mask a fasta[.gz] reference genome against fasta.[gz] target genomes"""
    ref_path, target_path = Path(reference), Path(target)
    out_dir.mkdir(exist_ok=True, parents=True)
    kmers_path = out_dir / "kmers.fasta.gz"
    mask_path = out_dir / "mask.bed"
    masked_ref_path = out_dir / f"{out_dir.name}.fa"
    masked_ref_index_path = out_dir / out_dir.name

    if ref_path.suffix == ".gz":  # Decompress reference if necessary
        new_ref_path = out_dir / ref_path.stem
        logging.info(f"Decompressing reference into {new_ref_path}")
        with gzip.open(ref_path, "rb") as in_fh:
            with open(new_ref_path, "wb") as out_fh:
                shutil.copyfileobj(in_fh, out_fh)
        ref_path = new_ref_path

    logging.info(
        f"k-merising target genome(s) {target_path} ({kmer_length=}, {kmer_step=})"
    )
    util.kmerise(
        in_path=target_path,
        out_path=out_dir / "kmers.fasta.gz",
        k=kmer_length,
        step=kmer_step,
    )

    mask_cmd = (
        f"minimap2 -N 10000 -p 0.0001 -t {threads} '{ref_path}' '{kmers_path}'"
        f" | awk -v OFS='\t' '{{print $6, $8, $9}}'"
        f" | sort -k1,1 -k2,2n"
        f" | bedtools merge -i stdin > '{mask_path}'"
    )
    logging.info(f"Masking ({mask_cmd})")
    mask_cmd_run = util.run(mask_cmd)
    if mask_cmd_run.stderr:
        logging.info(mask_cmd_run.stderr.strip())

    count_cmd = f"awk '{{sum += $3 - $2}} END {{print sum}}' {mask_path}"
    logging.info(f"Counting masked sites ({mask_cmd})")
    count_cmd_run = util.run(count_cmd)
    n_masked_positions = int(count_cmd_run.stdout.strip())
    if count_cmd_run.stderr:
        logging.info(count_cmd_run.stderr.strip())
    logging.info(f"Masked {n_masked_positions} positions")

    apply_cmd = (
        f"bedtools maskfasta"
        f" -fi '{ref_path}' -bed '{mask_path}' -fo '{masked_ref_path}'"
    )
    logging.info(f"Applying mask ({apply_cmd=})")
    apply_cmd_run = util.run(apply_cmd)
    if apply_cmd_run.stderr:
        logging.info(apply_cmd_run.stderr)

    build_masked_index_cmd = f"bowtie2-build --threads '{threads}' '{masked_ref_path}' '{masked_ref_index_path}'"
    logging.info(f"Indexing masked reference ({build_masked_index_cmd})")
    build_masked_index_cmd_run = util.run(build_masked_index_cmd)
    if build_masked_index_cmd_run.stderr:
        logging.info(build_masked_index_cmd_run.stderr.strip())

    logging.info(f"Masked {n_masked_positions} positions")
    logging.info(f"Masked reference genome: {masked_ref_path}")
    logging.info(f"Masked reference Bowtie2 index: {masked_ref_index_path}")

    return masked_ref_path, masked_ref_index_path, n_masked_positions
