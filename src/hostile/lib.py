import logging
import gzip
import shutil

from dataclasses import dataclass
from pathlib import Path

from hostile import util, __version__
from hostile.aligner import ALIGNER, get_mmi_path


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
    reads_in: int
    reads_out: int
    reads_removed: int
    reads_removed_proportion: float
    fastq2_in_name: str | None = None
    fastq2_in_path: str | None = None
    fastq1_out_name: str | None = None
    fastq1_out_path: str | None = None
    fastq2_out_name: str | None = None
    fastq2_out_path: str | None = None


def gather_stats(
    fastqs: list[Path | None],
    aligner: str,
    index: str,
    invert: bool,
    rename: bool,
    reorder: bool,
    casava: bool,
    stdin: bool,
    stdout: bool,
    output: Path,
) -> list[dict[str, str | int | float | list[str]]]:
    stats = []
    logging.debug(f"gather_stats() {fastqs=}")
    for fastq1 in fastqs:
        fastq1_stem = util.fastq_path_to_stem(fastq1)
        fastq1_out_path = output / f"{fastq1_stem}.clean.fastq.gz"
        n_reads_in_path = output / (fastq1_stem + ".reads_in.txt")
        n_reads_out_path = output / (fastq1_stem + ".reads_out.txt")
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
            for k, v in {
                "invert": invert,
                "rename": rename,
                "reorder": reorder,
                "casava": casava,
                "stdin": stdin,
                "stdout": stdout,
            }.items()
            if v
        ]
        report = SampleReport(
            version=__version__,
            aligner=aligner,
            index=index,
            options=options,
            fastq1_in_name=Path(fastq1).name if not stdin else None,
            fastq1_in_path=str(fastq1) if not stdin else None,
            fastq1_out_name=fastq1_out_path.name if not stdout else None,
            fastq1_out_path=str(fastq1_out_path) if not stdout else None,
            reads_in=n_reads_in,
            reads_out=n_reads_out,
            reads_removed=n_reads_removed,
            reads_removed_proportion=proportion_removed,
        ).__dict__
        stats.append({k: v for k, v in report.items() if v is not None})
    return stats


def gather_stats_paired(
    fastqs: list[tuple[Path | None, Path | None]],
    aligner: str,
    index: str,
    invert: bool,
    rename: bool,
    reorder: bool,
    casava: bool,
    stdin: bool,
    stdout: bool,
    output: Path,
) -> list[dict[str, str | int | float | list[str]]]:
    stats = []
    logging.debug(f"gather_stats_paired() {fastqs=}")
    for fastq1, fastq2 in fastqs:
        fastq1_stem = util.fastq_path_to_stem(fastq1)
        fastq2_stem = util.fastq_path_to_stem(fastq2)
        fastq1_out_path = output / f"{fastq1_stem}.clean_1.fastq.gz"
        fastq2_out_path = output / f"{fastq2_stem}.clean_2.fastq.gz"
        n_reads_in_path = output / (fastq1_stem + ".reads_in.txt")
        n_reads_out_path = output / (fastq1_stem + ".reads_out.txt")
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
            for k, v in {
                "invert": invert,
                "rename": rename,
                "reorder": reorder,
                "casava": casava,
                "stdin": stdin,
                "stdout": stdout,
            }.items()
            if v
        ]
        report = SampleReport(
            version=__version__,
            aligner=aligner,
            index=index,
            options=options,
            fastq1_in_name=Path(fastq1).name if not stdin else None,
            fastq2_in_name=Path(fastq2).name if not stdin else None,
            fastq1_in_path=str(fastq1) if not stdin else None,
            fastq2_in_path=str(fastq2) if not stdin else None,
            fastq1_out_name=fastq1_out_path.name if not stdout else None,
            fastq2_out_name=fastq2_out_path.name if not stdout else None,
            fastq1_out_path=str(fastq1_out_path) if not stdout else None,
            fastq2_out_path=str(fastq2_out_path) if not stdout else None,
            reads_in=n_reads_in,
            reads_out=n_reads_out,
            reads_removed=n_reads_removed,
            reads_removed_proportion=proportion_removed,
        ).__dict__
        stats.append({k: v for k, v in report.items() if v is not None})

    return stats


def clean_fastqs(
    fastqs: list[Path],
    aligner: ALIGNER = ALIGNER.minimap2,
    index: str = util.DEFAULT_INDEX_NAME,
    invert: bool = False,
    rename: bool = False,
    reorder: bool = False,
    casava: bool = False,
    output: Path = util.CWD,
    aligner_args: str = "",
    threads: int = util.CPU_COUNT,
    force: bool = False,
    airplane: bool = False,
):
    stdin = str(fastqs[0]) == "-"
    stdout = str(output) == "-"
    output = Path(util.CWD) if stdout else Path(output)
    aligner_threads, compression_threads = util.allocate_threads(threads, stdout=stdout)
    logging.debug(
        f"clean_fastqs() {threads=} {aligner_threads=} {compression_threads=}"
        f" {util.CACHE_DIR=} {util.INDEX_REPOSITORY_URL=}"
    )
    if aligner == ALIGNER.bowtie2:
        logging.info(
            f"Hostile v{__version__}. Mode: short read {'from stdin ' if stdin else ''}(Bowtie2)"
        )
    elif aligner == ALIGNER.minimap2:
        logging.info(
            f"Hostile v{__version__}. Mode: long read {'from stdin ' if stdin else ''}(Minimap2)"
        )
    if not stdin:
        fastqs = [Path(path).absolute() for path in fastqs]
        if not all(fastq.is_file() for fastq in fastqs):
            logging.debug(f"{fastqs=}")
            raise FileNotFoundError("One or more fastq files do not exist")
    index_path = aligner.value.check_index(index, airplane=airplane)
    backend_cmds = [
        aligner.value.gen_clean_cmd(
            fastq=fastq,
            index_path=index_path,
            invert=invert,
            rename=rename,
            reorder=reorder,
            casava=casava,
            stdin=stdin,
            stdout=stdout,
            output=output,
            aligner_args=aligner_args,
            aligner_threads=aligner_threads,
            compression_threads=compression_threads,
            force=force,
        )
        for fastq in fastqs
    ]
    logging.debug(f"{backend_cmds=}")
    logging.info("Cleaning…")
    if stdin:
        util.run_bash(backend_cmds[0], stdin=True)
        fastqs[0] = "stdin"
    else:
        util.run_bash_parallel(backend_cmds, description="Cleaning")
    stats = gather_stats(
        fastqs=fastqs,
        aligner=aligner.name,
        index=index,
        invert=invert,
        rename=rename,
        reorder=reorder,
        casava=casava,
        stdin=stdin,
        stdout=stdout,
        output=output,
    )
    util.fix_empty_fastqs(stats)
    logging.info("Cleaning complete")
    return stats


def clean_paired_fastqs(
    fastqs: list[tuple[Path, Path]],
    aligner: ALIGNER = ALIGNER.bowtie2,
    index: str = util.DEFAULT_INDEX_NAME,
    invert: bool = False,
    rename: bool = False,
    reorder: bool = False,
    casava: bool = False,
    output: Path = util.CWD,
    aligner_args: str = "",
    threads: int = util.CPU_COUNT,
    force: bool = False,
    airplane: bool = False,
):
    stdin = str(fastqs[0][0]) == "-"
    stdout = str(output) == "-"
    output = Path(util.CWD) if stdout else Path(output)
    aligner_threads, compression_threads = util.allocate_threads(threads, stdout=stdout)
    logging.debug(
        f"clean_paired_fastqs() {threads=} {aligner_threads=} {compression_threads=}"
        f" {util.CACHE_DIR=} {util.INDEX_REPOSITORY_URL=}"
    )
    if aligner == ALIGNER.bowtie2:
        logging.info(
            f"Hostile v{__version__}. Mode: paired short read {'from stdin ' if stdin else ''}(Bowtie2)"
        )
    elif aligner == ALIGNER.minimap2:
        logging.info(
            f"Hostile v{__version__}. Mode: paired short read {'from stdin ' if stdin else ''}(Minimap2)"
        )
    if not stdin:
        fastqs = [
            (Path(path1).absolute(), Path(path2).absolute()) for path1, path2 in fastqs
        ]
        if not all(path.is_file() for fastq_pair in fastqs for path in fastq_pair):
            raise FileNotFoundError("One or more fastq files do not exist")
    index_path = aligner.value.check_index(index, airplane=airplane)
    backend_cmds = [
        aligner.value.gen_paired_clean_cmd(
            fastq1=fastq_pair[0],
            fastq2=fastq_pair[1],
            index_path=index_path,
            invert=invert,
            rename=rename,
            reorder=reorder,
            casava=casava,
            stdin=stdin,
            stdout=stdout,
            output=output,
            aligner_args=aligner_args,
            aligner_threads=aligner_threads,
            compression_threads=compression_threads,
            force=force,
        )
        for fastq_pair in fastqs
    ]
    logging.debug(f"{backend_cmds=}")
    logging.info("Cleaning…")
    if stdin:
        util.run_bash(backend_cmds[0], stdin=True)
        fastqs[0] = ("stdin", "stdin")
    else:
        util.run_bash_parallel(backend_cmds, description="Cleaning")
    stats = gather_stats_paired(
        fastqs=fastqs,
        aligner=aligner.name,
        index=index,
        invert=invert,
        rename=rename,
        reorder=reorder,
        casava=casava,
        stdin=stdin,
        stdout=stdout,
        output=output,
    )
    util.fix_empty_fastqs(stats)
    logging.info("Cleaning complete")
    return stats


def mask(
    reference: Path,
    target: Path,
    output=Path("masked"),
    kmer_length: int = 150,
    kmer_step: int = 10,
    threads: int = util.CPU_COUNT,
) -> tuple[Path, int, int]:
    """Mask a fasta[.gz] reference genome against fasta.[gz] target genomes"""
    ref_path, target_path = Path(reference), Path(target)
    output.mkdir(exist_ok=True, parents=True)
    kmers_path = output / "kmers.fasta.gz"
    mask_path = output / "mask.bed"
    masked_ref_path = output / f"{output.name}.fa"
    masked_ref_index_path = output / output.name

    if ref_path.suffix == ".gz":  # Decompress reference if necessary
        new_ref_path = output / ref_path.stem
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
        out_path=output / "kmers.fasta.gz",
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

    masked_ref_mmi_path = get_mmi_path(masked_ref_path)
    build_mm2_masked_index_cmd = (
        f"minimap2 -d '{masked_ref_mmi_path}' '{masked_ref_path}'"
    )
    logging.info(f"Building Minimap2 index ({build_mm2_masked_index_cmd})")
    build_masked_index_cmd_run = util.run(build_mm2_masked_index_cmd)

    build_bt2_masked_index_cmd = f"bowtie2-build --threads '{threads}' '{masked_ref_path}' '{masked_ref_index_path}'"
    logging.info(f"Building Bowtie2 index ({build_bt2_masked_index_cmd})")
    build_masked_index_cmd_run = util.run(build_bt2_masked_index_cmd)
    if build_masked_index_cmd_run.stderr:
        logging.info(build_masked_index_cmd_run.stderr.strip())

    logging.info(f"Masked {n_masked_positions} positions")
    logging.info(f"Masked reference: {masked_ref_path}")
    logging.info(f"Masked Minimap2 index: {masked_ref_mmi_path}")
    logging.info(f"Masked Bowtie2 index: {masked_ref_index_path}")

    return masked_ref_path, masked_ref_index_path, n_masked_positions


def fetch_index(
    name: str = util.DEFAULT_INDEX_NAME,
    minimap2: bool = False,
    bowtie2: bool = False,
) -> None:
    if minimap2 or (not minimap2 and not bowtie2):
        logging.info(f"Looking for Minimap2 index {name}")
        ALIGNER.minimap2.value.check_index(name)
    if bowtie2 or (not minimap2 and not bowtie2):
        logging.info(f"Looking for Bowtie2 index {name}")
        ALIGNER.bowtie2.value.check_index(name)


def list_indexes(airplane: bool = False):
    total_size_gb = sum(
        f.stat().st_size / (1024**3) for f in util.CACHE_DIR.glob("**/*") if f.is_file()
    )
    logging.info(f"Remote index: {util.INDEX_REPOSITORY_URL}/manifest.json")
    logging.info(f"Local cache: '{util.CACHE_DIR}' ({total_size_gb:.1f}GB total)")
    if not airplane:
        manifest = util.fetch_manifest()
        for name in manifest.keys():
            print(f"Remote\t{name}")
    unique_prefixes = set()
    for f in util.CACHE_DIR.iterdir():
        if f.is_file():
            name = f.name
            if name.endswith(".fa.gz"):
                prefix = name[:-6]
            elif name.endswith(".bt2"):
                prefix = name[:-6]
                if prefix.endswith(".rev"):
                    prefix = prefix[:-4]
            else:
                continue
            unique_prefixes.add(prefix)
    for prefix in unique_prefixes:
        suffixes = []
        fa_exists = (util.CACHE_DIR / f"{prefix}.fa.gz").exists()
        mmi_exists = (util.CACHE_DIR / f"{prefix}.mmi").exists()
        bt2_exists = (util.CACHE_DIR / f"{prefix}.1.bt2").exists()
        if fa_exists:
            suffixes.append(f"Minimap2{' with MMI' if mmi_exists else ''}")
        if bt2_exists:
            suffixes.append("Bowtie2")
        suffix_str = ", ".join(suffixes)
        print(f"Local\t{prefix} ({suffix_str})")


def delete_index(name: str = "", all: bool = False, mmi: bool = False) -> None:
    total_size_gb = sum(
        f.stat().st_size / (1024**3) for f in util.CACHE_DIR.glob("**/*") if f.is_file()
    )
    logging.info(
        f"Local cache: '{util.CACHE_DIR}' ({total_size_gb:.1f}GB total before deletion)"
    )
    for f in util.CACHE_DIR.iterdir():
        if f.is_file() and mmi and f.name.endswith(".mmi"):
            f.unlink()
            logging.info(f"Deleted {f}")
    for f in util.CACHE_DIR.iterdir():
        if f.is_file():
            if all or (name and f.name.startswith(name)):
                f.unlink()
                logging.info(f"Deleted {f}")
    if not name and not all and not mmi:
        logging.error("Provide an index name (--name NAME) or --all")
