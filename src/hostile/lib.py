import logging
import gzip
import shutil

from dataclasses import dataclass
from pathlib import Path

import dnaio

from hostile import util, __version__
from hostile.aligner import ALIGNER


logging.basicConfig(
    format="%(asctime)s %(levelname)s: %(message)s",
    datefmt="%H:%M:%S",
    level=logging.INFO,
)

# logging.getLogger("httpx").setLevel(logging.WARNING)


@dataclass
class SampleReport:
    version: str
    aligner: str
    index: str
    options: str
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
    index: Path | None,
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
        index_fmt = (
            index
            if index
            else Path(ALIGNER[aligner].value.data_dir)
            / Path(ALIGNER[aligner].value.idx_name)
        )
        options = [
            k
            for k, v in {"rename": rename, "reorder": reorder, "invert": invert}.items()
            if v
        ]
        report = SampleReport(
            version=__version__,
            aligner=aligner,
            index=str(index_fmt),
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
    index: Path | None,
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
        index_fmt = (
            index
            if index
            else Path(ALIGNER[aligner].value.data_dir)
            / Path(ALIGNER[aligner].value.idx_name)
        )
        options = [
            k
            for k, v in {"rename": rename, "reorder": reorder, "invert": invert}.items()
            if v
        ]
        stats.append(
            SampleReport(
                version=__version__,
                aligner=aligner,
                index=str(index_fmt),
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
    index: str = "human-t2t-hla",
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
        logging.info("Using Bowtie2")
    elif aligner == ALIGNER.minimap2:
        logging.info("Using Minimap2's long read preset")
    fastqs = [Path(path).absolute() for path in fastqs]
    if not all(fastq.is_file() for fastq in fastqs):
        raise FileNotFoundError("One or more fastq files do not exist")
    Path(out_dir).mkdir(exist_ok=True, parents=True)
    backend_cmds = [
        aligner.value.gen_clean_cmd(
            fastq=fastq,
            out_dir=out_dir,
            index=index,
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
    logging.info("Finished cleaning")
    return stats


def clean_paired_fastqs(
    fastqs: list[tuple[Path, Path]],
    index: str = "human-t2t-hla",
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
        logging.info("Using Bowtie2 (paired reads)")
    elif aligner == ALIGNER.minimap2:
        logging.info("Using Minimap2's short read preset (paired reads)")
    fastqs = [
        (Path(path1).absolute(), Path(path2).absolute()) for path1, path2 in fastqs
    ]
    if not all(path.is_file() for fastq_pair in fastqs for path in fastq_pair):
        raise FileNotFoundError("One or more fastq files do not exist")
    Path(out_dir).mkdir(exist_ok=True, parents=True)
    backend_cmds = [
        aligner.value.gen_paired_clean_cmd(
            fastq1=pair[0],
            fastq2=pair[1],
            out_dir=out_dir,
            index=index,
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
    logging.info("Finished cleaning")
    return stats


def mask(
    reference: Path,
    target: Path,
    out_dir=Path("masked"),
    k: int = 150,
    i: int = 50,
    threads: int = util.CPU_COUNT,
) -> tuple[Path, int, int]:
    """Mask a fasta[.gz] reference genome against fasta.[gz] target genomes"""
    ref_path, target_path = Path(reference), Path(target)
    out_dir.mkdir(exist_ok=True, parents=True)
    ref_index_path = out_dir / "existing"
    kmers_path = out_dir / "kmers.fasta.gz"
    alignments_path = out_dir / "alignments.sam"
    mask_path = out_dir / "mask.bed"
    masked_ref_path = out_dir / "masked.fa"
    masked_ref_index_path = out_dir / "masked"
    masked_alignments_path = out_dir / "masked-alignments.sam"

    if ref_path.suffix == ".gz":  # Decompress reference if necessary
        new_ref_path = out_dir / ref_path.stem
        logging.info(f"Decompressing reference into {new_ref_path}")
        with gzip.open(ref_path, "rb") as in_fh:
            with open(new_ref_path, "wb") as out_fh:
                shutil.copyfileobj(in_fh, out_fh)
        ref_path = new_ref_path

    build_existing_cmd = (
        f"bowtie2-build --threads '{threads}' '{ref_path}' '{out_dir}/existing'"
    )
    logging.info(f"Indexing existing reference ({build_existing_cmd})")
    build_existing_cmd_run = util.run(build_existing_cmd)
    if build_existing_cmd_run.stderr.strip():
        logging.info(build_existing_cmd_run.stderr.strip())

    logging.info(f"k-merising target genome(s) {target_path} ({k=}, {i=})")
    kmerise(path=target_path, out_dir=out_dir, k=k, i=i)

    align_cmd = (
        f"bowtie2 -a -p '{threads}'"
        f" -x '{ref_index_path}'"
        f" -f '{kmers_path}'"
        f" > '{alignments_path}'"
    )
    logging.info(f"Aligning target k-mers to existing reference ({align_cmd})")
    align_cmd_run = util.run(align_cmd)
    if align_cmd_run.stderr:
        logging.info(align_cmd_run.stderr.strip())

    count_alignments_cmd = (  # Exclude unmapped reads and secondary alignments
        f"samtools view -c -F 0x904 {alignments_path}"
    )
    logging.info(
        f"Counting target k-mers aligned to existing reference ({count_alignments_cmd})"
    )
    count_alignments_cmd_run = util.run(count_alignments_cmd)
    if count_alignments_cmd_run.stderr:
        logging.info(count_alignments_cmd_run.stderr)
    n_alignments = int(count_alignments_cmd_run.stdout.strip())
    logging.info(f"{n_alignments} k-mers aligned before masking")

    make_cmd = (
        f"samtools view -F 4 -bS '{alignments_path}'"
        f" | samtools sort -"
        f" | bedtools bamtobed -i stdin"
        f" | bedtools merge -i stdin"
        f" > '{mask_path}'"
    )
    logging.info(f"Making mask ({make_cmd=})")
    make_cmd_run = util.run(make_cmd)
    if make_cmd_run.stderr:
        logging.info(make_cmd_run.stderr)

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

    align_masked_cmd = (
        f"bowtie2 -a -p '{threads}'"
        f" -x '{masked_ref_index_path}'"
        f" -f '{kmers_path}'"
        f" > '{masked_alignments_path}'"
    )
    logging.info(f"Aligning target k-mers to masked reference ({align_masked_cmd})")
    align_masked_cmd_run = util.run(align_masked_cmd)
    if align_masked_cmd_run.stderr:
        logging.info(align_masked_cmd_run.stderr.strip())

    count_masked_alignments_cmd = f"samtools view -c -F 0x904 {masked_alignments_path}"  # Exclude secondaries (0x100), supplementaries (0x800), and unmapped (0x4)
    logging.info(
        f"Counting target k-mers aligned to masked reference ({count_masked_alignments_cmd})"
    )
    count_masked_alignments_cmd_run = util.run(count_masked_alignments_cmd)
    if count_masked_alignments_cmd_run.stderr:
        logging.info(count_masked_alignments_cmd_run.stderr)
    n_masked_alignments = int(count_masked_alignments_cmd_run.stdout.strip())
    logging.info(
        f"{n_masked_alignments} k-mers aligned after masking ({n_alignments} aligned before masking)"
    )
    logging.info(
        f"Masked genome path (for use with long reads / Minimap2): {masked_ref_path.absolute()}"
    )
    logging.info(
        f"Masked Bowtie2 index path (for use with short reads): {masked_ref_index_path.absolute()} (multiple files)"
    )

    return masked_ref_path, n_alignments, n_masked_alignments


def get_default_reference_filenames() -> list[Path]:
    return [ALIGNER.minimap2.value.ref_archive_fn, ALIGNER.bowtie2.value.idx_archive_fn]


def fetch_reference(filename: str) -> None:
    util.download(url=f"{util.BUCKET_URL}/{filename}", path=Path(filename))
    if filename.endswith(".tar"):
        logging.info("Extracting…")
        util.untar_file(input_path=Path(filename), output_path=Path("."))
        logging.info(f"Downloaded and extracted {filename}")
    else:
        logging.info(f"Downloaded {filename}")


def kmerise(path: Path, out_dir: Path, k: int, i: int) -> Path:
    out_path = out_dir / "kmers.fasta.gz"
    with dnaio.open(path) as reader, dnaio.open(out_path, mode="w") as writer:
        for r in reader:
            for offset in range(0, len(r.sequence) - k + 1, i):
                kmer = r.sequence[offset : offset + k]
                name = r.name.partition(" ")[0]
                kmer_id = f"{name}_{offset}"
                writer.write(dnaio.SequenceRecord(kmer_id, kmer))
    return out_path.absolute()
