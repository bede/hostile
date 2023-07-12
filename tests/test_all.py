import gzip
import shutil
import subprocess
from pathlib import Path

import pytest

from hostile import lib

data_dir = Path("tests/data")


def run(cmd, cwd=data_dir):  # Helper for CLI testing
    return subprocess.run(
        cmd, cwd=cwd, shell=True, check=True, text=True, capture_output=True
    )


def get_first_line_of_gzip_file(file_path):
    with gzip.open(file_path, "rt") as fh:
        return fh.readline().strip()


def test_version_cli():
    run("hostile --version")


def test_minimal_fastq():
    lib.clean_fastqs(
        fastqs=[data_dir / "h37rv_10.r1.fastq.gz"],
        aligner=lib.ALIGNER.bowtie2,
        index=data_dir / "MN908947/MN908947",
        out_dir=Path("test_minimal_fastq"),
    )
    shutil.rmtree("test_minimal_fastq")


def test_minimal_paired_fastqs():
    lib.clean_paired_fastqs(
        fastqs=[(data_dir / "h37rv_10.r1.fastq.gz", data_dir / "h37rv_10.r2.fastq.gz")],
        aligner=lib.ALIGNER.bowtie2,
        index=data_dir / "MN908947/MN908947",
        out_dir=Path("test_minimal_fastqs"),
    )
    shutil.rmtree("test_minimal_fastqs")


def test_minimal_uncompressed_paired_fastqs():
    lib.clean_paired_fastqs(
        fastqs=[(data_dir / "h37rv_10.r1.fastq", data_dir / "h37rv_10.r2.fastq")],
        aligner=lib.ALIGNER.bowtie2,
        index=data_dir / "MN908947/MN908947",
        out_dir=Path("test_minimal_fastqs"),
    )
    shutil.rmtree("test_minimal_fastqs")


def test_minimal_paired_fastqs_cli():
    run(
        f"hostile clean --index MN908947/MN908947 --fastq1 h37rv_10.r1.fastq.gz --fastq2 h37rv_10.r2.fastq.gz --out-dir test_minimal_fastqs"
    )
    shutil.rmtree(f"{data_dir}/test_minimal_fastqs")


def test_custom_index():
    lib.clean_paired_fastqs(
        fastqs=[(data_dir / "h37rv_10.r1.fastq.gz", data_dir / "h37rv_10.r2.fastq.gz")],
        index=data_dir / "MN908947/MN908947",
        out_dir=Path("test_minimal_fastqs"),
    )
    shutil.rmtree("test_minimal_fastqs")


def test_both_aligners_paired_and_unpaired():
    stats = lib.clean_paired_fastqs(
        fastqs=[(data_dir / "h37rv_10.r1.fastq.gz", data_dir / "h37rv_10.r2.fastq.gz")],
        aligner=lib.ALIGNER.bowtie2,
        index=data_dir / "MN908947/MN908947",
        out_dir=Path("tst"),
        force=True,
    )
    assert (
        stats[0]["aligner"] == "bowtie2"
        and stats[0]["fastq2_out_name"] == "h37rv_10.r2.clean_2.fastq.gz"
    )

    stats = lib.clean_paired_fastqs(
        fastqs=[(data_dir / "h37rv_10.r1.fastq.gz", data_dir / "h37rv_10.r2.fastq.gz")],
        aligner=lib.ALIGNER.minimap2,
        index=data_dir / "MN908947/MN908947.fasta.gz",
        out_dir=Path("tst"),
        force=True,
    )
    assert (
        stats[0]["aligner"] == "minimap2"
        and stats[0]["fastq2_out_name"] == "h37rv_10.r2.clean_2.fastq.gz"
    )

    stats = lib.clean_fastqs(
        fastqs=[data_dir / "h37rv_10.r1.fastq.gz"],
        aligner=lib.ALIGNER.bowtie2,
        index=data_dir / "MN908947/MN908947",
        out_dir=Path("tst"),
        force=True,
    )
    assert (
        stats[0]["aligner"] == "bowtie2"
        and stats[0]["fastq1_out_name"] == "h37rv_10.r1.clean.fastq.gz"
    )

    stats = lib.clean_fastqs(
        fastqs=[data_dir / "h37rv_10.r1.fastq.gz"],
        aligner=lib.ALIGNER.minimap2,
        index=data_dir / "MN908947/MN908947.fasta.gz",
        out_dir=Path("tst"),
        force=True,
    )
    assert (
        stats[0]["aligner"] == "minimap2"
        and stats[0]["fastq1_out_name"] == "h37rv_10.r1.clean.fastq.gz"
    )

    shutil.rmtree(Path("tst"))


def test_rename():
    stats = lib.clean_paired_fastqs(
        fastqs=[
            (
                data_dir / "tuberculosis_1_2.fastq.gz",
                data_dir / "tuberculosis_1_1.fastq.gz",
            )
        ],
        aligner=lib.ALIGNER.bowtie2,
        index=data_dir / "MN908947/MN908947",
        rename=True,
        out_dir=Path("tst"),
    )
    first_line = get_first_line_of_gzip_file(
        (Path("tst") / "tuberculosis_1_2.clean_1.fastq.gz").resolve()
    )
    assert first_line == "@1 /1"
    shutil.rmtree(Path("tst"))


def test_with_and_without_force():
    stats = lib.clean_paired_fastqs(
        fastqs=[
            (
                data_dir / "tuberculosis_1_2.fastq.gz",
                data_dir / "tuberculosis_1_1.fastq.gz",
            )
        ],
        aligner=lib.ALIGNER.bowtie2,
        index=data_dir / "MN908947/MN908947",
        rename=True,
        out_dir=Path("tst"),
    )
    with pytest.raises(FileExistsError):
        stats = lib.clean_paired_fastqs(
            fastqs=[
                (
                    data_dir / "tuberculosis_1_2.fastq.gz",
                    data_dir / "tuberculosis_1_1.fastq.gz",
                )
            ],
            aligner=lib.ALIGNER.bowtie2,
            index=data_dir / "MN908947/MN908947",
            rename=True,
            out_dir=Path("tst"),
        )
    shutil.rmtree(Path("tst"))


def test_no_rename():
    stats = lib.clean_paired_fastqs(
        fastqs=[
            (
                data_dir / "tuberculosis_1_2.fastq.gz",
                data_dir / "tuberculosis_1_1.fastq.gz",
            )
        ],
        aligner=lib.ALIGNER.bowtie2,
        index=data_dir / "MN908947/MN908947",
        out_dir=Path("tst"),
    )
    first_line = get_first_line_of_gzip_file(
        (Path("tst") / "tuberculosis_1_2.clean_1.fastq.gz").resolve()
    )
    assert first_line == "@NC_000962.3_3000195_3000563_0_1_0_0_1:0:0_0:0:0_0/1"
    shutil.rmtree(Path("tst"))


def test_broken_fastq_path():
    with pytest.raises(FileNotFoundError):
        stats = lib.clean_fastqs(
            fastqs=[Path("invalid_path.fastq.gz")],
            aligner=lib.ALIGNER.bowtie2,
            index=data_dir / "MN908947/MN908947",
        )


def test_mask():
    lib.mask(
        reference=data_dir / "MN908947/MN908947.fasta.gz",
        target=data_dir / "MN908947/partial-for-mask-testing.fa.gz",
    )
    assert Path("masked/mask.bed").exists() and Path("masked/masked.fa").exists()
    shutil.rmtree("masked")
