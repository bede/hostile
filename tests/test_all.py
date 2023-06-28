import shutil
import subprocess
from pathlib import Path

from hostile import lib

data_dir = Path("tests/data")


def run(cmd, cwd=data_dir):  # Helper for CLI testing
    return subprocess.run(
        cmd, cwd=cwd, shell=True, check=True, text=True, capture_output=True
    )


def test_version_cli():
    run("hostile --version")


def test_minimal_fastq():
    lib.clean_fastqs(
        fastqs=[data_dir / "h37rv_10.r1.fastq.gz"],
        out_dir=Path("test_minimal_fastq"),
    )
    shutil.rmtree("test_minimal_fastq")


def test_minimal_paired_fastqs():
    lib.clean_paired_fastqs(
        fastqs=[(data_dir / "h37rv_10.r1.fastq.gz", data_dir / "h37rv_10.r2.fastq.gz")],
        out_dir=Path("test_minimal_fastqs"),
    )
    shutil.rmtree("test_minimal_fastqs")


def test_minimal_uncompressed_paired_fastqs():
    lib.clean_paired_fastqs(
        fastqs=[(data_dir / "h37rv_10.r1.fastq", data_dir / "h37rv_10.r2.fastq")],
        out_dir=Path("test_minimal_fastqs"),
    )
    shutil.rmtree("test_minimal_fastqs")


def test_minimal_paired_fastqs_cli():
    run(
        f"hostile clean --fastq1 h37rv_10.r1.fastq.gz --fastq2 h37rv_10.r2.fastq.gz --out-dir test_minimal_fastqs"
    )
    shutil.rmtree(f"{data_dir}/test_minimal_fastqs")


def test_many_minimal_paired_fastqs_cli():
    run(
        f"hostile clean-many h37rv_10.r1.fastq.gz,h37rv_10.r2.fastq.gz h37rv_10_2.r1.fastq.gz,h37rv_10_2.r2.fastq.gz --out-dir test_minimal_fastqs"
    )
    shutil.rmtree(f"{data_dir}/test_minimal_fastqs")


def test_custom_index():
    lib.clean_paired_fastqs(
        fastqs=[(data_dir / "h37rv_10.r1.fastq.gz", data_dir / "h37rv_10.r2.fastq.gz")],
        custom_index="tests/data/MN908947/MN908947",
        out_dir=Path("test_minimal_fastqs"),
    )
    shutil.rmtree("test_minimal_fastqs")
