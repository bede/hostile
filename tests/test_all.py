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


def test_minimal_fastqs():
    lib.dehost_fastqs(fastq1=data_dir / "h37rv_10.r1.fastq.gz", fastq2=data_dir / "h37rv_10.r1.fastq.gz", ref=data_dir / "human.1k.fa.gz", out_dir="test_minimal_fastqs")
    shutil.rmtree("test_minimal_fastqs")


def test_minimal_fastqs_cli():
    run(f"hostile --fastq1 h37rv_10.r1.fastq.gz --fastq2 h37rv_10.r2.fastq.gz --ref human.1k.fa.gz --out-dir test_minimal_fastqs")
    shutil.rmtree(f"{data_dir}/test_minimal_fastqs")
