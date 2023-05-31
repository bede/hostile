import shutil
import subprocess
from pathlib import Path

from hostile import lib

data_dir = Path("tests/data")


def run(cmd, cwd=data_dir):  # Helper for CLI testing
    return subprocess.run(
        cmd, cwd=cwd, shell=True, check=True, text=True, capture_output=True
    )


def test_cli_version():
    run("hostile --version")


def test_minimal_fastqs():
    lib.decontaminate_paired_method_1(fastq1=data_dir / "h37rv_10.r1.fastq.gz", fastq2=data_dir / "h37rv_10.r1.fastq.gz", ref=data_dir / "human.1k.fa.gz", out_dir="test_minimal_fastqs")
    # shutil.rmtree("test_minimal_fastqs")
