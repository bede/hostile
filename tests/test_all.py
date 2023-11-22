import gzip
import shutil
import subprocess
from pathlib import Path

import pytest

from hostile import lib

data_dir = Path("tests/data")
out_dir = Path("test_data")


def run(cmd: str, cwd: Path = Path()):  # Helper for CLI testing
    return subprocess.run(
        cmd, cwd=cwd, shell=True, check=True, text=True, capture_output=True
    )


def get_nth_line_of_gzip_file(file_path, line_number=1):
    with gzip.open(file_path, "rt") as fh:
        for i, line in enumerate(fh, start=1):
            if i == line_number:
                return line.strip()
    return None


def test_version_cli():
    run("hostile --version")


def test_minimal_fastq():
    stats = lib.clean_fastqs(
        fastqs=[data_dir / "tuberculosis_1_1.fastq.gz"],
        aligner=lib.ALIGNER.bowtie2,
        index=data_dir / "sars-cov-2/sars-cov-2",
        out_dir=out_dir,
        force=True,
    )
    assert stats[0]["rename"] == False
    assert stats[0]["reads_out"] == 1
    shutil.rmtree(out_dir, ignore_errors=True)


def test_multiple_fastqs_bowtie2():
    stats = lib.clean_fastqs(
        fastqs=[
            data_dir / "sars-cov-2_1_1.fastq",
            data_dir / "human_1_1.fastq.gz",
            data_dir / "tuberculosis_1_1.fastq",
        ],
        aligner=lib.ALIGNER.bowtie2,
        index=data_dir / "sars-cov-2/sars-cov-2",
        out_dir=out_dir,
        force=True,
    )
    assert stats[0]["reads_out"] == 0
    assert stats[1]["reads_out"] == 1
    assert stats[2]["reads_out"] == 1
    shutil.rmtree(out_dir, ignore_errors=True)


def test_multiple_fastqs_minimap2():
    stats = lib.clean_fastqs(
        fastqs=[
            data_dir / "sars-cov-2_1_1.fastq",
            data_dir / "human_1_1.fastq.gz",
            data_dir / "tuberculosis_1_1.fastq",
        ],
        aligner=lib.ALIGNER.minimap2,
        index=data_dir / "sars-cov-2/sars-cov-2.fasta.gz",
        out_dir=out_dir,
        force=True,
    )
    assert stats[0]["reads_out"] == 0
    assert stats[1]["reads_out"] == 1
    assert stats[2]["reads_out"] == 1
    shutil.rmtree(out_dir, ignore_errors=True)


def test_multiple_paired_fastqs_bowtie2():
    stats = lib.clean_paired_fastqs(
        fastqs=[
            (data_dir / "sars-cov-2_1_1.fastq", data_dir / "sars-cov-2_1_2.fastq"),
            (data_dir / "human_1_1.fastq.gz", data_dir / "human_1_2.fastq.gz"),
            (data_dir / "tuberculosis_1_1.fastq", data_dir / "tuberculosis_1_2.fastq"),
        ],
        aligner=lib.ALIGNER.bowtie2,
        index=data_dir / "sars-cov-2/sars-cov-2",
        out_dir=out_dir,
        force=True,
    )
    assert stats[0]["rename"] == False
    assert stats[0]["reads_out"] == 0
    assert stats[1]["reads_out"] == 2
    assert stats[2]["reads_out"] == 2
    shutil.rmtree(out_dir, ignore_errors=True)


def test_multiple_paired_fastqs_minimap2():
    stats = lib.clean_paired_fastqs(
        fastqs=[
            (data_dir / "sars-cov-2_1_1.fastq", data_dir / "sars-cov-2_1_2.fastq"),
            (data_dir / "human_1_1.fastq.gz", data_dir / "human_1_2.fastq.gz"),
            (data_dir / "tuberculosis_1_1.fastq", data_dir / "tuberculosis_1_2.fastq"),
        ],
        aligner=lib.ALIGNER.minimap2,
        index=data_dir / "sars-cov-2/sars-cov-2.fasta.gz",
        out_dir=out_dir,
        force=True,
    )
    assert stats[0]["reads_out"] == 0
    assert stats[1]["reads_out"] == 2
    assert stats[2]["reads_out"] == 2
    shutil.rmtree(out_dir, ignore_errors=True)


def test_minimal_paired_fastqs():
    stats = lib.clean_paired_fastqs(
        fastqs=[
            (
                data_dir / "tuberculosis_1_1.fastq.gz",
                data_dir / "tuberculosis_1_2.fastq.gz",
            )
        ],
        aligner=lib.ALIGNER.bowtie2,
        index=data_dir / "sars-cov-2/sars-cov-2",
        out_dir=out_dir,
        force=True,
    )
    assert stats[0]["reads_out"] == 2
    shutil.rmtree(out_dir, ignore_errors=True)


def test_minimal_uncompressed_paired_fastqs():
    shutil.rmtree(out_dir, ignore_errors=True)
    lib.clean_paired_fastqs(
        fastqs=[
            (data_dir / "tuberculosis_1_1.fastq", data_dir / "tuberculosis_1_2.fastq")
        ],
        aligner=lib.ALIGNER.bowtie2,
        index=data_dir / "sars-cov-2/sars-cov-2",
        out_dir=out_dir,
        force=True,
    )
    shutil.rmtree(out_dir, ignore_errors=True)


def test_minimal_paired_fastqs_cli():
    run(
        f"hostile clean --index {data_dir}/sars-cov-2/sars-cov-2 --fastq1 {data_dir}/tuberculosis_1_1.fastq.gz --fastq2 {data_dir}/tuberculosis_1_2.fastq.gz --out-dir {out_dir} --force"
    )
    shutil.rmtree(out_dir)


def test_custom_index():
    lib.clean_paired_fastqs(
        fastqs=[
            (
                data_dir / "tuberculosis_1_1.fastq.gz",
                data_dir / "tuberculosis_1_2.fastq.gz",
            )
        ],
        index=data_dir / "sars-cov-2/sars-cov-2",
        out_dir=out_dir,
        force=True,
    )
    shutil.rmtree(out_dir, ignore_errors=True)


def test_both_aligners_paired_and_unpaired():
    stats = lib.clean_paired_fastqs(
        fastqs=[
            (
                data_dir / "tuberculosis_1_1.fastq.gz",
                data_dir / "tuberculosis_1_2.fastq.gz",
            )
        ],
        aligner=lib.ALIGNER.bowtie2,
        index=data_dir / "sars-cov-2/sars-cov-2",
        out_dir=out_dir,
        force=True,
    )
    assert (
        stats[0]["aligner"] == "bowtie2"
        and stats[0]["fastq2_out_name"] == "tuberculosis_1_2.clean_2.fastq.gz"
    )

    stats = lib.clean_paired_fastqs(
        fastqs=[
            (
                data_dir / "tuberculosis_1_1.fastq.gz",
                data_dir / "tuberculosis_1_2.fastq.gz",
            )
        ],
        aligner=lib.ALIGNER.minimap2,
        index=data_dir / "sars-cov-2/sars-cov-2.fasta.gz",
        out_dir=out_dir,
        force=True,
    )
    assert (
        stats[0]["aligner"] == "minimap2"
        and stats[0]["fastq2_out_name"] == "tuberculosis_1_2.clean_2.fastq.gz"
    )

    stats = lib.clean_fastqs(
        fastqs=[data_dir / "tuberculosis_1_1.fastq.gz"],
        aligner=lib.ALIGNER.bowtie2,
        index=data_dir / "sars-cov-2/sars-cov-2",
        out_dir=out_dir,
        force=True,
    )
    assert (
        stats[0]["aligner"] == "bowtie2"
        and stats[0]["fastq1_out_name"] == "tuberculosis_1_1.clean.fastq.gz"
    )

    stats = lib.clean_fastqs(
        fastqs=[data_dir / "tuberculosis_1_1.fastq.gz"],
        aligner=lib.ALIGNER.minimap2,
        index=data_dir / "sars-cov-2/sars-cov-2.fasta.gz",
        out_dir=out_dir,
        force=True,
    )
    assert (
        stats[0]["aligner"] == "minimap2"
        and stats[0]["fastq1_out_name"] == "tuberculosis_1_1.clean.fastq.gz"
    )
    shutil.rmtree(out_dir, ignore_errors=True)


def test_rename():
    stats = lib.clean_fastqs(
        fastqs=[data_dir / "tuberculosis_1_1.fastq.gz"],
        aligner=lib.ALIGNER.bowtie2,
        index=data_dir / "sars-cov-2/sars-cov-2",
        rename=True,
        out_dir=out_dir,
    )
    first_line = get_nth_line_of_gzip_file(out_dir / "tuberculosis_1_1.clean.fastq.gz")
    assert first_line == "@1"
    assert stats[0]["rename"] == True
    shutil.rmtree(out_dir, ignore_errors=True)


def test_rename_two_records():
    stats = lib.clean_fastqs(
        fastqs=[data_dir / "tuberculosis_2.fastq"],
        aligner=lib.ALIGNER.bowtie2,
        index=data_dir / "sars-cov-2/sars-cov-2",
        rename=True,
        out_dir=out_dir,
    )
    first_line = get_nth_line_of_gzip_file(out_dir / "tuberculosis_2.clean.fastq.gz")
    fifth_line = get_nth_line_of_gzip_file(
        out_dir / "tuberculosis_2.clean.fastq.gz", line_number=5
    )
    assert first_line == "@1"
    assert fifth_line == "@2"
    assert stats[0]["rename"] == True
    shutil.rmtree(out_dir, ignore_errors=True)


def test_paired_rename():
    stats = lib.clean_paired_fastqs(
        fastqs=[
            (
                data_dir / "tuberculosis_1_2.fastq.gz",
                data_dir / "tuberculosis_1_1.fastq.gz",
            )
        ],
        aligner=lib.ALIGNER.bowtie2,
        index=data_dir / "sars-cov-2/sars-cov-2",
        rename=True,
        out_dir=out_dir,
    )
    first_line = get_nth_line_of_gzip_file(
        out_dir / "tuberculosis_1_2.clean_1.fastq.gz"
    )
    assert first_line == "@1 /1"
    assert stats[0]["rename"] == True
    shutil.rmtree(out_dir, ignore_errors=True)


def test_with_and_without_force():
    shutil.rmtree(out_dir, ignore_errors=True)
    stats = lib.clean_paired_fastqs(
        fastqs=[
            (
                data_dir / "tuberculosis_1_2.fastq.gz",
                data_dir / "tuberculosis_1_1.fastq.gz",
            )
        ],
        aligner=lib.ALIGNER.bowtie2,
        index=data_dir / "sars-cov-2/sars-cov-2",
        rename=True,
        out_dir=out_dir,
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
            index=data_dir / "sars-cov-2/sars-cov-2",
            rename=True,
            out_dir=out_dir,
        )
    shutil.rmtree(out_dir, ignore_errors=True)


def test_no_rename():
    stats = lib.clean_paired_fastqs(
        fastqs=[
            (
                data_dir / "tuberculosis_1_2.fastq.gz",
                data_dir / "tuberculosis_1_1.fastq.gz",
            )
        ],
        aligner=lib.ALIGNER.bowtie2,
        index=data_dir / "sars-cov-2/sars-cov-2",
        out_dir=out_dir,
        force=True,
    )
    first_line = get_nth_line_of_gzip_file(
        (out_dir / "tuberculosis_1_2.clean_1.fastq.gz").resolve()
    )
    assert first_line == "@Mycobacterium_tuberculosis/1"
    shutil.rmtree(out_dir, ignore_errors=True)


def test_broken_fastq_path():
    with pytest.raises(FileNotFoundError):
        stats = lib.clean_fastqs(
            fastqs=[Path("invalid_path.fastq.gz")],
            aligner=lib.ALIGNER.bowtie2,
            index=data_dir / "sars-cov-2/sars-cov-2",
            out_dir=out_dir,
        )
    shutil.rmtree(out_dir, ignore_errors=True)


def test_no_reads_remaining_after_decontamination():
    stats = lib.clean_fastqs(
        fastqs=[
            data_dir / "sars-cov-2_1_1.fastq",
        ],
        index=data_dir / "sars-cov-2/sars-cov-2.fasta.gz",
        out_dir=out_dir,
        force=True,
    )
    assert stats[0]["reads_out"] == 0
    shutil.rmtree(out_dir, ignore_errors=True)


def test_no_reads_remaining_after_decontamination_paired():
    stats = lib.clean_paired_fastqs(
        fastqs=[
            (
                data_dir / "sars-cov-2_1_1.fastq",
                data_dir / "sars-cov-2_1_2.fastq",
            )
        ],
        index=data_dir / "sars-cov-2/sars-cov-2",
        out_dir=out_dir,
        force=True,
    )
    assert stats[0]["reads_out"] == 0
    shutil.rmtree(out_dir, ignore_errors=True)


def test_decontamination_performance_sars2_bowtie2():
    stats = lib.clean_paired_fastqs(
        fastqs=[
            (
                data_dir / "sars-cov-2_100_1.fastq.gz",
                data_dir / "sars-cov-2_100_2.fastq.gz",
            )
        ],
        index=data_dir / "sars-cov-2/sars-cov-2",
        out_dir=out_dir,
    )
    assert stats[0]["reads_out"] == 6
    shutil.rmtree(out_dir, ignore_errors=True)


def test_decontamination_performance_sars2_minimap2():
    stats = lib.clean_paired_fastqs(
        fastqs=[
            (
                data_dir / "sars-cov-2_100_1.fastq.gz",
                data_dir / "sars-cov-2_100_2.fastq.gz",
            )
        ],
        index=data_dir / "sars-cov-2/sars-cov-2.fasta.gz",
        aligner=lib.ALIGNER.minimap2,
        out_dir=out_dir,
    )
    assert stats[0]["reads_out"] == 0
    shutil.rmtree(out_dir, ignore_errors=True)


def test_mask():
    lib.mask(
        reference=data_dir / "sars-cov-2/sars-cov-2.fasta.gz",
        target=data_dir / "sars-cov-2/partial-for-mask-testing.fa.gz",
    )
    assert Path("masked/mask.bed").exists() and Path("masked/masked.fa").exists()
    shutil.rmtree("masked")


def test_sort():
    stats = lib.clean_fastqs(
        fastqs=[data_dir / "sars-cov-2_100_1.fastq.gz"],
        aligner=lib.ALIGNER.minimap2,
        index=data_dir / "sars-cov-2/sars-cov-2.fasta.gz",
        sort_by_name=True,
        out_dir=out_dir,
        force=True,
    )
    first_line_1 = get_nth_line_of_gzip_file(
        out_dir / "sars-cov-2_100_1.clean.fastq.gz", line_number=1
    )
    assert stats[0]["reads_out"] == 1
    assert first_line_1 == "@NB552678:8:HGKJNAFX3:1:11104:3356:2796"
    shutil.rmtree(out_dir, ignore_errors=True)


def test_sort_rename():
    stats = lib.clean_fastqs(
        fastqs=[data_dir / "sars-cov-2_100_1.fastq.gz"],
        aligner=lib.ALIGNER.minimap2,
        index=data_dir / "sars-cov-2/sars-cov-2.fasta.gz",
        rename=True,
        sort_by_name=True,
        out_dir=out_dir,
        force=True,
    )
    first_line_1 = get_nth_line_of_gzip_file(
        out_dir / "sars-cov-2_100_1.clean.fastq.gz", line_number=1
    )
    assert stats[0]["reads_out"] == 1
    assert first_line_1 == "@1"
    shutil.rmtree(out_dir, ignore_errors=True)


def test_paired_sort():
    stats = lib.clean_paired_fastqs(
        fastqs=[
            (
                data_dir / "sars-cov-2_100_1.fastq.gz",
                data_dir / "sars-cov-2_100_2.fastq.gz",
            ),
        ],
        aligner=lib.ALIGNER.bowtie2,
        index=data_dir / "sars-cov-2/sars-cov-2",
        sort_by_name=True,
        out_dir=out_dir,
        force=True,
    )
    first_line_1 = get_nth_line_of_gzip_file(
        out_dir / "sars-cov-2_100_1.clean_1.fastq.gz", line_number=1
    )
    first_line_2 = get_nth_line_of_gzip_file(
        out_dir / "sars-cov-2_100_2.clean_2.fastq.gz", line_number=1
    )
    assert stats[0]["reads_out"] == 6
    assert first_line_1 == "@NB552678:8:HGKJNAFX3:1:11101:21702:4929/1"
    assert first_line_2 == "@NB552678:8:HGKJNAFX3:1:11101:21702:4929/2"
    shutil.rmtree(out_dir, ignore_errors=True)


def test_paired_sort_rename():
    stats = lib.clean_paired_fastqs(
        fastqs=[
            (
                data_dir / "sars-cov-2_100_1.fastq.gz",
                data_dir / "sars-cov-2_100_2.fastq.gz",
            ),
        ],
        aligner=lib.ALIGNER.bowtie2,
        index=data_dir / "sars-cov-2/sars-cov-2",
        rename=True,
        sort_by_name=True,
        out_dir=out_dir,
        force=True,
    )
    first_line_1 = get_nth_line_of_gzip_file(
        out_dir / "sars-cov-2_100_1.clean_1.fastq.gz", line_number=1
    )
    first_line_2 = get_nth_line_of_gzip_file(
        out_dir / "sars-cov-2_100_2.clean_2.fastq.gz", line_number=1
    )
    fifth_line_1 = get_nth_line_of_gzip_file(
        out_dir / "sars-cov-2_100_1.clean_1.fastq.gz", line_number=5
    )
    fifth_line_2 = get_nth_line_of_gzip_file(
        out_dir / "sars-cov-2_100_2.clean_2.fastq.gz", line_number=5
    )
    assert stats[0]["reads_out"] == 6
    assert first_line_1 == "@1 /1"
    assert first_line_2 == "@1 /2"
    assert fifth_line_1 == "@2 /1"
    assert fifth_line_2 == "@2 /2"
    shutil.rmtree(out_dir, ignore_errors=True)


def test_minimap2_aligner_args():
    stats = lib.clean_fastqs(
        fastqs=[data_dir / "sars-cov-2_100_1.fastq.gz"],
        aligner=lib.ALIGNER.minimap2,
        index=data_dir / "sars-cov-2/sars-cov-2.fasta.gz",
        sort_by_name=True,
        out_dir=out_dir,
        aligner_args="-x asm5",  # Lets everything through
        force=True,
    )
    assert stats[0]["reads_out"] == 50
    shutil.rmtree(out_dir, ignore_errors=True)


def test_bowtie2_aligner_args():
    stats = lib.clean_paired_fastqs(
        fastqs=[
            (
                data_dir / "sars-cov-2_100_1.fastq.gz",
                data_dir / "sars-cov-2_100_2.fastq.gz",
            ),
        ],
        aligner=lib.ALIGNER.bowtie2,
        index=data_dir / "sars-cov-2/sars-cov-2",
        rename=True,
        sort_by_name=True,
        out_dir=out_dir,
        aligner_args="--ignore-quals",
        force=True,
    )
    assert stats[0]["reads_out"] == 8
    shutil.rmtree(out_dir, ignore_errors=True)
