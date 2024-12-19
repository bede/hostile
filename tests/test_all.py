import gzip
import os
import shutil
import subprocess
from pathlib import Path


import pytest

from hostile import aligner, lib, util

data_dir = Path("tests/data")


# Unit


def test_get_mmi_path():
    assert aligner.get_mmi_path("/path/to/human-t2t-hla.fa.gz") == Path(
        "/path/to/human-t2t-hla.mmi"
    )
    assert aligner.get_mmi_path("/path/to/human-t2t-hla.mmi") == Path(
        "/path/to/human-t2t-hla.mmi"
    )


def test_fastq_path_to_stem():
    assert util.fastq_path_to_stem("long.fastq.gz") == "long"
    assert util.fastq_path_to_stem("-") == "stdin"


# System


def run(cmd: str, cwd: Path = Path(), env: dict | None = None):
    return subprocess.run(
        cmd,
        cwd=cwd,
        shell=True,
        check=True,
        text=True,
        capture_output=True,
        env=env,
    )


def get_nth_line_of_gzip_file(file_path, line_number=1):
    with gzip.open(file_path, "rt") as fh:
        for i, line in enumerate(fh, start=1):
            if i == line_number:
                return line.strip()
    return None


def test_version_cli():
    run("hostile --version")


def test_minimal_fastq(tmp_path):
    stats = lib.clean_fastqs(
        fastqs=[data_dir / "tuberculosis_1_1.fastq.gz"],
        aligner=lib.ALIGNER.bowtie2,
        index=data_dir / "sars-cov-2/sars-cov-2",
        output=tmp_path,
        force=True,
    )
    assert "rename" not in stats[0]["options"]
    assert stats[0]["reads_out"] == 1


def test_multiple_fastqs_bowtie2(tmp_path):
    stats = lib.clean_fastqs(
        fastqs=[
            data_dir / "sars-cov-2_1_1.fastq",
            data_dir / "human_1_1.fastq.gz",
            data_dir / "tuberculosis_1_1.fastq",
        ],
        aligner=lib.ALIGNER.bowtie2,
        index=data_dir / "sars-cov-2/sars-cov-2",
        output=tmp_path,
        force=True,
    )
    assert stats[0]["reads_out"] == 0
    assert stats[1]["reads_out"] == 1
    assert stats[2]["reads_out"] == 1


def test_multiple_fastqs_minimap2(tmp_path):
    stats = lib.clean_fastqs(
        fastqs=[
            data_dir / "sars-cov-2_1_1.fastq",
            data_dir / "human_1_1.fastq.gz",
            data_dir / "tuberculosis_1_1.fastq",
        ],
        aligner=lib.ALIGNER.minimap2,
        index=data_dir / "sars-cov-2/sars-cov-2.fasta.gz",
        output=tmp_path,
        force=True,
    )
    assert stats[0]["reads_out"] == 0
    assert stats[1]["reads_out"] == 1
    assert stats[2]["reads_out"] == 1


def test_multiple_paired_fastqs_bowtie2(tmp_path):
    stats = lib.clean_paired_fastqs(
        fastqs=[
            (data_dir / "sars-cov-2_1_1.fastq", data_dir / "sars-cov-2_1_2.fastq"),
            (data_dir / "human_1_1.fastq.gz", data_dir / "human_1_2.fastq.gz"),
            (data_dir / "tuberculosis_1_1.fastq", data_dir / "tuberculosis_1_2.fastq"),
        ],
        aligner=lib.ALIGNER.bowtie2,
        index=data_dir / "sars-cov-2/sars-cov-2",
        output=tmp_path,
        force=True,
    )
    assert "rename" not in stats[0]["options"]
    assert stats[0]["reads_out"] == 0
    assert stats[1]["reads_out"] == 2
    assert stats[2]["reads_out"] == 2


def test_multiple_paired_fastqs_minimap2(tmp_path):
    stats = lib.clean_paired_fastqs(
        fastqs=[
            (data_dir / "sars-cov-2_1_1.fastq", data_dir / "sars-cov-2_1_2.fastq"),
            (data_dir / "human_1_1.fastq.gz", data_dir / "human_1_2.fastq.gz"),
            (data_dir / "tuberculosis_1_1.fastq", data_dir / "tuberculosis_1_2.fastq"),
        ],
        aligner=lib.ALIGNER.minimap2,
        index=data_dir / "sars-cov-2/sars-cov-2.fasta.gz",
        output=tmp_path,
        force=True,
    )
    assert stats[0]["reads_out"] == 0
    assert stats[1]["reads_out"] == 2
    assert stats[2]["reads_out"] == 2


def test_minimal_paired_fastqs(tmp_path):
    stats = lib.clean_paired_fastqs(
        fastqs=[
            (
                data_dir / "tuberculosis_1_1.fastq.gz",
                data_dir / "tuberculosis_1_2.fastq.gz",
            )
        ],
        aligner=lib.ALIGNER.bowtie2,
        index=data_dir / "sars-cov-2/sars-cov-2",
        output=tmp_path,
        force=True,
    )
    assert stats[0]["reads_out"] == 2


def test_minimal_uncompressed_paired_fastqs(tmp_path):
    lib.clean_paired_fastqs(
        fastqs=[
            (data_dir / "tuberculosis_1_1.fastq", data_dir / "tuberculosis_1_2.fastq")
        ],
        aligner=lib.ALIGNER.bowtie2,
        index=data_dir / "sars-cov-2/sars-cov-2",
        output=tmp_path,
        force=True,
    )


def test_minimal_paired_fastqs_cli(tmp_path):
    run(
        f"hostile clean --index {data_dir}/sars-cov-2/sars-cov-2 --fastq1 {data_dir}/tuberculosis_1_1.fastq.gz --fastq2 {data_dir}/tuberculosis_1_2.fastq.gz --output {tmp_path} --force"
    )


def test_custom_index(tmp_path):
    lib.clean_paired_fastqs(
        fastqs=[
            (
                data_dir / "tuberculosis_1_1.fastq.gz",
                data_dir / "tuberculosis_1_2.fastq.gz",
            )
        ],
        index=data_dir / "sars-cov-2/sars-cov-2",
        output=tmp_path,
        force=True,
    )


def test_both_aligners_paired_and_unpaired(tmp_path):
    stats = lib.clean_paired_fastqs(
        fastqs=[
            (
                data_dir / "tuberculosis_1_1.fastq.gz",
                data_dir / "tuberculosis_1_2.fastq.gz",
            )
        ],
        aligner=lib.ALIGNER.bowtie2,
        index=data_dir / "sars-cov-2/sars-cov-2",
        output=tmp_path,
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
        output=tmp_path,
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
        output=tmp_path,
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
        output=tmp_path,
        force=True,
    )
    assert (
        stats[0]["aligner"] == "minimap2"
        and stats[0]["fastq1_out_name"] == "tuberculosis_1_1.clean.fastq.gz"
    )


def test_rename(tmp_path):
    stats = lib.clean_fastqs(
        fastqs=[data_dir / "tuberculosis_1_1.fastq.gz"],
        aligner=lib.ALIGNER.bowtie2,
        index=data_dir / "sars-cov-2/sars-cov-2",
        rename=True,
        output=tmp_path,
    )
    first_line = get_nth_line_of_gzip_file(tmp_path / "tuberculosis_1_1.clean.fastq.gz")
    assert first_line == "@1"
    assert "rename" in stats[0]["options"]


def test_rename_two_records(tmp_path):
    stats = lib.clean_fastqs(
        fastqs=[data_dir / "tuberculosis_1_12.fastq"],
        aligner=lib.ALIGNER.bowtie2,
        index=data_dir / "sars-cov-2/sars-cov-2",
        rename=True,
        output=tmp_path,
    )
    first_line = get_nth_line_of_gzip_file(
        tmp_path / "tuberculosis_1_12.clean.fastq.gz"
    )
    fifth_line = get_nth_line_of_gzip_file(
        tmp_path / "tuberculosis_1_12.clean.fastq.gz", line_number=5
    )
    assert first_line == "@1"
    assert fifth_line == "@2"
    assert "rename" in stats[0]["options"]


def test_paired_rename(tmp_path):
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
        output=tmp_path,
    )
    first_line = get_nth_line_of_gzip_file(
        tmp_path / "tuberculosis_1_2.clean_1.fastq.gz"
    )
    assert first_line == "@1/1"
    assert "rename" in stats[0]["options"]


def test_with_and_without_force(tmp_path):
    lib.clean_paired_fastqs(
        fastqs=[
            (
                data_dir / "tuberculosis_1_2.fastq.gz",
                data_dir / "tuberculosis_1_1.fastq.gz",
            )
        ],
        aligner=lib.ALIGNER.bowtie2,
        index=data_dir / "sars-cov-2/sars-cov-2",
        rename=True,
        output=tmp_path,
    )
    with pytest.raises(FileExistsError):
        lib.clean_paired_fastqs(
            fastqs=[
                (
                    data_dir / "tuberculosis_1_2.fastq.gz",
                    data_dir / "tuberculosis_1_1.fastq.gz",
                )
            ],
            aligner=lib.ALIGNER.bowtie2,
            index=data_dir / "sars-cov-2/sars-cov-2",
            rename=True,
            output=tmp_path,
        )


def test_no_rename(tmp_path):
    lib.clean_paired_fastqs(
        fastqs=[
            (
                data_dir / "tuberculosis_1_2.fastq.gz",
                data_dir / "tuberculosis_1_1.fastq.gz",
            )
        ],
        aligner=lib.ALIGNER.bowtie2,
        index=data_dir / "sars-cov-2/sars-cov-2",
        output=tmp_path,
        force=True,
    )
    first_line = get_nth_line_of_gzip_file(
        (tmp_path / "tuberculosis_1_2.clean_1.fastq.gz").absolute()
    )
    assert first_line == "@Mycobacterium_tuberculosis/1"


def test_broken_fastq_path(tmp_path):
    with pytest.raises(FileNotFoundError):
        lib.clean_fastqs(
            fastqs=[Path("invalid_path.fastq.gz")],
            aligner=lib.ALIGNER.bowtie2,
            index=data_dir / "sars-cov-2/sars-cov-2",
            output=tmp_path,
        )


def test_no_reads_remaining_after_decontamination(tmp_path):
    stats = lib.clean_fastqs(
        fastqs=[
            data_dir / "sars-cov-2_1_1.fastq",
        ],
        index=data_dir / "sars-cov-2/sars-cov-2.fasta.gz",
        output=tmp_path,
        force=True,
    )
    assert stats[0]["reads_out"] == 0


def test_no_reads_remaining_after_decontamination_paired(tmp_path):
    stats = lib.clean_paired_fastqs(
        fastqs=[
            (
                data_dir / "sars-cov-2_1_1.fastq",
                data_dir / "sars-cov-2_1_2.fastq",
            )
        ],
        index=data_dir / "sars-cov-2/sars-cov-2",
        output=tmp_path,
        force=True,
    )
    assert stats[0]["reads_out"] == 0


def test_decontamination_performance_sars2_bowtie2(tmp_path):
    stats = lib.clean_paired_fastqs(
        fastqs=[
            (
                data_dir / "sars-cov-2_100_1.fastq.gz",
                data_dir / "sars-cov-2_100_2.fastq.gz",
            )
        ],
        index=data_dir / "sars-cov-2/sars-cov-2",
        output=tmp_path,
    )
    assert stats[0]["reads_out"] == 6


def test_decontamination_performance_sars2_minimap2(tmp_path):
    stats = lib.clean_paired_fastqs(
        fastqs=[
            (
                data_dir / "sars-cov-2_100_1.fastq.gz",
                data_dir / "sars-cov-2_100_2.fastq.gz",
            )
        ],
        index=data_dir / "sars-cov-2/sars-cov-2.fasta.gz",
        aligner=lib.ALIGNER.minimap2,
        output=tmp_path,
    )
    assert stats[0]["reads_out"] == 0


def test_mask(tmp_path):
    masked_ref_path, _, _ = lib.mask(
        reference=data_dir / "sars-cov-2/sars-cov-2.fasta.gz",
        target=data_dir / "sars-cov-2/partial-for-mask-testing.fa.gz",
        output=tmp_path,
    )
    assert (Path(tmp_path) / "mask.bed").exists()
    assert (Path(tmp_path) / "mask.bed").exists() and masked_ref_path.exists()


def test_mask_performance(tmp_path):
    masked_ref_path, masked_ref_index_path, n_masked_positions = lib.mask(
        reference=data_dir / "mask/t2t-chm13v2.0-chr21-subset.fa",
        target=data_dir / "mask/gallid-herpesvirus-2.fa",
    )
    assert Path("masked/masked.fa").exists()
    assert str(masked_ref_path) == "masked/masked.fa"
    assert str(masked_ref_index_path) == "masked/masked"
    assert n_masked_positions == 2255
    shutil.rmtree("masked")


def test_sort(tmp_path):
    stats = lib.clean_fastqs(
        fastqs=[data_dir / "sars-cov-2_100_1.fastq.gz"],
        aligner=lib.ALIGNER.minimap2,
        index=data_dir / "sars-cov-2/sars-cov-2.fasta.gz",
        reorder=True,
        output=tmp_path,
        force=True,
    )
    first_line_1 = get_nth_line_of_gzip_file(
        tmp_path / "sars-cov-2_100_1.clean.fastq.gz", line_number=1
    )
    assert stats[0]["reads_out"] == 1
    assert first_line_1 == "@NB552678:8:HGKJNAFX3:1:11104:3356:2796"


def test_sort_rename(tmp_path):
    stats = lib.clean_fastqs(
        fastqs=[data_dir / "sars-cov-2_100_1.fastq.gz"],
        aligner=lib.ALIGNER.minimap2,
        index=data_dir / "sars-cov-2/sars-cov-2.fasta.gz",
        rename=True,
        reorder=True,
        output=tmp_path,
        force=True,
    )
    first_line_1 = get_nth_line_of_gzip_file(
        tmp_path / "sars-cov-2_100_1.clean.fastq.gz", line_number=1
    )
    assert stats[0]["reads_out"] == 1
    assert first_line_1 == "@1"


def test_paired_sort(tmp_path):
    stats = lib.clean_paired_fastqs(
        fastqs=[
            (
                data_dir / "sars-cov-2_100_1.fastq.gz",
                data_dir / "sars-cov-2_100_2.fastq.gz",
            ),
        ],
        aligner=lib.ALIGNER.bowtie2,
        index=data_dir / "sars-cov-2/sars-cov-2",
        reorder=True,
        output=tmp_path,
        force=True,
    )
    first_line_1 = get_nth_line_of_gzip_file(
        tmp_path / "sars-cov-2_100_1.clean_1.fastq.gz", line_number=1
    )
    first_line_2 = get_nth_line_of_gzip_file(
        tmp_path / "sars-cov-2_100_2.clean_2.fastq.gz", line_number=1
    )
    assert stats[0]["reads_out"] == 6
    assert first_line_1 == "@NB552678:8:HGKJNAFX3:1:11101:21702:4929/1"
    assert first_line_2 == "@NB552678:8:HGKJNAFX3:1:11101:21702:4929/2"


def test_paired_sort_rename(tmp_path):
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
        reorder=True,
        output=tmp_path,
        force=True,
    )
    first_line_1 = get_nth_line_of_gzip_file(
        tmp_path / "sars-cov-2_100_1.clean_1.fastq.gz", line_number=1
    )
    first_line_2 = get_nth_line_of_gzip_file(
        tmp_path / "sars-cov-2_100_2.clean_2.fastq.gz", line_number=1
    )
    fifth_line_1 = get_nth_line_of_gzip_file(
        tmp_path / "sars-cov-2_100_1.clean_1.fastq.gz", line_number=5
    )
    fifth_line_2 = get_nth_line_of_gzip_file(
        tmp_path / "sars-cov-2_100_2.clean_2.fastq.gz", line_number=5
    )
    assert stats[0]["reads_out"] == 6
    assert first_line_1 == "@1/1"
    assert first_line_2 == "@1/2"
    assert fifth_line_1 == "@2/1"
    assert fifth_line_2 == "@2/2"


def test_minimap2_aligner_args(tmp_path):
    stats = lib.clean_fastqs(
        fastqs=[data_dir / "sars-cov-2_100_1.fastq.gz"],
        aligner=lib.ALIGNER.minimap2,
        index=data_dir / "sars-cov-2/sars-cov-2.fasta.gz",
        reorder=True,
        output=tmp_path,
        aligner_args="-x asm5",  # Lets everything through
        force=True,
    )
    assert stats[0]["reads_out"] == 50


def test_bowtie2_aligner_args(tmp_path):
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
        reorder=True,
        output=tmp_path,
        aligner_args="--ignore-quals",
        force=True,
    )
    assert stats[0]["reads_out"] == 8


def test_invert_single(tmp_path):
    stats = lib.clean_fastqs(
        fastqs=[data_dir / "sars-cov-2_100_1.fastq.gz"],
        aligner=lib.ALIGNER.minimap2,
        index=data_dir / "sars-cov-2/sars-cov-2.fasta.gz",
        output=tmp_path,
        force=True,
    )
    assert stats[0]["reads_out"] == 1

    stats_invert = lib.clean_fastqs(
        fastqs=[data_dir / "sars-cov-2_100_1.fastq.gz"],
        aligner=lib.ALIGNER.minimap2,
        index=data_dir / "sars-cov-2/sars-cov-2.fasta.gz",
        output=tmp_path,
        invert=True,
        force=True,
    )
    assert stats_invert[0]["reads_in"] - stats_invert[0]["reads_out"] == 1


def test_invert_paired(tmp_path):
    stats = lib.clean_paired_fastqs(
        fastqs=[
            (
                data_dir / "sars-cov-2_100_1.fastq.gz",
                data_dir / "sars-cov-2_100_2.fastq.gz",
            ),
        ],
        aligner=lib.ALIGNER.bowtie2,
        index=data_dir / "sars-cov-2/sars-cov-2",
        output=tmp_path,
        force=True,
    )
    assert stats[0]["reads_out"] == 6

    stats_invert = lib.clean_paired_fastqs(
        fastqs=[
            (
                data_dir / "sars-cov-2_100_1.fastq.gz",
                data_dir / "sars-cov-2_100_2.fastq.gz",
            ),
        ],
        aligner=lib.ALIGNER.bowtie2,
        index=data_dir / "sars-cov-2/sars-cov-2",
        output=tmp_path,
        invert=True,
        force=True,
    )
    assert stats_invert[0]["reads_in"] - stats_invert[0]["reads_out"] == 6


def test_minimap2_reordering_linux(tmp_path):
    lib.clean_paired_fastqs(
        fastqs=[
            (
                data_dir / "sars-cov-2_1_1.fastq",
                data_dir / "sars-cov-2_1_2.fastq",
            ),
        ],
        aligner=lib.ALIGNER.minimap2,
        index=data_dir / "sars-cov-2/sars-cov-2.fasta.gz",
        reorder=True,
        output=tmp_path,
        force=True,
    )


def test_rename_invert_single(tmp_path):
    stats = lib.clean_fastqs(
        fastqs=[data_dir / "sars-cov-2_1_1.fastq"],
        aligner=lib.ALIGNER.minimap2,
        index=data_dir / "sars-cov-2/sars-cov-2.fasta.gz",
        output=tmp_path,
        invert=True,
        rename=True,
        force=True,
    )
    assert stats[0]["reads_out"] == 1


def test_rename_invert_paired(tmp_path):
    stats = lib.clean_paired_fastqs(
        fastqs=[
            (
                data_dir / "sars-cov-2_1_1.fastq",
                data_dir / "sars-cov-2_1_2.fastq",
            ),
        ],
        aligner=lib.ALIGNER.bowtie2,
        index=data_dir / "sars-cov-2/sars-cov-2",
        output=tmp_path,
        invert=True,
        rename=True,
        force=True,
    )
    assert stats[0]["reads_out"] == 2


def test_stats_options():
    stats = lib.clean_fastqs(
        fastqs=[data_dir / "sars-cov-2_1_1.fastq"],
        aligner=lib.ALIGNER.minimap2,
        index=data_dir / "sars-cov-2/sars-cov-2.fasta.gz",
        output="-",
        invert=True,
        rename=True,
        reorder=True,
        force=True,
    )
    assert {"rename", "reorder", "invert", "stdout"} == set(stats[0]["options"])


def test_fixing_empty_fastqs_single(tmp_path):
    stats = lib.clean_fastqs(
        fastqs=[data_dir / "sars-cov-2_1_1.fastq"],
        aligner=lib.ALIGNER.minimap2,
        index=data_dir / "sars-cov-2/sars-cov-2.fasta.gz",
        output=tmp_path,
        force=True,
    )
    run(f"gzip -dc {stats[0]['fastq1_out_path']}")


def test_fixing_empty_fastqs_paired(tmp_path):
    stats = lib.clean_paired_fastqs(
        fastqs=[
            (
                data_dir / "sars-cov-2_1_1.fastq",
                data_dir / "sars-cov-2_1_2.fastq",
            ),
        ],
        aligner=lib.ALIGNER.bowtie2,
        index=data_dir / "sars-cov-2/sars-cov-2",
        output=tmp_path,
        force=True,
    )
    run(f"gzip -dc {stats[0]['fastq1_out_path']}")
    run(f"gzip -dc {stats[0]['fastq2_out_path']}")


def test_mismatched_number_of_reads_bowtie2(tmp_path):
    """This has caused sinister errors in the wild, yet is handled gracefully here"""
    with pytest.raises(RuntimeError):
        lib.clean_paired_fastqs(
            fastqs=[
                (
                    data_dir / "sars-cov-2_100_1.fastq.gz",
                    data_dir / "sars-cov-2_1_2.fastq",
                ),
            ],
            aligner=lib.ALIGNER.bowtie2,
            index=data_dir / "sars-cov-2/sars-cov-2",
            output="-",
        )


def test_airplane_invalid_mm2_standard_index_name(tmp_path):
    with pytest.raises(FileNotFoundError):
        lib.clean_fastqs(
            fastqs=[data_dir / "sars-cov-2_1_1.fastq"],
            index="invalid_index_name",
            aligner=lib.ALIGNER.minimap2,
            output=tmp_path,
            airplane=True,
        )


def test_airplane_invalid_bt2_standard_index_name(tmp_path):
    with pytest.raises(FileNotFoundError):
        lib.clean_fastqs(
            fastqs=[data_dir / "sars-cov-2_1_1.fastq"],
            index="invalid_index_name",
            aligner=lib.ALIGNER.bowtie2,
            output=tmp_path,
            airplane=True,
        )


def test_override_cache_dir(tmp_path):
    env = os.environ.copy()
    env["HOSTILE_CACHE_DIR"] = "custom_directory"
    result = run(
        f"hostile clean --debug --aligner minimap2 --index {data_dir}/sars-cov-2/sars-cov-2.fasta.gz --fastq1 {data_dir}/tuberculosis_1_1.fastq.gz --output {tmp_path} --force",
        env=env,
    )
    assert "custom_directory" in result.stderr
    shutil.rmtree("custom_directory")


def test_override_cache_dir_paired(tmp_path):
    env = os.environ.copy()
    env["HOSTILE_CACHE_DIR"] = "custom_directory"
    result = run(
        f"hostile clean --debug --aligner bowtie2 --index {data_dir}/sars-cov-2/sars-cov-2 --fastq1 {data_dir}/tuberculosis_1_1.fastq.gz --fastq2 {data_dir}/tuberculosis_1_2.fastq.gz --output {tmp_path} --force",
        env=env,
    )
    assert "custom_directory" in result.stderr
    shutil.rmtree("custom_directory")


def test_override_repository_url(tmp_path):
    env = os.environ.copy()
    env["HOSTILE_REPOSITORY_URL"] = "http://example.com"
    result = run(
        f"hostile clean --debug --aligner minimap2 --index {data_dir}/sars-cov-2/sars-cov-2.fasta.gz --fastq1 {data_dir}/tuberculosis_1_1.fastq.gz --output {tmp_path} --force",
        env=env,
    )
    assert "http://example.com" in result.stderr


def test_override_repository_url_paired(tmp_path):
    env = os.environ.copy()
    env["HOSTILE_REPOSITORY_URL"] = "http://example.com"
    result = run(
        f"hostile clean --debug --aligner bowtie2 --index {data_dir}/sars-cov-2/sars-cov-2 --fastq1 {data_dir}/tuberculosis_1_1.fastq.gz --fastq2 {data_dir}/tuberculosis_1_2.fastq.gz --output {tmp_path} --force",
        env=env,
    )
    assert "http://example.com" in result.stderr


def test_stdout_single_bt2():
    result = run(
        f"hostile clean --aligner bowtie2 --index {data_dir}/sars-cov-2/sars-cov-2 --fastq1 {data_dir}/tuberculosis_1_1.fastq --force -o -"
    )
    assert "@Mycobacterium_tuberculosis" in result.stdout
    assert result.stdout.count("\n") == 4


def test_stdout_single_mm2():
    result = run(
        f"hostile clean --aligner minimap2 --index {data_dir}/sars-cov-2/sars-cov-2.fasta.gz --fastq1 {data_dir}/tuberculosis_1_1.fastq --force -o -"
    )
    assert "@Mycobacterium_tuberculosis" in result.stdout
    assert result.stdout.count("\n") == 4


def test_stdout_paired_bt2():
    result = run(
        f"hostile clean --aligner bowtie2 --index {data_dir}/sars-cov-2/sars-cov-2 --fastq1 {data_dir}/tuberculosis_1_1.fastq --fastq2 {data_dir}/tuberculosis_1_2.fastq --force -o -"
    )
    assert "@Mycobacterium_tuberculosis/1" in result.stdout
    assert "@Mycobacterium_tuberculosis/2" in result.stdout
    assert result.stdout.count("\n") == 8


def test_stdout_paired_mm2():
    result = run(
        f"hostile clean --aligner minimap2 --index {data_dir}/sars-cov-2/sars-cov-2.fasta.gz --fastq1 {data_dir}/tuberculosis_1_1.fastq --fastq2 {data_dir}/tuberculosis_1_2.fastq --force -o -"
    )
    assert "@Mycobacterium_tuberculosis/1" in result.stdout
    assert "@Mycobacterium_tuberculosis/2" in result.stdout
    assert result.stdout.count("\n") == 8


def test_log_keys(tmp_path):
    stats = lib.clean_fastqs(
        fastqs=[data_dir / "sars-cov-2_100_1.fastq.gz"],
        aligner=lib.ALIGNER.minimap2,
        index=data_dir / "sars-cov-2/sars-cov-2.fasta.gz",
        output=tmp_path,
    )
    stats_paired = lib.clean_paired_fastqs(
        fastqs=[
            (
                data_dir / "sars-cov-2_1_1.fastq",
                data_dir / "sars-cov-2_1_2.fastq",
            ),
        ],
        aligner=lib.ALIGNER.bowtie2,
        index=data_dir / "sars-cov-2/sars-cov-2",
        output=tmp_path,
    )
    assert set(stats[0].keys()) == {
        "reads_in",
        "reads_out",
        "index",
        "options",
        "aligner",
        "fastq1_in_name",
        "reads_removed",
        "reads_removed_proportion",
        "version",
        "fastq1_in_path",
        "fastq1_out_name",
        "fastq1_out_path",
    }
    assert set(stats_paired[0].keys()) == {
        "reads_in",
        "reads_out",
        "index",
        "options",
        "aligner",
        "fastq2_in_name",
        "reads_removed",
        "fastq1_in_name",
        "reads_removed_proportion",
        "version",
        "fastq1_in_path",
        "fastq2_in_path",
        "fastq1_out_name",
        "fastq2_out_name",
        "fastq1_out_path",
        "fastq2_out_path",
    }


def test_log_keys_stdout():
    stats = lib.clean_fastqs(
        fastqs=[data_dir / "sars-cov-2_100_1.fastq.gz"],
        aligner=lib.ALIGNER.minimap2,
        index=data_dir / "sars-cov-2/sars-cov-2.fasta.gz",
        output="-",
    )
    stats_paired = lib.clean_paired_fastqs(
        fastqs=[
            (
                data_dir / "sars-cov-2_1_1.fastq",
                data_dir / "sars-cov-2_1_2.fastq",
            ),
        ],
        aligner=lib.ALIGNER.bowtie2,
        index=data_dir / "sars-cov-2/sars-cov-2",
        output="-",
    )
    assert set(stats[0].keys()) == {
        "reads_in",
        "reads_out",
        "index",
        "options",
        "aligner",
        "fastq1_in_name",
        "reads_removed",
        "reads_removed_proportion",
        "version",
        "fastq1_in_path",
    }
    assert set(stats_paired[0].keys()) == {
        "reads_in",
        "reads_out",
        "index",
        "options",
        "aligner",
        "fastq2_in_name",
        "reads_removed",
        "fastq1_in_name",
        "reads_removed_proportion",
        "version",
        "fastq1_in_path",
        "fastq2_in_path",
    }


def test_casava_single():
    run_cmd = run(
        f"hostile clean --aligner bowtie2 --index {data_dir}/sars-cov-2/sars-cov-2 --fastq1 {data_dir}/tuberculosis_1_1.fastq -o - --casava"
    )
    stdout_lines = run_cmd.stdout.split("\n")
    assert stdout_lines[0] == "@Mycobacterium_tuberculosis 0:N:0:0"


def test_casava_single_rename():
    run_cmd = run(
        f"hostile clean --aligner bowtie2 --index {data_dir}/sars-cov-2/sars-cov-2 --fastq1 {data_dir}/tuberculosis_1_1.fastq -o - -c --rename"
    )
    stdout_lines = run_cmd.stdout.split("\n")
    assert stdout_lines[0] == "@1 0:N:0:0"


def test_casava_paired():
    run_cmd = run(
        f"hostile clean --index {data_dir}/sars-cov-2/sars-cov-2 --fastq1 {data_dir}/tuberculosis_1_1.fastq --fastq2 {data_dir}/tuberculosis_1_2.fastq -o - --casava"
    )
    stdout_lines = run_cmd.stdout.split("\n")
    assert stdout_lines[0] == "@Mycobacterium_tuberculosis 1:N:0:0"
    assert stdout_lines[4] == "@Mycobacterium_tuberculosis 2:N:0:0"


def test_casava_paired_rename():
    run_cmd = run(
        f"hostile clean --index {data_dir}/sars-cov-2/sars-cov-2 --fastq1 {data_dir}/tuberculosis_1_1.fastq --fastq2 {data_dir}/tuberculosis_1_2.fastq -o - -c --rename"
    )
    stdout_lines = run_cmd.stdout.split("\n")
    assert stdout_lines[0] == "@1 1:N:0:0"
    assert stdout_lines[4] == "@1 2:N:0:0"


def test_single_minimap2_stdout():
    run_cmd = run(
        f"hostile clean --aligner minimap2 --index {data_dir}/sars-cov-2/sars-cov-2.fasta.gz --fastq1 {data_dir}/tuberculosis_1_1.fastq -o -"
    )
    stdout_lines = run_cmd.stdout.split("\n")
    assert stdout_lines[0] == "@Mycobacterium_tuberculosis"
    assert len(stdout_lines) == 5


def test_stdin_single_minimap2_stdout():
    run_cmd = run(
        f"cat {data_dir}/tuberculosis_1_1.fastq | hostile clean --aligner minimap2 --index {data_dir}/sars-cov-2/sars-cov-2.fasta.gz --fastq1 - -o -"
    )
    stdout_lines = run_cmd.stdout.split("\n")
    assert stdout_lines[0] == "@Mycobacterium_tuberculosis"
    assert len(stdout_lines) == 5


def test_stdin_single_minimap2(tmp_path):
    run(
        f"cat {data_dir}/tuberculosis_1_1.fastq | hostile clean --aligner minimap2 --index {data_dir}/sars-cov-2/sars-cov-2.fasta.gz --fastq1 - --output {tmp_path}"
    )
    with gzip.open(tmp_path / "stdin.clean.fastq.gz", "rt") as fh:
        contents = fh.read()
    lines = contents.split("\n")
    assert lines[0] == "@Mycobacterium_tuberculosis"
    assert len(lines) == 5


def test_single_bowtie2_stdout():
    run_cmd = run(
        f"hostile clean --aligner bowtie2 --index {data_dir}/sars-cov-2/sars-cov-2 --fastq1 {data_dir}/tuberculosis_1_1.fastq -o -"
    )
    stdout_lines = run_cmd.stdout.split("\n")
    assert stdout_lines[0] == "@Mycobacterium_tuberculosis"
    assert len(stdout_lines) == 5


def test_stdin_single_bowtie2_stdout():
    run_cmd = run(
        f"cat {data_dir}/tuberculosis_1_1.fastq | hostile clean --aligner bowtie2 --index {data_dir}/sars-cov-2/sars-cov-2 --fastq1 - -o -"
    )
    stdout_lines = run_cmd.stdout.split("\n")
    assert stdout_lines[0] == "@Mycobacterium_tuberculosis"
    assert len(stdout_lines) == 5


def test_stdin_single_bowtie2(tmp_path):
    run(
        f"cat {data_dir}/tuberculosis_1_1.fastq | hostile clean --aligner bowtie2 --index {data_dir}/sars-cov-2/sars-cov-2 --fastq1 - --output {tmp_path}"
    )
    with gzip.open(tmp_path / "stdin.clean.fastq.gz", "rt") as fh:
        contents = fh.read()
    lines = contents.split("\n")
    assert lines[0] == "@Mycobacterium_tuberculosis"
    assert len(lines) == 5


def test_paired_interleaved_bowtie2_stdout():
    run_cmd = run(
        f"hostile clean --aligner bowtie2 --index {data_dir}/sars-cov-2/sars-cov-2 --fastq1 {data_dir}/tuberculosis_1_1.fastq --fastq2 {data_dir}/tuberculosis_1_2.fastq -o -"
    )
    stdout_lines = run_cmd.stdout.split("\n")
    assert stdout_lines[0] == "@Mycobacterium_tuberculosis/1"
    assert stdout_lines[4] == "@Mycobacterium_tuberculosis/2"
    assert len(stdout_lines) == 9


def test_stdin_paired_interleaved_bowtie2_stdout():
    run_cmd = run(
        f"cat {data_dir}/tuberculosis_1_12.fastq | hostile clean --aligner bowtie2 --index {data_dir}/sars-cov-2/sars-cov-2 --fastq1 - --fastq2 - -o -"
    )
    stdout_lines = run_cmd.stdout.split("\n")
    assert stdout_lines[0] == "@Mycobacterium_tuberculosis/1"
    assert stdout_lines[4] == "@Mycobacterium_tuberculosis/2"
    assert len(stdout_lines) == 9


def test_stdin_paired_interleaved_bowtie2(tmp_path):
    run(
        f"cat {data_dir}/tuberculosis_1_12.fastq | hostile clean --aligner bowtie2 --index {data_dir}/sars-cov-2/sars-cov-2 --fastq1 - --fastq2 - --output {tmp_path}"
    )
    with gzip.open(tmp_path / "stdin.clean_1.fastq.gz", "rt") as fh:
        contents = fh.read()
    lines = contents.split("\n")
    assert lines[0] == "@Mycobacterium_tuberculosis/1"
    assert len(lines) == 5
    with gzip.open(tmp_path / "stdin.clean_2.fastq.gz", "rt") as fh:
        contents = fh.read()
    lines = contents.split("\n")
    assert lines[0] == "@Mycobacterium_tuberculosis/2"
    assert len(lines) == 5


def test_paired_interleaved_minimap2_stdout():
    run_cmd = run(
        f"hostile clean --aligner minimap2 --index {data_dir}/sars-cov-2/sars-cov-2.fasta.gz --fastq1 {data_dir}/tuberculosis_1_1.fastq --fastq2 {data_dir}/tuberculosis_1_2.fastq -o -"
    )
    stdout_lines = run_cmd.stdout.split("\n")
    assert stdout_lines[0] == "@Mycobacterium_tuberculosis/1"
    assert stdout_lines[4] == "@Mycobacterium_tuberculosis/2"
    assert len(stdout_lines) == 9


def test_stdin_paired_interleaved_minimap2_stdout():
    run_cmd = run(
        f"cat {data_dir}/tuberculosis_1_12.fastq | hostile clean --aligner minimap2 --index {data_dir}/sars-cov-2/sars-cov-2.fasta.gz --fastq1 - --fastq2 - -o -"
    )
    stdout_lines = run_cmd.stdout.split("\n")
    assert stdout_lines[0] == "@Mycobacterium_tuberculosis/1"
    assert stdout_lines[4] == "@Mycobacterium_tuberculosis/2"
    assert len(stdout_lines) == 9


def test_stdin_paired_interleaved_minimap2(tmp_path):
    run(
        f"cat {data_dir}/tuberculosis_1_12.fastq | hostile clean --aligner minimap2 --index {data_dir}/sars-cov-2/sars-cov-2.fasta.gz --fastq1 - --fastq2 - --output {tmp_path}"
    )
    with gzip.open(tmp_path / "stdin.clean_1.fastq.gz", "rt") as fh:
        contents = fh.read()
    lines = contents.split("\n")
    assert lines[0] == "@Mycobacterium_tuberculosis/1"
    assert len(lines) == 5
    with gzip.open(tmp_path / "stdin.clean_2.fastq.gz", "rt") as fh:
        contents = fh.read()
    lines = contents.split("\n")
    assert lines[0] == "@Mycobacterium_tuberculosis/2"
    assert len(lines) == 5
