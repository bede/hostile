import logging
import shutil
import subprocess
import tempfile

from dataclasses import dataclass
from enum import Enum
from pathlib import Path

from hostile import util


def get_mmi_path(index_path: Path) -> Path:
    return Path(
        str(index_path)
        .removesuffix(".fa")
        .removesuffix(".fasta")
        .removesuffix(".fa.gz")
        .removesuffix(".fasta.gz")
        .removesuffix(".mmi")
        + ".mmi"
    )


@dataclass
class Aligner:
    name: str
    bin_path: Path
    data_dir: Path
    single_cmd: str
    single_unindexed_cmd: str
    paired_cmd: str
    paired_unindexed_cmd: str
    interleaved_cmd: str
    interleaved_unindexed_cmd: str

    def __post_init__(self):
        Path(self.data_dir).mkdir(exist_ok=True, parents=True)

    def check_index(
        self, index: str, airplane: bool = False, build_mmi: bool = False
    ) -> Path:
        """Test aligner and check/download a ref/index if necessary, returning genome or index path"""
        try:
            util.run(f"{self.bin_path} --version", cwd=self.data_dir)
        except subprocess.CalledProcessError:
            raise RuntimeError(f"Failed to execute {self.bin_path}")

        if self.name == "Bowtie2":
            if Path(f"{index}.1.bt2").is_file():
                index_path = Path(index)
                logging.info(f"Found custom index {index_path}")
            elif (self.data_dir / f"{index}.1.bt2").is_file():
                index_path = self.data_dir / index
                logging.info(f"Found cached standard index {index}")
            elif not airplane and util.fetch_manifest(util.INDEX_REPOSITORY_URL).get(
                index
            ):
                file_name = f"{index}.tar"
                file_url = f"{util.INDEX_REPOSITORY_URL}/{file_name}"
                logging.info(f"Fetching standard index {index} ({file_url})")
                manifest = util.fetch_manifest(util.INDEX_REPOSITORY_URL)
                with tempfile.NamedTemporaryFile() as temporary_file:
                    tmp_path = Path(temporary_file.name)
                    util.download(f"{file_url}", tmp_path)
                    expected_sha256 = manifest[index]["assets"][file_name]["sha256"]
                    logging.info(f"Verifying checksum {expected_sha256}…")
                    observed_sha256 = util.sha256(tmp_path)
                    if observed_sha256 != expected_sha256:
                        raise ValueError(f"Checksum mismatch for {file_name}")
                    logging.info(f"Extracting {file_name}…")
                    util.untar_file(tmp_path, self.data_dir)
                index_path = self.data_dir / index
                logging.info(f"Downloaded standard index {index_path}")
            else:
                message = f"{index} is neither a valid custom Bowtie2 index path nor a valid standard index name. Mode: short read (Bowtie2)"
                if airplane:
                    message += ". Disable airplane mode to enable discovery of standard indexes"
                raise FileNotFoundError(message)

        elif self.name == "Minimap2":
            if Path(f"{index}").is_file():
                index_path = Path(index)
                if get_mmi_path(index_path).is_file():
                    logging.info(f"Found cached standard index {index} (MMI available)")
                else:
                    logging.info(
                        f"Found cached standard index {index} (building MMI cache; do not interrupt)"
                    )
            elif (self.data_dir / f"{index}.fa.gz").is_file():
                index_path = self.data_dir / f"{index}.fa.gz"
                if get_mmi_path(index_path).is_file():
                    logging.info(f"Found cached standard index {index} (MMI available)")
                else:
                    logging.info(
                        f"Found cached standard index {index} (building MMI cache; do not interrupt)"
                    )
            elif not airplane and util.fetch_manifest(util.INDEX_REPOSITORY_URL).get(
                index
            ):
                file_name = f"{index}.fa.gz"
                file_url = f"{util.INDEX_REPOSITORY_URL}/{file_name}"
                logging.info(f"Fetching standard index {index} ({file_url})")
                manifest = util.fetch_manifest(util.INDEX_REPOSITORY_URL)
                with tempfile.NamedTemporaryFile() as temporary_file:
                    tmp_path = Path(temporary_file.name)
                    util.download(f"{file_url}", tmp_path)
                    expected_sha256 = manifest[index]["assets"][file_name]["sha256"]
                    logging.info(f"Verifying checksum {expected_sha256}…")
                    observed_sha256 = util.sha256(tmp_path)
                    if observed_sha256 != expected_sha256:
                        raise ValueError(f"Checksum mismatch for {file_name}")
                    shutil.copy(tmp_path, self.data_dir / f"{index}.fa.gz")
                index_path = self.data_dir / f"{index}.fa.gz"
                logging.info(f"Downloaded standard index {index_path}")
                mmi_path = get_mmi_path(index_path)
                logging.info(f"Building MMI cache; do not interrupt {mmi_path}…")
                util.run(f'minimap2 -x map-ont -d "{mmi_path}" "{index_path}"')
            else:
                message = f"{index} is neither a valid custom FASTA path, nor a valid standard index name. Mode: long read (Minimap2)"
                if airplane:
                    message += ". Disable airplane mode to enable discovery of standard indexes"
                raise FileNotFoundError(message)

        return index_path

    def gen_clean_cmd(
        self,
        fastq: Path | str,
        index_path: Path,
        invert: bool,
        rename: bool,
        reorder: bool,
        casava: bool,
        stdin: bool,
        stdout: bool,
        output: Path | str,
        aligner_args: str,
        aligner_threads: int,
        compression_threads: int,
        force: bool,
    ) -> str:
        fastq, output = Path(fastq), Path(output)
        output.mkdir(exist_ok=True, parents=True)
        fastq_stem = util.fastq_path_to_stem(fastq)
        fastq_out_path = output / f"{fastq_stem}.clean.fastq.gz"
        count_before_path = output / f"{fastq_stem}.reads_in.txt"
        count_after_path = output / f"{fastq_stem}.reads_out.txt"

        if not stdout and not force and fastq_out_path.exists():
            raise FileExistsError(
                "Output file already exists. Use --force to overwrite"
            )

        filter_cmd = (
            " | samtools view -hF 4 -" if invert else " | samtools view -hf 4 -"
        )
        reorder_cmd = " | samtools sort -n -O sam -@ 6 -m 1G" if reorder else ""
        rename_cmd = (
            # Preserve header (^@) lines but do not start count until first non ^@ line
            ' | awk \'BEGIN {{ FS=OFS="\\t"; line_count=0 }} /^@/ {{ print $0; next }}'
            " {{ $1=int(line_count+1); print $0; line_count++ }}'"
            if rename
            else ""
        )

        mmi_path = get_mmi_path(index_path)

        cmd_template = {
            "{BIN_PATH}": str(self.bin_path),
            "{INDEX_PATH}": str(index_path),
            "{MMI_PATH}": str(mmi_path),
            "{FASTQ1}": str(fastq),
            "{ALIGNER_ARGS}": str(aligner_args),
            "{ALIGNER_THREADS}": str(aligner_threads),
        }

        if self.name == "Minimap2" and not mmi_path.is_file():
            alignment_cmd = self.single_unindexed_cmd
        else:
            alignment_cmd = self.single_cmd

        for k in cmd_template.keys():
            alignment_cmd = alignment_cmd.replace(k, cmd_template[k])

        if casava:
            header_fmt = "-i --index-format 'i*'"
        else:
            header_fmt = ""

        if stdout:
            fastq_cmd = f"samtools fastq --threads 0 {header_fmt} -c 6 -0 -"
        else:
            fastq_cmd = f"samtools fastq --threads {compression_threads} -c 6 {header_fmt} -0 '{fastq_out_path}'"

        cmd = (
            # Align, stream reads to stdout in SAM format
            f"{alignment_cmd}"
            # Count reads in stream before filtering (2048 + 256 = 2304)
            f" | tee >(samtools view -F 2304 -c - > '{count_before_path}')"
            # Discard mapped reads (or inverse)
            f"{filter_cmd}"
            # Count reads in stream after filtering (2048 + 256 = 2304)
            f" | tee >(samtools view -F 2304 -c - > '{count_after_path}')"
            # Optionally sort reads by name
            f"{reorder_cmd}"
            # Optionally replace read headers with integers
            f"{rename_cmd}"
            # Stream remaining records into fastq
            f" | {fastq_cmd}"
        )

        return cmd

    def gen_paired_clean_cmd(
        self,
        fastq1: Path | str,
        fastq2: Path | str,
        index_path: Path,
        invert: bool,
        rename: bool,
        reorder: bool,
        casava: bool,
        stdin: bool,
        stdout: bool,
        output: Path | str,
        aligner_args: str,
        aligner_threads: int,
        compression_threads: int,
        force: bool,
    ) -> str:
        fastq1, fastq2, output = Path(fastq1), Path(fastq2), Path(output)
        output.mkdir(exist_ok=True, parents=True)
        fastq1_stem = util.fastq_path_to_stem(fastq1)
        fastq2_stem = util.fastq_path_to_stem(fastq2)
        fastq1_out_path = output / f"{fastq1_stem}.clean_1.fastq.gz"
        fastq2_out_path = output / f"{fastq2_stem}.clean_2.fastq.gz"
        count_before_path = output / f"{fastq1_stem}.reads_in.txt"
        count_after_path = output / f"{fastq1_stem}.reads_out.txt"

        if (
            not stdout
            and not force
            and (fastq1_out_path.exists() or fastq2_out_path.exists())
        ):
            raise FileExistsError(
                "Output files already exist. Use --force to overwrite"
            )

        filter_cmd = (
            " | samtools view -h -e 'flag.unmap == 0 || flag.munmap == 0' -"
            if invert
            else " | samtools view -hf 12 -"
        )
        reorder_cmd = ""
        if self.name == "Bowtie2" and reorder:
            if util.get_platform() == "darwin":
                reorder_cmd = " | samtools sort -n -O sam -@ 6 -m 1G"
            else:  # Under Linux, Bowtie2's --reorder option is efficient
                reorder_cmd = ""
                aligner_args += " --reorder"
        rename_cmd = (
            # Preserve header (^@) lines but do not start counting until first non ^@ line
            ' | awk \'BEGIN {{ FS=OFS="\\t"; start=0; line_count=1 }} /^@/ {{ print $0; next }}'
            " !start && !/^@/ {{ start=1 }} start {{ $1=int((line_count+1)/2);"
            " print $0; line_count++ }}'"
            if rename
            else ""
        )

        mmi_path = get_mmi_path(index_path)

        cmd_template = {
            "{BIN_PATH}": str(self.bin_path),
            "{INDEX_PATH}": str(index_path),
            "{MMI_PATH}": str(mmi_path),
            "{FASTQ1}": str(fastq1),
            "{FASTQ2}": str(fastq2),
            "{ALIGNER_ARGS}": str(aligner_args),
            "{ALIGNER_THREADS}": str(aligner_threads),
        }

        if self.name == "Minimap2":
            logging.warning(
                "Minimap2 is not recommended for decontaminating short (paired) reads"
            )

        if self.name == "Minimap2":
            if not mmi_path.is_file():  # No MMI, make one
                if stdin:  # Interleaved stdin
                    alignment_cmd = self.interleaved_unindexed_cmd
                else:  # Separate fastq1 and fastq2 file input
                    alignment_cmd = self.paired_unindexed_cmd
            else:  # MMI exists
                if stdin:  # Interleaved stdin
                    alignment_cmd = self.interleaved_cmd
                else:  # Separate fastq1 and fastq2 file input
                    alignment_cmd = self.paired_cmd
        else:  # Bowtie2
            if stdin:  # Interleaved stdin
                alignment_cmd = self.interleaved_cmd
            else:  # Separate fastq1 and fastq2 file input
                alignment_cmd = self.paired_cmd

        for k in cmd_template.keys():
            alignment_cmd = alignment_cmd.replace(k, cmd_template[k])

        if casava:
            header_fmt = "-n -i --index-format 'i*'"
        else:
            header_fmt = "-N"

        if stdout:
            fastq_cmd = f"samtools fastq --threads 0 -c 6 {header_fmt} -0 -"
        else:
            fastq_cmd = (
                f"samtools fastq --threads {compression_threads} -c 6 {header_fmt}"
                f" -1 '{fastq1_out_path}' -2 '{fastq2_out_path}'"
                f" -0 /dev/null -s /dev/null"
            )
        cmd = (
            f"{alignment_cmd}"
            f" | tee >(samtools view -F 2304 -c - > '{count_before_path}')"
            f"{filter_cmd}"
            f" | tee >(samtools view -F 2304 -c - > '{count_after_path}')"
            f"{reorder_cmd}"
            f"{rename_cmd}"
            f" | {fastq_cmd}"
        )

        return cmd


ALIGNER = Enum(
    "Aligner",
    {
        "bowtie2": Aligner(
            name="Bowtie2",
            bin_path=Path("bowtie2"),
            data_dir=util.CACHE_DIR,
            single_cmd=(
                "'{BIN_PATH}' -x '{INDEX_PATH}' -U '{FASTQ1}'"
                " -k 1 --mm -p {ALIGNER_THREADS} {ALIGNER_ARGS}"
            ),
            paired_cmd=(
                "{BIN_PATH} -x '{INDEX_PATH}' -1 '{FASTQ1}' -2 '{FASTQ2}'"
                " -k 1 --mm -p {ALIGNER_THREADS} {ALIGNER_ARGS}"
            ),
            interleaved_cmd=(
                "{BIN_PATH} -x '{INDEX_PATH}' --interleaved -"
                " -k 1 --mm -p {ALIGNER_THREADS} {ALIGNER_ARGS}"
            ),
            single_unindexed_cmd="",
            paired_unindexed_cmd="",
            interleaved_unindexed_cmd="",
        ),
        "minimap2": Aligner(
            name="Minimap2",
            bin_path=Path("minimap2"),
            data_dir=util.CACHE_DIR,
            single_cmd=(
                "'{BIN_PATH}' -ax map-ont --secondary no -t {ALIGNER_THREADS}"
                " {ALIGNER_ARGS} '{MMI_PATH}' '{FASTQ1}'"
            ),
            single_unindexed_cmd=(
                "'{BIN_PATH}' -ax map-ont --secondary no -t {ALIGNER_THREADS}"
                " {ALIGNER_ARGS} -d '{MMI_PATH}' '{INDEX_PATH}' '{FASTQ1}'"
            ),
            paired_cmd=(
                "'{BIN_PATH}' -ax sr --secondary no -t {ALIGNER_THREADS} {ALIGNER_ARGS}"
                " '{MMI_PATH}' '{FASTQ1}' '{FASTQ2}'"
            ),
            paired_unindexed_cmd=(
                "'{BIN_PATH}' -ax sr --secondary no -t {ALIGNER_THREADS} {ALIGNER_ARGS}"
                " -d '{MMI_PATH}' '{INDEX_PATH}' '{FASTQ1}' '{FASTQ2}'"
            ),
            interleaved_cmd=(
                "'{BIN_PATH}' -ax sr --secondary no -t {ALIGNER_THREADS}"
                " {ALIGNER_ARGS} '{MMI_PATH}' '{FASTQ1}'"
            ),
            interleaved_unindexed_cmd=(
                "'{BIN_PATH}' -ax sr --secondary no -t {ALIGNER_THREADS}"
                " {ALIGNER_ARGS} -d '{MMI_PATH}' '{INDEX_PATH}' '{FASTQ1}'"
            ),
        ),
    },
)
