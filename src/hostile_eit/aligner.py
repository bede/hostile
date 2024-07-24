import logging
import shutil
import subprocess
import tempfile

from dataclasses import dataclass
from enum import Enum
from pathlib import Path

from hostile_eit import util


@dataclass
class Aligner:
    name: str
    short_name: str
    bin_path: Path
    data_dir: Path
    cmd: str
    paired_cmd: str

    def __post_init__(self):
        Path(self.data_dir).mkdir(exist_ok=True, parents=True)

    def check_index(self, index: str, offline: bool = False) -> Path:
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
            elif not offline and util.fetch_manifest(util.BUCKET_URL).get(index):
                file_name = f"{index}.tar"
                file_url = f"{util.BUCKET_URL}/{file_name}"
                logging.info(f"Fetching standard index {index} ({file_url})")
                manifest = util.fetch_manifest(util.BUCKET_URL)
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
                message = f"{index} is neither a valid custom index path nor a valid standard index name"
                if offline:
                    message += (
                        ". Disable offline mode to enable discovery of standard indexes"
                    )
                raise FileNotFoundError(message)
        elif self.name == "Minimap2":
            if Path(f"{index}").is_file():
                index_path = Path(index)
                logging.info(f"Found custom index {index}")
            elif (self.data_dir / f"{index}.fa.gz").is_file():
                index_path = self.data_dir / f"{index}.fa.gz"
                logging.info(f"Found cached standard index {index}")
            elif not offline and util.fetch_manifest(util.BUCKET_URL).get(index):
                file_name = f"{index}.fa.gz"
                file_url = f"{util.BUCKET_URL}/{file_name}"
                logging.info(f"Fetching standard index {index} ({file_url})")
                manifest = util.fetch_manifest(util.BUCKET_URL)
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
            else:
                message = f"{index} is neither a valid custom index path nor a valid standard index name"
                if offline:
                    message += (
                        ". Disable offline mode to enable discovery of standard indexes"
                    )
                raise FileNotFoundError(message)
        return index_path

    def gen_clean_cmd(
        self,
        fastq: Path,
        out_dir: Path,
        index_path: Path,
        invert: bool,
        rename: bool,
        reorder: bool,
        aligner_args: str,
        threads: int,
        force: bool,
        offline: bool,
    ) -> str:
        fastq, out_dir = Path(fastq), Path(out_dir)
        out_dir.mkdir(exist_ok=True, parents=True)
        fastq_stem = util.fastq_path_to_stem(fastq)
        fastq_out_path = out_dir / f"{fastq_stem}.clean.fastq.gz"
        count_before_path = out_dir / f"{fastq_stem}.reads_in.txt"
        count_after_path = out_dir / f"{fastq_stem}.reads_out.txt"
        if not force and fastq_out_path.exists():
            raise FileExistsError(
                "Output file already exists. Use --force to overwrite"
            )
        filter_cmd = " | samtools view -hF 4 -" if invert else " | samtools view -f 4 -"
        reorder_cmd = " | samtools sort -n -O sam -@ 6 -m 1G" if reorder else ""
        rename_cmd = (
            # Preserve header (^@) lines but do not start counting until first non ^@ line
            ' | awk \'BEGIN {{ FS=OFS="\\t"; line_count=0 }} /^@/ {{ print $0; next }}'
            ' {{ $1=int(line_count+1)" "; print $0; line_count++ }}\''
            if rename
            else ""
        )
        cmd_template = {  # Templating for Aligner.cmd
            "{BIN_PATH}": str(self.bin_path),
            "{INDEX_PATH}": str(index_path),
            "{FASTQ}": str(fastq),
            "{ALIGNER_ARGS}": str(aligner_args),
            "{THREADS}": str(threads),
        }
        alignment_cmd = self.cmd
        for k in cmd_template.keys():
            alignment_cmd = alignment_cmd.replace(k, cmd_template[k])
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
            # Stream remaining records into fastq files
            f" | samtools fastq --threads 4 -c 6 -0 '{fastq_out_path}'"
        )
        return cmd

    def gen_paired_clean_cmd(
        self,
        fastq1: Path,
        fastq2: Path,
        out_dir: Path,
        index_path: Path,
        invert: bool,
        rename: bool,
        reorder: bool,
        aligner_args: str,
        threads: int,
        force: bool,
        offline: bool,
    ) -> str:
        fastq1, fastq2, out_dir = Path(fastq1), Path(fastq2), Path(out_dir)
        out_dir.mkdir(exist_ok=True, parents=True)
        fastq1_stem = util.fastq_path_to_stem(fastq1)
        fastq2_stem = util.fastq_path_to_stem(fastq2)
        fastq1_out_path = out_dir / f"{fastq1_stem}.clean_1.fastq.gz"
        fastq2_out_path = out_dir / f"{fastq2_stem}.clean_2.fastq.gz"
        count_before_path = out_dir / f"{fastq1_stem}.reads_in.txt"
        count_after_path = out_dir / f"{fastq1_stem}.reads_out.txt"
        if not force and (fastq1_out_path.exists() or fastq2_out_path.exists()):
            raise FileExistsError(
                "Output files already exist. Use --force to overwrite"
            )
        filter_cmd = (
            " | samtools view -hF 12 -" if invert else " | samtools view -f 12 -"
        )
        reorder_cmd = ""
        if self.name == "Bowtie2" and reorder:
            if (
                util.get_platform() == "darwin"
            ):  # Under MacOS, Bowtie2's native --reorder is very slow
                reorder_cmd = " | samtools sort -n -O sam -@ 6 -m 1G" if reorder else ""
            else:  # Under Linux, Bowtie2's --reorder option works very well
                reorder_cmd = ""
                aligner_args += " --reorder"
        rename_cmd = (
            # Preserve header (^@) lines but do not start counting until first non ^@ line
            ' | awk \'BEGIN {{ FS=OFS="\\t"; start=0; line_count=1 }} /^@/ {{ print $0; next }}'
            ' !start && !/^@/ {{ start=1 }} start {{ $1=int((line_count+1)/2)" ";'
            " print $0; line_count++ }}'"
            if rename
            else ""
        )
        cmd_template = {  # Templating for Aligner.cmd
            "{BIN_PATH}": str(self.bin_path),
            "{INDEX_PATH}": str(index_path),
            "{FASTQ1}": str(fastq1),
            "{FASTQ2}": str(fastq2),
            "{ALIGNER_ARGS}": str(aligner_args),
            "{THREADS}": str(threads),
        }
        alignment_cmd = self.paired_cmd
        for k in cmd_template.keys():
            alignment_cmd = alignment_cmd.replace(k, cmd_template[k])
        cmd = (
            # Align, stream reads to stdout in SAM format
            f"{alignment_cmd}"
            # Count reads in stream before filtering (2048 + 256 = 2304)
            f" | tee >(samtools view -F 2304 -c - > '{count_before_path}')"
            # Discard mapped reads and reads with mapped mates (or inverse)
            f"{filter_cmd}"
            # Count reads in stream after filtering (2048 + 256 = 2304)
            f" | tee >(samtools view -F 2304 -c - > '{count_after_path}')"
            # Optionally sort reads by name
            f"{reorder_cmd}"
            # Optionally replace paired read headers with integers
            f"{rename_cmd}"
            # Stream remaining records into fastq files
            f" | samtools fastq --threads 4 -c 6 -N -1 '{fastq1_out_path}' -2 '{fastq2_out_path}'"
        )
        return cmd


ALIGNER = Enum(
    "Aligner",
    {
        "bowtie2": Aligner(
            name="Bowtie2",
            short_name="bt2",
            bin_path=Path("bowtie2"),
            data_dir=util.CACHE_DIR,
            cmd=(
                "{BIN_PATH} -x '{INDEX_PATH}' -U '{FASTQ}'"
                " -k 1 --mm -p {THREADS} {ALIGNER_ARGS}"
            ),
            paired_cmd=(
                "{BIN_PATH} -x '{INDEX_PATH}' -1 '{FASTQ1}' -2 '{FASTQ2}'"
                " -k 1 --mm -p {THREADS} {ALIGNER_ARGS}"
            ),
        ),
        "minimap2": Aligner(
            name="Minimap2",
            short_name="mm2",
            bin_path=Path("minimap2"),
            data_dir=util.CACHE_DIR,
            cmd="{BIN_PATH} -ax map-ont -m 40 --secondary no -t {THREADS} {ALIGNER_ARGS} '{INDEX_PATH}' '{FASTQ}'",
            paired_cmd="{BIN_PATH} -ax sr -m 40 --secondary no -t {THREADS} {ALIGNER_ARGS} '{INDEX_PATH}' '{FASTQ1}' '{FASTQ2}'",
        ),
    },
)
