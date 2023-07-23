import logging
import shutil
import subprocess
import tempfile

from dataclasses import dataclass
from pathlib import Path

from hostile import util


@dataclass
class Aligner:
    name: str
    short_name: str
    bin_path: Path
    cdn_base_url: str
    data_dir: Path
    cmd: str
    paired_cmd: str
    idx_archive_fn: str = ""
    ref_archive_fn: str = ""
    idx_name: str = ""
    idx_paths: tuple[Path, ...] = tuple()

    def __post_init__(self):
        self.ref_archive_url = f"{self.cdn_base_url}/{self.ref_archive_fn}"
        self.idx_archive_url = f"{self.cdn_base_url}/{self.idx_archive_fn}"
        self.ref_archive_path = self.data_dir / self.ref_archive_fn
        self.idx_archive_path = self.data_dir / self.idx_archive_fn
        self.idx_path = self.data_dir / self.idx_name
        Path(self.data_dir).mkdir(exist_ok=True, parents=True)

    def check(self, using_custom_index: bool):
        """Test aligner and check/download a ref/index if necessary"""
        if not using_custom_index:  # Check for and if necessary fetch a genome/index
            if self.name == "Bowtie2":
                if not all(path.exists() for path in self.idx_paths):
                    self.data_dir.mkdir(exist_ok=True, parents=True)
                    logging.info(f"Fetching human index ({self.idx_archive_url})")
                    with tempfile.NamedTemporaryFile() as temporary_file:
                        tmp_path = Path(temporary_file.name)
                        util.download(self.idx_archive_url, tmp_path)
                        logging.info("Extracting index…")
                        util.untar_file(tmp_path, self.data_dir)
                    logging.info(f"Saved human index ({self.idx_path})")
                else:
                    logging.info(f"Found cached index ({self.idx_path})")
            elif self.name == "Minimap2":
                if not self.ref_archive_path.exists():
                    with tempfile.NamedTemporaryFile(delete=False) as temporary_file:
                        tmp_path = Path(temporary_file.name)
                        util.download(self.ref_archive_url, tmp_path)
                        shutil.move(tmp_path, self.ref_archive_path)
                    logging.info(f"Saved human reference ({self.ref_archive_path})")
                else:
                    logging.info(f"Found cached reference ({self.ref_archive_path})")
        try:
            util.run(f"{self.bin_path} --version", cwd=self.data_dir)
        except subprocess.CalledProcessError:
            logging.warning(f"Failed to execute {self.bin_path}")
            raise RuntimeError(f"Failed to execute {self.bin_path}")

    def gen_clean_cmd(
        self,
        fastq: Path,
        out_dir: Path,
        index: Path | None,
        rename: bool,
        threads: int,
        force: bool,
    ) -> str:
        fastq, out_dir = Path(fastq), Path(out_dir)
        out_dir.mkdir(exist_ok=True, parents=True)
        fastq_stem = util.fastq_path_to_stem(fastq)
        fastq_out_path = out_dir / f"{fastq_stem}.clean.fastq.gz"
        count_before_path = out_dir / f"{fastq_stem}.reads_in.txt"
        count_after_path = out_dir / f"{fastq_stem}.reads_out.txt"
        if not force and fastq_out_path.exists():
            raise FileExistsError(
                f"Output file already exists. Use --force to overwrite"
            )
        if index:
            self.idx_path = Path(index)
            self.ref_archive_path = Path(index)
            logging.info(f"Using custom index {index}")
        cmd_template = {  # Templating for Aligner.cmd
            "{BIN_PATH}": str(self.bin_path),
            "{REF_ARCHIVE_PATH}": str(self.ref_archive_path),
            "{INDEX_PATH}": str(self.idx_path),
            "{FASTQ}": str(fastq),
            "{THREADS}": str(threads),
        }
        alignment_cmd = self.cmd
        for k in cmd_template.keys():
            alignment_cmd = alignment_cmd.replace(k, cmd_template[k])
        rename_cmd = (
            ' | awk \'BEGIN{{FS=OFS="\\t"}} {{$1=int((NR+1)/2)" "; print $0}}\''
            if rename
            else ""
        )
        cmd = (
            # Align, stream reads to stdout in SAM format
            f"{alignment_cmd}"
            # Count reads in stream before filtering (2048 + 256 = 2304)
            f" | tee >(samtools view -F 2304 -c - > '{count_before_path}')"
            # Discard mapped reads
            f" | samtools view -f 4 -"
            # Count reads in stream after filtering
            f" | tee >(samtools view -F 256 -c - > '{count_after_path}')"
            # Optionally replace paired read headers with integers
            f"{rename_cmd}"
            # Stream remaining records into fastq files
            f" | samtools fastq --threads 2 -c 6 -0 '{fastq_out_path}'"
        )
        return cmd

    def gen_paired_clean_cmd(
        self,
        fastq1: Path,
        fastq2: Path,
        out_dir: Path,
        index: Path | None,
        rename: bool,
        threads: int,
        force: bool,
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
                f"Output files already exist. Use --force to overwrite"
            )
        if index:
            self.idx_path = Path(index)
            self.ref_archive_path = Path(index)
            logging.info(f"Using custom index ({index})")
        cmd_template = {  # Templating for Aligner.cmd
            "{BIN_PATH}": str(self.bin_path),
            "{REF_ARCHIVE_PATH}": str(self.ref_archive_path),
            "{INDEX_PATH}": str(self.idx_path),
            "{FASTQ1}": str(fastq1),
            "{FASTQ2}": str(fastq2),
            "{THREADS}": str(threads),
        }
        alignment_cmd = self.paired_cmd
        for k in cmd_template.keys():
            alignment_cmd = alignment_cmd.replace(k, cmd_template[k])
        rename_cmd = (
            f' | awk \'BEGIN{{FS=OFS="\\t"}} {{$1=int((NR+1)/2)" "; print $0}}\''
            if rename
            else ""
        )
        cmd = (
            # Align, stream reads to stdout in SAM format
            f"{alignment_cmd}"
            # Count reads in stream before filtering (2048 + 256 = 2304)
            f" | tee >(samtools view -F 2304 -c - > '{count_before_path}')"
            # Discard mapped reads and reads with mapped mates
            f" | samtools view -f 12 -"
            # Count reads in stream after filtering
            f" | tee >(samtools view -F 256 -c - > '{count_after_path}')"
            # Optionally replace paired read headers with integers
            f"{rename_cmd}"
            # Stream remaining records into fastq files
            f" | samtools fastq --threads 2 -c 6 -N -1 '{fastq1_out_path}' -2 '{fastq2_out_path}'"
        )
        return cmd
