import logging
import subprocess

from dataclasses import dataclass
from pathlib import Path

from hostile import util


@dataclass
class Aligner:
    name: str
    short_name: str
    bin_path: Path
    cdn_base_url: str
    working_dir: Path
    cmd: str
    idx_archive_fn: str = ""
    ref_archive_fn: str = ""
    idx_name: str = ""
    idx_paths: tuple[Path] = tuple()

    def __post_init__(self):
        self.ref_archive_url = f"{self.cdn_base_url}/{self.ref_archive_fn}"
        self.idx_archive_url = f"{self.cdn_base_url}/{self.idx_archive_fn}"
        self.ref_archive_path = self.working_dir / self.ref_archive_fn
        self.idx_archive_path = self.working_dir / self.idx_archive_fn
        self.idx_path = self.working_dir / self.idx_name

    def check(self):
        logging.info(f"Using {self.name}")
        if self.name == "Bowtie2":
            if not all(path.exists() for path in self.idx_paths):
                self.working_dir.mkdir(exist_ok=True, parents=True)
                logging.info(f"Fetching human index")
                util.download(self.idx_archive_url, self.idx_archive_path)
                util.untar_file(self.idx_archive_path, self.working_dir)
                self.idx_archive_path.unlink()
                logging.info(f"Saved human index ({self.idx_path})")
            else:
                logging.info(f"Using cached human index ({self.idx_path})")
        elif self.name == "Minimap2":
            if not self.ref_archive_path.exists():
                util.download(self.idx_archive_url, self.idx_archive_path)
                logging.info(f"Saved human reference ({self.ref_archive_path})")
            else:
                logging.info(f"Using cached human reference ({self.ref_archive_path})")
        try:
            util.run(f"{self.bin_path} --help", cwd=self.working_dir)
        except subprocess.CalledProcessError:
            raise RuntimeError(f"Failed to execute {self.bin_path}")

    def gen_paired_dehost_cmd(
        self, fastq1: Path, fastq2: Path, out_dir: Path, threads: int = 2
    ) -> str:
        fastq1, fastq2, out_dir = Path(fastq1), Path(fastq2), Path(out_dir)
        out_dir.mkdir(exist_ok=True, parents=True)
        fastq1_stem = util.fastq_path_to_stem(fastq1)
        fastq2_stem = util.fastq_path_to_stem(fastq2)
        fastq1_out_path = out_dir / f"{fastq1_stem}.dehosted_1.fastq.gz"
        fastq2_out_path = out_dir / f"{fastq2_stem}.dehosted_2.fastq.gz"
        count_before_path = out_dir / f"{fastq1_stem}.reads_in.txt"
        count_after_path = out_dir / f"{fastq1_stem}.reads_out.txt"
        cmd_template = {  # Templating for Aligner.cmd
            "{BIN_PATH}": str(self.bin_path),
            "{REF_ARCHIVE_PATH}": str(self.ref_archive_path),
            "{INDEX_PATH}": str(self.idx_path),
            "{FASTQ1}": str(fastq1),
            "{FASTQ2}": str(fastq2),
            "{THREADS}": str(threads),
        }
        for k in cmd_template.keys():
            self.cmd = self.cmd.replace(k, cmd_template[k])
        cmd = (
            # Align, stream reads to stdout in SAM format
            f"{self.cmd}"
            # Count reads in stream before filtering
            f" | tee >(samtools view -F 256 -c - > '{count_before_path}')"
            # Discard mapped reads and reads with mapped mates
            f" | samtools view --threads {int(threads/2)} -f 12 -"
            # Count reads in stream after filtering
            f" | tee >(samtools view -F 256 -c - > '{count_after_path}')"
            # Replace paired read headers with integers
            f' | awk \'BEGIN{{FS=OFS="\\t"}} {{$1=int((NR+1)/2)" "; print $0}}\''
            # Stream remaining records into fastq files
            f" | samtools fastq --threads {int(threads/2)} -c 6 -N -1 '{fastq1_out_path}' -2 '{fastq2_out_path}'"
        )
        return cmd
