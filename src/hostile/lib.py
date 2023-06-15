import hashlib
import logging
import multiprocessing
import subprocess
import tarfile
from pathlib import Path

from dataclasses import dataclass
from enum import Enum

import httpx

from platformdirs import user_data_dir
from tqdm import tqdm


class AlignerError(Exception):
    def __init__(self, message: str):
        super().__init__(message)


ALIGNERS = Enum("Aligner", {"bowtie2": "bowtie2", "minimap2": "minimap2"})


logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)

CWD = Path.cwd()
XDG_DATA_DIR = Path(user_data_dir("hostile", "Bede Constantinides"))
THREADS = multiprocessing.cpu_count()


def run(cmd, cwd=CWD):  # Helper for CLI testing
    return subprocess.run(
        cmd, cwd=cwd, shell=True, check=True, text=True, capture_output=True
    )


@dataclass
class Backend:
    name: str
    short_name: str
    bin_path: Path
    cdn_base_url: str
    cmd: str
    idx_archive_fn: str = ""
    ref_archive_fn: str = ""
    idx_name: str = ""
    idx_paths: tuple[Path] = tuple()

    def __post_init__(self):
        self.ref_archive_url = f"{self.cdn_base_url}/{self.ref_archive_fn}"
        self.idx_archive_url = f"{self.cdn_base_url}/{self.idx_archive_fn}"
        self.ref_archive_path = XDG_DATA_DIR / self.ref_archive_fn
        self.idx_archive_path = XDG_DATA_DIR / self.idx_archive_fn
        # self.ref_path = XDG_DATA_DIR / self.idx_name
        self.idx_path = XDG_DATA_DIR / self.idx_name

    def dehost_paired_fastq(
        self, fastq1: Path, fastq2: Path, out_dir: Path = CWD, threads: int = 2
    ) -> dict[str, str]:
        fastq1, fastq2, out_dir = Path(fastq1), Path(fastq2), Path(out_dir)
        out_dir.mkdir(exist_ok=True)
        fastq1_stem = fastq1.name.removesuffix(fastq1.suffixes[-1]).removesuffix(
            fastq1.suffixes[-2]
        )
        fastq2_stem = fastq2.name.removesuffix(fastq2.suffixes[-1]).removesuffix(
            fastq2.suffixes[-2]
        )
        fastq1_out_path = out_dir / f"{fastq1_stem}.dehosted_1.fastq.gz"
        fastq2_out_path = out_dir / f"{fastq2_stem}.dehosted_2.fastq.gz"
        cmd_template = {  # Templating for Backend.cmd
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
            f"{self.cmd} | samtools view --threads {threads/2} -f 12 -"
            f' | awk \'BEGIN{{FS=OFS="\\t"}} {{$1=int((NR+1)/2)" "; print $0}}\''
            f" | samtools fastq --threads {threads/2} -c 6 -N -1 '{fastq1_out_path}' -2 '{fastq2_out_path}'"
        )
        logging.info(f"Decontaminating {fastq1} and {fastq2}")
        logging.debug(f"{cmd}")
        run(cmd, cwd=CWD)
        checksums = {
            p.name: sha256sum(p)
            for p in (fastq1, fastq1_out_path, fastq2, fastq2_out_path)
        }
        logging.info(f"Decontaminated {fastq1} and {fastq2}")
        return checksums

    def check(self, threads: int) -> None:
        logging.info(f"Using {self.name} ({threads} threads)")
        if self.name == "Bowtie2":
            if not all(path.exists() for path in self.idx_paths):
                XDG_DATA_DIR.mkdir(exist_ok=True)
                logging.info(f"Fetching human index")
                download(self.idx_archive_url, self.idx_archive_path)
                untar_file(self.idx_archive_path, XDG_DATA_DIR)
                self.idx_archive_path.unlink()
                logging.info(f"Saved human index ({self.idx_path})")
            else:
                logging.info(f"Using cached human index ({self.idx_path})")
        elif self.name == "Minimap2":
            if not self.ref_archive_path.exists():
                download(self.idx_archive_url, self.idx_archive_path)
                logging.info(f"Saved human reference ({self.ref_archive_path})")
            else:
                logging.info(f"Using cached human reference ({self.ref_archive_path})")
        try:
            run(f"{self.bin_path} --help")
        except subprocess.CalledProcessError:
            raise AlignerError(f"Failed to execute {self.bin_path}")


backends = {
    "bowtie2": Backend(
        name="Bowtie2",
        short_name="bt2",
        bin_path=Path("/Users/bede/Downloads/bowtie2-2.5.1-macos-arm64/bowtie2"),
        cdn_base_url=f"http://178.79.139.243/hostile",
        cmd=(
            "{BIN_PATH} -x '{INDEX_PATH}' -1 '{FASTQ1}' -2 '{FASTQ2}'"
            " -k 1 --mm -p {THREADS}"
        ),
        idx_archive_fn="human-bowtie2.tar",
        idx_name="human-bowtie2",
        idx_paths=(
            XDG_DATA_DIR / "human-bowtie2.1.bt2",
            XDG_DATA_DIR / "human-bowtie2.2.bt2",
            XDG_DATA_DIR / "human-bowtie2.3.bt2",
            XDG_DATA_DIR / "human-bowtie2.4.bt2",
            XDG_DATA_DIR / "human-bowtie2.rev.1.bt2",
            XDG_DATA_DIR / "human-bowtie2.rev.2.bt2",
        ),
    ),
    "minimap2": Backend(
        name="Minimap2",
        short_name="mm2",
        bin_path=Path("minimap2"),
        cdn_base_url=f"http://178.79.139.243/hostile",
        cmd="{BIN_PATH} -ax sr -m 40 -t {THREADS} '{REF_ARCHIVE_PATH}' '{FASTQ1}' '{FASTQ2}'",
        ref_archive_fn="human.fa.gz",
        idx_name="human.fa.gz",
    ),
}


def untar_file(input_path, output_path):
    with tarfile.open(input_path) as fh:
        fh.extractall(path=output_path)


def sha256sum(filename) -> str:
    with open(filename, "rb", buffering=0) as fh:
        return hashlib.file_digest(fh, "sha256").hexdigest()


def download(url: str, path: Path) -> None:
    with open(path, "wb") as fh:
        with httpx.stream("GET", url) as response:
            total = int(response.headers["Content-Length"])
            with tqdm(
                total=total, unit_scale=True, unit_divisor=1024, unit="B"
            ) as progress:
                num_bytes_downloaded = response.num_bytes_downloaded
                for chunk in response.iter_bytes():
                    fh.write(chunk)
                    progress.update(
                        response.num_bytes_downloaded - num_bytes_downloaded
                    )
                    num_bytes_downloaded = response.num_bytes_downloaded


def dehost_paired_fastqs(
    fastqs: list[tuple[Path, Path]],
    out_dir: Path = CWD,
    threads: int = THREADS,
    aligner: ALIGNERS = ALIGNERS.bowtie2,
) -> dict[str, str]:
    backend = backends[aligner.name]
    backend.check(threads)
    checksums = {}
    for fastq_pair in fastqs:
        checksums.update(
            backend.dehost_paired_fastq(
                fastq_pair[0], fastq_pair[1], out_dir=out_dir, threads=threads
            )
        )
    return checksums
