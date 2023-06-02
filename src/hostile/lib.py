import hashlib
import logging
import subprocess
import tarfile
from pathlib import Path

import httpx
from platformdirs import user_data_dir
from tqdm import tqdm


logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)

CWD = Path.cwd()
XDG_DATA_DIR = Path(user_data_dir("hostile", "Bede Constantinides"))
BT2_INDEX_NAME = "human-bowtie2"
BT2_ARCHIVE_FILENAME = f"{BT2_INDEX_NAME}.tar"
BT2_ARCHIVE_URL = f"http://178.79.139.243/hostile/{BT2_ARCHIVE_FILENAME}"
BT2_ARCHIVE_PATH = Path(XDG_DATA_DIR / BT2_ARCHIVE_FILENAME)
BT2_INDEX_PATH = XDG_DATA_DIR / BT2_INDEX_NAME
BT2_INDEX_PATHS = [
    XDG_DATA_DIR / "human-bowtie2.1.bt2",
    XDG_DATA_DIR / "human-bowtie2.2.bt2",
    XDG_DATA_DIR / "human-bowtie2.3.bt2",
    XDG_DATA_DIR / "human-bowtie2.4.bt2",
    XDG_DATA_DIR / "human-bowtie2.rev.1.bt2",
    XDG_DATA_DIR / "human-bowtie2.rev.2.bt2"
]
logging.debug(f"{CWD=} {XDG_DATA_DIR=}")


def run(cmd, cwd=CWD):  # Helper for CLI testing
    return subprocess.run(
        cmd, cwd=cwd, shell=True, check=True, text=True, capture_output=True
)


def untar_file(input_path, output_path):
    with tarfile.open(input_path) as fh:
        fh.extractall(path=output_path)


def sha256sum(filename) -> str:
    with open(filename, 'rb', buffering=0) as fh:
        return hashlib.file_digest(fh, 'sha256').hexdigest()


def download(url: str, path: Path) -> None:
    with open(path, 'wb') as fh:
        with httpx.stream("GET", url) as response:
            total = int(response.headers["Content-Length"])
            with tqdm(total=total, unit_scale=True, unit_divisor=1024, unit="B") as progress:
                num_bytes_downloaded = response.num_bytes_downloaded
                for chunk in response.iter_bytes():
                    fh.write(chunk)
                    progress.update(response.num_bytes_downloaded - num_bytes_downloaded)
                    num_bytes_downloaded = response.num_bytes_downloaded


def check_bowtie2_index() -> None:
    if not all(path.exists() for path in BT2_INDEX_PATHS):
        XDG_DATA_DIR.mkdir(exist_ok=True)
        logging.info(f"Fetching human Bowtie2 index")
        download(BT2_ARCHIVE_URL, BT2_ARCHIVE_PATH)
        untar_file(BT2_ARCHIVE_PATH, XDG_DATA_DIR)
        BT2_ARCHIVE_PATH.unlink()
        logging.info(f"Saved Bowtie2 index ({BT2_INDEX_PATH})")
    else:
        logging.info(f"Using cached Bowtie2 index ({BT2_INDEX_PATH})")


def decontaminate_paired_bowtie2(fastq1: Path, fastq2: Path, index: Path = BT2_INDEX_PATH, out_dir: Path = CWD, threads: int = 8):
    """First attempt; non streaming yet faster than naive streaming approach"""
    fastq1, fastq2, out_dir = Path(fastq1), Path(fastq2), Path(out_dir)
    out_dir.mkdir(exist_ok=True)
    fastq1_stem = fastq1.name.removesuffix(fastq1.suffixes[-1]).removesuffix(fastq1.suffixes[-2])
    fastq2_stem = fastq2.name.removesuffix(fastq2.suffixes[-1]).removesuffix(fastq2.suffixes[-2])
    fastq1_out_path = out_dir / f"{fastq1_stem}.dehosted_1.fastq.gz"
    fastq2_out_path = out_dir / f"{fastq2_stem}.dehosted_2.fastq.gz"

    cmd = (f"/Users/bede/Downloads/bowtie2-2.5.1-macos-arm64/bowtie2"
           f" -k 1 -p {threads} -x '{index}'"
           f" -1 '{fastq1}' -2 '{fastq2}'"
           f" | samtools view --threads 4 -f 12 -"
           f' | gawk \'BEGIN{{FS=OFS="\\t"}} {{$1=int((NR+1)/2)" "; print $0}}\''
           f" | samtools fastq --threads 4 -c 6 -N -1 '{fastq1_out_path}' -2 '{fastq2_out_path}'")
    logging.info("Decontaminating reads")
    logging.debug(f"{cmd}")
    run(cmd, cwd=CWD)
    # time hostile --fastq1 tests/data/mtb-jeff/WTCHG_885333_73205296_1.fastq.gz --fastq2 tests/data/mtb-jeff/WTCHG_885333_73205296_2.fastq.gz
    # 21s (200k reads/s)
    logging.info("Generating checksums")
    checksums = {p.name: sha256sum(p) for p in (fastq1, fastq1_out_path, fastq2, fastq2_out_path)}
    logging.info("Complete")
    return checksums


def dehost_fastqs(fastq1: Path, fastq2: Path | None, index: Path = BT2_INDEX_PATH, out_dir: Path = CWD, threads: int = 8) -> dict[str, str]:
    if not fastq2:
        raise NotImplementedError("Hostile currently supports paired reads only")
    checksums = decontaminate_paired_bowtie2(fastq1, fastq2, index=index, out_dir=out_dir)
    return checksums