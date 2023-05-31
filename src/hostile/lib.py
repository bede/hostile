import hashlib
import logging
import subprocess
from pathlib import Path

import dnaio
import httpx
import pyfastx
from platformdirs import user_data_dir
from tqdm import tqdm


logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)

CWD = Path.cwd()
XDG_DATA_DIR = Path(user_data_dir("gpas-cli", "GPAS"))
REF_FN = "human.fa.gz"
REF_URL = f"http://178.79.139.243/gpas/{REF_FN}"
REF_PATH = Path(XDG_DATA_DIR / REF_FN)


def run(cmd, cwd=CWD):  # Helper for CLI testing
    return subprocess.run(
        cmd, cwd=cwd, shell=True, check=True, text=True, capture_output=True
)


def sha256sum(filename):
    with open(filename, 'rb', buffering=0) as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


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


def check_ref_exists() -> None:
    if not REF_PATH.exists():
        XDG_DATA_DIR.mkdir(exist_ok=True)
        logging.info(f"Fetching custom human reference sequence")
        download(HUMAN_REF_URL, HUMAN_REF_PATH)
        logging.info(f"Saved {HUMAN_REF_PATH.resolve()}")
 
 # Checksums before and after decontamination

def decontaminate_paired_method_1(fastq1: Path, fastq2: Path, ref: Path = REF_PATH, out_dir: Path = CWD):
    """First attempt; non streaming yet faster than naive streaming approach"""
    fastq1, fastq2, ref, out_dir = Path(fastq1), Path(fastq2), Path(ref), Path(out_dir)
    out_dir.mkdir(exist_ok=True)
    fastq1_stem = fastq1.name.removesuffix(fastq1.suffixes[-1]).removesuffix(fastq1.suffixes[-2])
    fastq2_stem = fastq2.name.removesuffix(fastq2.suffixes[-1]).removesuffix(fastq2.suffixes[-2])
    fastq1_renamed_path = out_dir / f"{fastq1_stem}.renamed_1.fastq.gz"
    fastq2_renamed_path = out_dir / f"{fastq2_stem}.renamed_2.fastq.gz"
    fastq1_out_path = out_dir / f"{fastq1_stem}.dehosted_1.fastq.gz"
    fastq2_out_path = out_dir / f"{fastq2_stem}.dehosted_2.fastq.gz"

    cmd = (f"seqkit replace -p .+ -r '{{nr}} /1' '{fastq1}'" 
           f" | pigz > '{fastq1_renamed_path}' &&"
           f"seqkit replace -p .+ -r '{{nr}} /2' '{fastq2}'"
           f" | pigz > '{fastq2_renamed_path}' &&"
           f"minimap2 -ax sr -m 40 -t 8 '{ref}'"
           f" '{fastq1_renamed_path}' '{fastq2_renamed_path}'"
           f" | samtools view -f 12 -"
           f" | samtools fastq -c 6 -1 '{fastq1_out_path}' -2 '{fastq2_out_path}' &&"
           f" rm '{fastq1_renamed_path}' '{fastq2_renamed_path}'")
    # time hostile --fastq1 tests/data/mtb-jeff/WTCHG_885333_73205296_1.fastq.gz --fastq2 tests/data/mtb-jeff/WTCHG_885333_73205296_2.fastq.gz
    # 120s / 35k / second
    run(cmd, cwd=CWD)


def dehost_fastqs(fastq1: Path, fastq2: Path | None, ref: Path = REF_PATH, out_dir: Path = CWD):
    if not fastq2:
        raise NotImplementedError("Hostile currently supports paired reads only")
    
    check_ref_exists()
    logging.info(f"{fastq1=} {fastq2=}")
    decontaminate_paired_method_1(fastq1, fastq2, out_dir=out_dir, ref=ref)
