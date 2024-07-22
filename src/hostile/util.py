import concurrent.futures
import gzip
import hashlib
import logging
import multiprocessing
import os
import platform
import subprocess
import tarfile

import dnaio

from pathlib import Path
from platformdirs import user_data_dir

import httpx

from tqdm import tqdm


def choose_default_thread_count(cpu_count: int) -> int:
    """Choose a sensible number of threads for alignment"""
    cpu_count = int(cpu_count)
    if cpu_count <= 1:
        return 1
    elif 1 < cpu_count < 17:
        return int(cpu_count / 2)
    else:
        return 10


CWD = Path.cwd()
CACHE_DIR = (
    Path(os.environ.get("HOSTILE_CACHE_DIR", ""))
    if os.environ.get("HOSTILE_CACHE_DIR")
    else Path(user_data_dir("hostile-eit", "Bede Constantinides"))
)
CPU_COUNT = multiprocessing.cpu_count()
THREADS = choose_default_thread_count(CPU_COUNT)
BUCKET_URL = "https://objectstorage.uk-london-1.oraclecloud.com/n/lr3yhdniv6gu/b/human-genome-indices/o"
DEFAULT_INDEX_NAME = "human-t2t-hla"


def run(cmd: str, cwd: Path | None = None) -> subprocess.CompletedProcess:
    return subprocess.run(
        cmd, shell=True, cwd=cwd, check=True, text=True, capture_output=True
    )


def run_bash(cmd: str, cwd: Path | None = None) -> subprocess.CompletedProcess:
    """Needed because /bin/sh does not support process substitution used for tee"""
    cmd_fmt = f"set -o pipefail; {cmd}"
    return subprocess.run(
        ["/bin/bash", "-c", cmd_fmt],
        cwd=cwd,
        check=True,
        text=True,
        capture_output=True,
    )


def handle_alignment_exceptions(exception: subprocess.CalledProcessError) -> None:
    """Catch samtools view's non-zero exit if all input reads are contaminated"""
    logging.debug(f"stdout: {exception.stdout}")
    logging.debug(f"stderr: {exception.stderr}")
    alignment_successful = False
    stream_empty = False
    if 'Failed to read header for "-"' in exception.stderr:
        stream_empty = True
    if "overall alignment rate" in exception.stderr:  # Bowtie2
        alignment_successful = True
    if "Peak RSS" in exception.stderr:  # Minimap2
        alignment_successful = True
    if alignment_successful and stream_empty:  # Non zero exit but actually fine
        logging.debug("Alignment complete, empty SAM stream, continuing")
        pass
    else:
        logging.error(
            f"Hostile encountered a problem. Check available RAM and storage\n"
            f"pipeline stdout:\n{exception.stdout}\n"
            f"pipeline stderr:\n{exception.stderr}\n"
        )
        raise exception


def run_bash_parallel(
    cmds: list[str], description: str = "Processing"
) -> dict[int, subprocess.CompletedProcess]:
    with concurrent.futures.ThreadPoolExecutor(max_workers=1) as x:
        futures = [x.submit(run_bash, cmd) for cmd in cmds]
        results = {}
        for future in tqdm(
            concurrent.futures.as_completed(futures),
            total=len(futures),
            desc=description,
            disable=len(cmds) == 1,
        ):
            i = futures.index(future)
            try:
                results[i] = future.result()
            except subprocess.CalledProcessError as e:
                handle_alignment_exceptions(e)
        return results


def fastq_path_to_stem(fastq_path: Path) -> str:
    fastq_path = Path(fastq_path)
    stem = fastq_path.name.removesuffix(".gz")
    for suffix in (".fastq", ".fq"):
        stem = stem.removesuffix(suffix)
    return stem


def parse_count_file(path: Path) -> int:
    try:
        with open(path, "r") as fh:
            count = int(fh.read().strip())
    except ValueError:  # file is empty and count is zero
        logging.debug(f"Count file missing: {path}")
        count = 0
    logging.debug(f"{path=} {count=}")
    return count


def fetch_manifest(url: str = BUCKET_URL) -> dict:
    logging.debug("Fetching bucket contents")
    try:
        r = httpx.get(f"{url}/manifest.json")
        r.raise_for_status()
    except httpx.HTTPError:
        raise httpx.HTTPError(
            "Failed to fetch manifest.json from object storage."
            " Ensure you are connected to the internet,"
            " or provide a valid path to a local index"
        )
    return r.json()


def download(url: str, path: Path) -> None:
    try:
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
            response.raise_for_status()
    except httpx.HTTPError:
        raise httpx.HTTPError(
            f"Failed to download {url}."
            f" Ensure you are connected to the internet,"
            f" or provide a valid path to a local index"
        )


def untar_file(input_path, output_path):
    with tarfile.open(input_path) as fh:
        fh.extractall(path=output_path)


def get_platform() -> str:
    return platform.system().lower()


def write_empty_gzip_text_file(path: Path) -> None:
    with gzip.open(path, "wt") as fh:
        fh.write("")


def fix_empty_fastqs(stats) -> None:
    """Find for empty output FASTQs and overwrite them with valid empty gzipped files"""
    for stat in stats:
        if stat.get("reads_out") == 0:
            fastq1_path = stat.get("fastq1_out_path")
            fastq2_path = stat.get("fastq2_out_path")
            if fastq1_path and Path(fastq1_path).is_file():
                write_empty_gzip_text_file(fastq1_path)
            logging.debug(f"Fixing empty fastq: {fastq1_path=}")
            if fastq2_path and Path(fastq2_path).is_file():
                write_empty_gzip_text_file(fastq2_path)
            logging.debug(f"Fixing empty fastq: {fastq2_path=}")


def sha256(file_path: Path) -> str:
    hasher = hashlib.sha256()
    CHUNK_SIZE = 2**24  # 16 MiB
    with open(Path(file_path), "rb") as fh:
        while chunk := fh.read(CHUNK_SIZE):
            hasher.update(chunk)
    return hasher.hexdigest()


def kmerise(in_path: Path, out_path: Path, k: int, step: int) -> Path:
    in_path, out_path = Path(in_path), Path(out_path)
    with dnaio.open(in_path) as reader, dnaio.open(out_path, mode="w") as writer:
        for r in reader:
            for offset in range(0, len(r.sequence) - k + 1, step):
                kmer = r.sequence[offset : offset + k]
                name = r.name.partition(" ")[0]
                kmer_id = f"{name}_{offset}"
                writer.write(dnaio.SequenceRecord(kmer_id, kmer))
    return out_path.absolute()
