import concurrent.futures
import logging
import subprocess
import tarfile

from pathlib import Path

import httpx

from tqdm import tqdm


def run(cmd: str, cwd: Path | None = None) -> subprocess.CompletedProcess:
    return subprocess.run(
        cmd, shell=True, cwd=cwd, check=True, text=True, capture_output=True
    )


def run_bash(cmd: str, cwd: Path | None = None) -> subprocess.CompletedProcess:
    """Needed because /bin/sh does not support process substitution used for tee"""
    return subprocess.run(
        ["/bin/bash", "-c", cmd], cwd=cwd, check=True, text=True, capture_output=True
    )


def handle_alignment_exceptions(exception: subprocess.CalledProcessError) -> None:
    """Catch samtools view's non-zero exit if all input reads are contaminated"""
    alignment_successful = False
    stream_empty = False
    if 'Failed to read header for "-"' in exception.stderr:
        stream_empty = True
    if "overall alignment rate" in exception.stderr:  # Bowtie2
        alignment_successful = True
    if "Peak RSS" in exception.stderr:  # Minimap2
        alignment_successful = True
    logging.debug(f"{stream_empty=} {alignment_successful=}")
    if alignment_successful and stream_empty:  # Non zero exit but actually fine
        pass
    else:
        print(f"Hostile encountered a problem. Stderr below")
        print(f"{exception.stderr}")
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


def untar_file(input_path, output_path):
    with tarfile.open(input_path) as fh:
        fh.extractall(path=output_path)


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
