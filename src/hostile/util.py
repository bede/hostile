import concurrent.futures
import subprocess
import tarfile

from functools import partial
from pathlib import Path

import httpx

from tqdm import tqdm


def run(cmd: str, cwd: Path) -> subprocess.CompletedProcess:
    return subprocess.run(
        cmd, shell=True, cwd=cwd, check=True, text=True, capture_output=True
    )


def run_bash(cmd: str, cwd: Path) -> subprocess.CompletedProcess:
    """Needed because /bin/sh does not support process substitution used for tee"""
    return subprocess.run(
        ["/bin/bash", "-c", cmd], cwd=cwd, check=True, text=True, capture_output=True
    )


def run_bash_parallel(
    cmds: dict[str, str], cwd: Path, description: str = "Processing tasks"
) -> dict[str, subprocess.CompletedProcess]:
    with concurrent.futures.ThreadPoolExecutor(max_workers=1) as x:
        futures = {
            x.submit(partial(run_bash, cwd=cwd), cmd): k for k, cmd in cmds.items()
        }
        results = {}
        for future in tqdm(
            concurrent.futures.as_completed(futures),
            total=len(futures),
            desc=description,
        ):
            key = futures[future]
            try:
                results[key] = future.result()
            except Exception as e:
                print(f"Exception occurred during executing command {cmds[key]}: {e}")
        return results


def fastq_path_to_stem(fastq_path: Path) -> str:
    fastq_path = Path(fastq_path)
    return fastq_path.name.removesuffix(fastq_path.suffixes[-1]).removesuffix(
        fastq_path.suffixes[-2]
    )


def parse_count_file(path: Path) -> int:
    with open(path, "r") as fh:
        return int(fh.read().strip())


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
