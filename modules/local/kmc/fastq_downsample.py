#!/usr/bin/env python3
import argparse
import logging
import subprocess
from pathlib import Path
from multiprocessing import Pool
from datetime import datetime
from typing import List
from dataclasses import dataclass


@dataclass
class DownsampleResult:
    success: bool
    output_path: str
    original_fastq: Path
    error_message: str = ""


@dataclass
class DownsampleTask:
    fastq: Path
    per_file_read_count: float
    output_dir: Path


def _process_single_fastq(task: DownsampleTask) -> DownsampleResult:
    """
    Downsample a single FASTQ file using seqkit and pigz.
    Args:
        task (DownsampleTask): Task parameters for downsampling.
    Returns:
        DownsampleResult: Result of the downsampling operation.
    """
    output_name = task.fastq.name.removesuffix(".fastq.gz") + ".downsample.fastq.gz"
    output_fastq = task.output_dir / output_name
    try:
        seqkit_cmd = [
            "seqkit", "head",
            "-n", str(int(task.per_file_read_count)),
            str(task.fastq),
        ]
        pigz_cmd = ["pigz", "-p", "1", "--fast", "--stdout"]
        logging.info(
            f"Downsampling {task.fastq}: {' '.join(seqkit_cmd)} | {' '.join(pigz_cmd)} > {output_fastq}"
        )
        with open(output_fastq, "wb") as outfile:
            seqkit_proc = subprocess.Popen(seqkit_cmd, stdout=subprocess.PIPE)
            pigz_proc = subprocess.Popen(
                pigz_cmd, stdin=seqkit_proc.stdout, stdout=outfile
            )
            # Allow seqkit to receive SIGPIPE if pigz exits early
            seqkit_proc.stdout.close()
            pigz_proc.wait()
            seqkit_proc.wait()

        if seqkit_proc.returncode != 0:
            return DownsampleResult(
                False,
                str(output_fastq),
                task.fastq,
                f"seqkit failed with return code {seqkit_proc.returncode}",
            )
        if pigz_proc.returncode != 0:
            return DownsampleResult(
                False,
                str(output_fastq),
                task.fastq,
                f"pigz failed with return code {pigz_proc.returncode}",
            )
        logging.info(f"Finished {task.fastq}.")
        return DownsampleResult(True, str(output_fastq), task.fastq)
    except Exception as exc:
        logging.error(f"Error processing {task.fastq}: {exc}")
        return DownsampleResult(False, str(output_fastq), task.fastq, str(exc))


def fastq_downsample(args: argparse.Namespace) -> None:
    """
    Downsample a list of FASTQ files in parallel.
    Args:
        args (argparse.Namespace): Parsed command-line arguments.
    """
    fastq_files: Path = args.input_file
    read_count: float = args.read_count
    threads: int = int(args.threads)
    output_dir: Path = args.output_dir
    output_file_list: Path = args.output_file_list

    output_dir.mkdir(parents=True, exist_ok=True)

    start_time = datetime.now()
    logging.info(
        f"Starting fastq downsampling at: {start_time.strftime('%Y-%m-%d %H:%M:%S')}"
    )

    with fastq_files.open("r") as f:
        fastqs: List[Path] = [Path(x.strip()) for x in f.readlines() if x.strip()]

    total_fastq_count = len(fastqs)
    if total_fastq_count == 0:
        raise ValueError("No FASTQ files found in the input file.")

    if read_count <= 1:
        raise ValueError(f"Invalid read_count value: {read_count}")

    per_file_read_count = int(read_count / total_fastq_count)
    logging.info(
        f"Total read count: {read_count}, "
        f"total fastq count: {total_fastq_count}, "
        f"per file read count: {per_file_read_count}."
    )

    tasks = [
        DownsampleTask(fastq, per_file_read_count, output_dir)
        for fastq in fastqs
    ]
    num_processes = min(threads, len(tasks))
    logging.info(f"Using {num_processes} processes for parallel downsampling.")

    outputs: List[str] = []
    failed_files: List[DownsampleResult] = []

    with Pool(processes=num_processes) as pool:
        results: List[DownsampleResult] = pool.map(_process_single_fastq, tasks)
        for result in results:
            if result.success:
                outputs.append(result.output_path)
            else:
                failed_files.append(result)

    if failed_files:
        logging.error(f"Failed to process {len(failed_files)} files.")
        for res in failed_files:
            logging.error(f"  {res.original_fastq}: {res.error_message}")
        raise RuntimeError("One or more FASTQ files failed to downsample.")

    with output_file_list.open("w") as f:
        f.write("\n".join(outputs) + "\n")

    end_time = datetime.now()
    duration = end_time - start_time
    logging.info(
        f"Finished fastq downsampling at: {end_time.strftime('%Y-%m-%d %H:%M:%S')}"
    )
    logging.info(f"Total processing time: {duration}")
    logging.info(
        f"Output file list written to {output_file_list} with {len(outputs)} downsampled file paths."
    )


def main() -> None:
    """
    Entry point for the fastq_downsample script.
    """
    logging.basicConfig(
        level=logging.INFO, format="%(filename)s-%(levelname)s: %(message)s"
    )
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-file", type=Path, required=True)
    parser.add_argument("--read-count", type=float, required=False, default=0)
    parser.add_argument("--threads", type=int, required=False, default=16)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--output-file-list", type=Path, required=True)
    args = parser.parse_args()
    try:
        fastq_downsample(args)
    except Exception as exc:
        logging.error(f"Downsampling failed: {exc}")
        exit(1)


if __name__ == "__main__":
    main()
