#!/usr/bin/env python3

import argparse
import math
from enum import Enum
from pathlib import Path
import random
import shutil
import sys
from typing import Protocol
import logging

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

__version__ = "latest"


class RunType(Enum):
    fastq = "fastq"
    basecall_rdb = "basecall_rdb"


class SubsampleStrategy(Protocol):
    def select(self, input_files: list[Path]) -> list[Path]: ...


class FirstFilesStrategy:
    def __init__(self, ratio: float):
        self.ratio = ratio

    def select(self, input_files: list[Path]) -> list[Path]:
        total_files = len(input_files)
        logging.info(f"Total number of files: {total_files}")

        n = int(math.ceil(total_files * self.ratio))
        logging.info(f"Number of files to select based on ratio {self.ratio}: {n}")

        return input_files[:n]


class LastFilesStrategy:
    def __init__(self, ratio: float):
        self.ratio = ratio

    def select(self, input_files: list[Path]) -> list[Path]:
        total_files = len(input_files)
        logging.info(f"Total number of files: {total_files}")

        n = int(math.ceil(total_files * self.ratio))
        logging.info(f"Number of files to select based on ratio {self.ratio}: {n}")

        return input_files[-n:]


class RandomSizeStrategy:
    def __init__(self, ratio: float):
        self.ratio = ratio

    def select(self, input_files: list[Path]) -> list[Path]:
        total_files = len(input_files)
        logging.info(f"Total number of files: {total_files}")

        random.seed(0)
        shuffled = list(input_files)
        random.shuffle(shuffled)

        total_size = sum(f.stat().st_size for f in shuffled)
        logging.info(f"Total size of input files: {total_size} bytes")

        selected: list[Path] = []
        selected_size = 0
        for f in shuffled:
            if total_size == 0 or (selected_size / total_size) >= self.ratio:
                break
            selected.append(f)
            selected_size += f.stat().st_size

        return selected


class RandomAbsoluteSizeStrategy:
    def __init__(self, absolute_size: int):
        self.absolute_size = max(0, absolute_size)

    def select(self, input_files: list[Path]) -> list[Path]:
        total_files = len(input_files)
        logging.info(f"Total number of files: {total_files}")

        random.seed(0)
        shuffled = list(input_files)
        random.shuffle(shuffled)

        target = self.absolute_size
        logging.info(f"Target absolute size: {target} bytes")

        selected: list[Path] = []
        selected_size = 0
        for f in shuffled:
            if selected_size >= target:
                break
            selected.append(f)
            selected_size += f.stat().st_size

        logging.info(f"Selected cumulative size: {selected_size} bytes")

        return selected


def to_output(input_dir: Path, input_file: Path, output_dir: Path) -> Path:
    return output_dir / input_file.relative_to(input_dir)


def create_link(src: Path, tgt: Path, copy_mode: str):
    tgt.parent.mkdir(parents=True, exist_ok=True)
    if copy_mode == "hardlink":
        tgt.hardlink_to(src)
    elif copy_mode == "symlink":
        tgt.symlink_to(src.resolve())
    elif copy_mode == "copy":
        shutil.copy2(src, tgt)


def create_mirror_link(
    input_dir: Path,
    input_file: Path,
    output_dir: Path,
    copy_mode: str,
):
    create_link(input_file, to_output(input_dir, input_file, output_dir), copy_mode)


def subsample_device_fastq(
    input_device_dir: Path,
    output_device_dir: Path,
    copy_mode: str,
    strategy: SubsampleStrategy,
):
    # Find all fastq.gz files in the input directory
    input_files = [Path(p) for p in sorted(input_device_dir.rglob("*.fastq.gz"))]
    selected_files = strategy.select(input_files)

    with open(f"{input_device_dir.name}.selected.txt", "w") as fp:
        fp.writelines([str(x) + "\n" for x in selected_files])

    for input_file in selected_files:
        create_mirror_link(input_device_dir, input_file, output_device_dir, copy_mode)


def subsample_device_basecall_rdb(
    input_device_dir: Path,
    output_device_dir: Path,
    copy_mode: str,
    strategy: SubsampleStrategy,
):

    # Find all Primary Analysis RDB files in the input directory.
    # For legacy versions we also need to copy the main rdb for each device as well.
    # Files should be in <input_device_dir>/device_*/main/
    patterns = ["ExtractMt_*.rdb", "PrimaryAnalysis_*.rdb"]
    input_files = [
        Path(p)
        for p in sorted(
            f for pattern in patterns for f in input_device_dir.rglob(pattern)
        )
    ]
    selected_files = strategy.select(input_files)

    # Legacy RDB files will need a metadata RDB file at device_*/main/main_000001.rdb file
    main_rdb = input_device_dir / "main/main_000001.rdb"
    if main_rdb.exists():
        # If found, we add the symlink
        create_mirror_link(input_device_dir, main_rdb, output_device_dir, copy_mode)

    with open(f"{input_device_dir.name}.selected.txt", "w") as fp:
        fp.writelines([str(x) + "\n" for x in selected_files])

    for input_file in selected_files:
        filename = str(input_file.name)
        if filename.startswith("PrimaryAnalysis_"):
            file_index = filename[len("PrimaryAnalysis_") : -len(".rdb")]
            files_to_link = [
                f"main/PrimaryAnalysis_{file_index}.rdb",
            ]
        elif filename.startswith("ExtractMt_"):
            file_index = filename[len("ExtractMt_") : -len(".rdb")]
            files_to_link = [
                f"main/ExtractMt_{file_index}.rdb",
            ]
        else:
            # we don't mirror other files
            continue

        for fn in files_to_link:
            create_mirror_link(
                input_device_dir,
                input_device_dir / fn,
                output_device_dir,
                copy_mode,
            )


def main(
    input_dir: Path,
    output_dir: Path,
    run_type: RunType,
    copy_mode: str,
    strategy: SubsampleStrategy,
):
    device_dirs = [Path(p) for p in sorted(input_dir.rglob("device_*"))]

    if run_type == RunType.fastq:
        if len(device_dirs) > 0:
            for device_dir in device_dirs:
                subsample_device_fastq(
                    device_dir,
                    to_output(input_dir, device_dir, output_dir),
                    copy_mode,
                    strategy,
                )
            return

        fastq_conversion_dirs = [
            Path(p)
            for p in sorted(input_dir.rglob("fastqConversion"))
            if Path(p).is_dir()
        ]
        if len(fastq_conversion_dirs) > 0:
            for fastq_conversion_dir in fastq_conversion_dirs:
                subsample_device_fastq(
                    fastq_conversion_dir,
                    to_output(input_dir, fastq_conversion_dir, output_dir),
                    copy_mode,
                    strategy,
                )
            return

        raise Exception(
            "No fastq input directories found. Expected either 'device_*' directories "
            "or 'fastqConversion' directories under "
            f"'{input_dir}'."
        )

    if run_type == RunType.basecall_rdb:
        if len(device_dirs) == 0:
            raise Exception(
                f"No device directories found, run directory is likely corrupted: {input_dir}"
            )

        for device_dir in device_dirs:
            subsample_device_basecall_rdb(
                device_dir,
                to_output(input_dir, device_dir, output_dir),
                copy_mode,
                strategy,
            )
        return

    raise ValueError(f"Unrecognized run_type {run_type}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Create a subset of files with links based on ratio or absolute size."
    )
    parser.add_argument(
        "--version",
        action="version",
        version=__version__,
    )
    parser.add_argument(
        "--run-type",
        type=lambda x: RunType[x],
        choices=RunType,
        required=True,
        help="The run dir type: fastq or basecall_rdb.",
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="The input directory containing *.fastq.gz files.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="The output directory where symbolic links will be created.",
    )
    parser.add_argument(
        "--ratio",
        type=float,
        required=False,
        help="Ratio (0-1) of files (count or cumulative size for random). Required for first/last; optional for random.",
    )
    parser.add_argument(
        "--absolute-size",
        type=int,
        required=False,
        help="Absolute cumulative size in bytes for random selection (alternative to --ratio, only valid with --strategy random).",
    )
    parser.add_argument(
        "--copy-mode",
        type=str,
        choices=["symlink", "hardlink", "copy"],
        default="symlink",
        help="The strategy to use for copying files.",
    )
    parser.add_argument(
        "--strategy",
        type=str,
        choices=["first", "last", "random-ratio", "random-absolute-size"],
        default="random-ratio",
        help="The strategy to use for subsampling files.",
    )

    args = parser.parse_args()

    run_type = args.run_type
    input_dir = args.input_dir
    output_dir = args.output_dir
    ratio = args.ratio
    absolute_size = args.absolute_size
    copy_mode = args.copy_mode
    strategy = args.strategy

    if strategy in ("first", "last") and ratio is None:
        logging.error("--ratio is required for strategies 'first' and 'last'.")
        sys.exit(1)

    if strategy == "random-ratio" and ratio is None:
        logging.error("--ratio is required for strategy 'random-ratio'.")
        sys.exit(1)

    if strategy == "random-absolute-size" and absolute_size is None:
        logging.error(
            "--absolute-size is required for strategy 'random-absolute-size'."
        )
        sys.exit(1)

    if ratio is not None and not (0 <= ratio <= 1):
        logging.error("Ratio must be between 0 and 1.")
        sys.exit(1)

    if strategy == "first":
        strategy = FirstFilesStrategy(ratio)
    elif strategy == "last":
        strategy = LastFilesStrategy(ratio)
    elif strategy == "random-ratio":
        strategy = RandomSizeStrategy(ratio)
    elif strategy == "random-absolute-size":
        strategy = RandomAbsoluteSizeStrategy(absolute_size)
    else:
        logging.error(f"Unsupported strategy: {strategy}")
        sys.exit(1)

    main(input_dir, output_dir, run_type, copy_mode, strategy)
