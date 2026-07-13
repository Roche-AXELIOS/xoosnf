#!/usr/bin/env python3

import argparse
import logging
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

__version__ = "latest"


def add_run_info_header(
    input_config: Path,
    output_config: Path,
    run_name: str,
    run_type: str,
) -> None:
    """Prepend a MultiQC report header block with run info to a config file."""
    header_lines = [
        "report_header_info:",
        f"  - Run Name: {run_name}",
        f"  - Run Type: {run_type}",
    ]

    existing = input_config.read_text()

    with output_config.open("w") as out:
        out.write("\n".join(header_lines) + "\n")
        out.write(existing)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Prepend run information to a MultiQC config header.",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=__version__,
    )
    parser.add_argument(
        "--input-config",
        type=Path,
        required=True,
        help="Path to the input MultiQC config YAML file.",
    )
    parser.add_argument(
        "--output-config",
        type=Path,
        required=True,
        help="Path to the output MultiQC config YAML file to write.",
    )
    parser.add_argument(
        "--run-name",
        required=True,
        help="The run name to include in the report header.",
    )
    parser.add_argument(
        "--run-type",
        required=True,
        help="The run type to include in the report header.",
    )
    args = parser.parse_args()

    if not args.input_config.is_file():
        logging.error(f"Input config not found: {args.input_config}")
        sys.exit(1)

    add_run_info_header(
        args.input_config,
        args.output_config,
        args.run_name,
        args.run_type,
    )
