#!/usr/bin/env python3
"""Stage 09: select 50 entries from each final bucket with protein-cluster constraints."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from pdbbind_pl.subset_selector import run_subset_selection


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        default=str(PROJECT_ROOT / "config" / "paths.yaml"),
        help="Path to the pipeline paths.yaml configuration file.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    summary = run_subset_selection(PROJECT_ROOT, Path(args.config).resolve())
    print(f"Selected subset rows from {summary['row_count']} manifest rows")
    print(f"Counts: {summary['selected_counts']}")
    print(f"Parquet: {summary['output_parquet']}")


if __name__ == "__main__":
    main()
