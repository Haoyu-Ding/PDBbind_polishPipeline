#!/usr/bin/env python3
"""Stage 05: deduplicate nonsugar entries by ligand and protein similarity."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from pdbbind_pl.dedup import run_deduplication


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
    summary = run_deduplication(PROJECT_ROOT, Path(args.config).resolve())
    print(f"Deduplicated {summary['row_count']} manifest rows")
    print(f"Parquet: {summary['output_parquet']}")
    print(f"Final dataset buckets: {summary['final_dataset_bucket_counts']}")


if __name__ == "__main__":
    main()
