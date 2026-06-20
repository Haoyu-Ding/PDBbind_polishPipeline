#!/usr/bin/env python3
"""Stage 10: select a 100-entry nonsugar subset with MW and protein-length constraints."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from pdbbind_pl.subset100_selector import run_subset100_selection


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
    summary = run_subset100_selection(PROJECT_ROOT, Path(args.config).resolve())
    print(f"Selected {summary['selected_count']} nonsugar entries")
    print(f"Rigid: {summary['selected_rigid_count']}")
    print(f"Flexible: {summary['selected_flexible_count']}")
    print(f"Unique protein clusters: {summary['unique_protein_cluster_count']}")
    print(f"Parquet: {summary['output_parquet']}")


if __name__ == "__main__":
    main()
