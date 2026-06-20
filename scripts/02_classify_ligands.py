#!/usr/bin/env python3
"""Stage 02: classify ligands into sugar, nonsugar, or ambiguous."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from pdbbind_pl.ligand_classifier import run_ligand_classification


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
    summary = run_ligand_classification(PROJECT_ROOT, Path(args.config).resolve())
    print(f"Classified {summary['row_count']} ligands")
    print(f"Parquet: {summary['output_parquet']}")
    print(f"Class counts: {summary['ligand_class_counts']}")


if __name__ == "__main__":
    main()
