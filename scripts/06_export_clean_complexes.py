#!/usr/bin/env python3
"""Stage 06: export cleaned single-chain single-ligand complex PDB files."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from pdbbind_pl.complex_builder import run_complex_export


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
    summary = run_complex_export(PROJECT_ROOT, Path(args.config).resolve())
    print(f"Processed {summary['row_count']} manifest rows")
    print(f"Exported complexes: {summary['exported_count']}")
    print(f"Skipped: {summary['skipped_count']}")
    print(f"Errors: {summary['error_count']}")
    print(f"Parquet: {summary['output_parquet']}")


if __name__ == "__main__":
    main()
