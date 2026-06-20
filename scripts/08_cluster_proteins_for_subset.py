#!/usr/bin/env python3
"""Stage 08: cluster validated exported proteins with MMseqs2 for subset selection."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from pdbbind_pl.protein_cluster import run_protein_clustering


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
    summary = run_protein_clustering(PROJECT_ROOT, Path(args.config).resolve())
    print(f"Clustered {summary['assigned_cluster_count']} validated entries")
    print(f"Unique clusters: {summary['unique_cluster_count']}")
    print(f"Parquet: {summary['output_parquet']}")


if __name__ == "__main__":
    main()
