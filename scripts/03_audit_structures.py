#!/usr/bin/env python3
"""Stage 03: audit protein and ligand structural constraints."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from pdbbind_pl.structure_audit import run_structure_audit


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
    summary = run_structure_audit(PROJECT_ROOT, Path(args.config).resolve())
    print(f"Audited {summary['row_count']} manifest rows")
    print(f"Hard-filter pass count: {summary['hard_filter_pass_count']}")
    print(f"CSV: {summary['output_csv']}")
    print(f"JSONL: {summary['output_jsonl']}")


if __name__ == "__main__":
    main()
