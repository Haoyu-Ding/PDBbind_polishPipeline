#!/usr/bin/env python3
"""Stage 00: extract or normalize the PDBbind protein-ligand index."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from pdbbind_pl.index_parser import INDEX_FIELD_ORDER, load_index_records
from pdbbind_pl.utils_io import load_simple_yaml, write_csv, write_json, write_jsonl, write_parquet


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
    config_path = Path(args.config).resolve()
    config = load_simple_yaml(config_path)

    dataset_cfg = config["dataset"]
    workspace_cfg = config["workspace"]

    records = load_index_records(
        index_archive=Path(dataset_cfg["index_archive"]),
        index_member=dataset_cfg["index_member"],
    )

    interim_dir = Path(workspace_cfg["interim_dir"])
    reports_dir = Path(workspace_cfg["reports_dir"])
    output_stem = interim_dir / "index_general_pl_normalized"

    write_csv(output_stem.with_suffix(".csv"), records, INDEX_FIELD_ORDER)
    write_jsonl(output_stem.with_suffix(".jsonl"), records)
    write_parquet(output_stem.with_suffix(".parquet"), records, INDEX_FIELD_ORDER)
    write_json(
        reports_dir / "index_general_pl_summary.json",
        build_summary(records, dataset_cfg["index_member"], output_stem),
    )

    print(f"Wrote {len(records)} normalized index rows")
    print(f"CSV: {output_stem.with_suffix('.csv')}")
    print(f"JSONL: {output_stem.with_suffix('.jsonl')}")
    print(f"Parquet: {output_stem.with_suffix('.parquet')}")


def build_summary(
    records: list[dict[str, object]],
    index_member: str,
    output_stem: Path,
) -> dict[str, object]:
    """Build a compact extraction summary for quick inspection."""

    return {
        "row_count": len(records),
        "source_index_member": index_member,
        "output_csv": str(output_stem.with_suffix(".csv")),
        "output_jsonl": str(output_stem.with_suffix(".jsonl")),
        "output_parquet": str(output_stem.with_suffix(".parquet")),
        "flag_counts": {
            "covalent_complex": sum(bool(record["flag_covalent_complex"]) for record in records),
            "incomplete_ligand": sum(bool(record["flag_incomplete_ligand"]) for record in records),
            "isomer_annotation": sum(bool(record["flag_isomer_annotation"]) for record in records),
            "redundant_annotation": sum(bool(record["flag_redundant_annotation"]) for record in records),
            "peptide_like_mer_annotation": sum(
                bool(record["flag_peptide_like_mer_annotation"]) for record in records
            ),
        },
        "method_counts": count_by(records, "structure_method_class"),
    }


def count_by(records: list[dict[str, object]], field: str) -> dict[str, int]:
    """Count categorical values in a list of records."""

    counts: dict[str, int] = {}
    for record in records:
        key = str(record[field])
        counts[key] = counts.get(key, 0) + 1
    return counts


if __name__ == "__main__":
    main()
