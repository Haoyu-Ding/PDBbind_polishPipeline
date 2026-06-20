#!/usr/bin/env python3
"""Stage 07: validate outputs and generate dataset reports."""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from pdbbind_pl.complex_builder import analyze_exported_complex
from pdbbind_pl.utils_io import load_simple_yaml, read_parquet_records, write_csv, write_json, write_jsonl, write_parquet


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
    summary = run_validation(PROJECT_ROOT, Path(args.config).resolve())
    print(f"Validated {summary['row_count']} manifest rows")
    print(f"Validated exported complexes: {summary['validated_exported_count']}")
    print(f"Passed validation: {summary['validation_status_counts'].get('ok', 0)}")
    print(f"Parquet: {summary['output_parquet']}")


def run_validation(project_root: Path, config_path: Path) -> dict[str, object]:
    """Validate exported complexes and generate a final QC report."""

    paths_cfg = load_simple_yaml(config_path)
    workspace_cfg = paths_cfg["workspace"]

    rows = read_parquet_records(Path(workspace_cfg["interim_dir"]) / "master_manifest_exported.parquet")
    output_rows: list[dict[str, object]] = [dict(row) for row in rows]

    validated_exported_count = 0
    for row in output_rows:
        if row.get("final_export_status") != "exported":
            row["validation_status"] = "not_exported"
            row["validation_notes"] = None
            continue

        complex_path = Path(str(row["final_complex_pdb_path"]))
        protein_path = Path(str(row["final_protein_chain_pdb_path"]))
        ligand_path = Path(str(row["final_ligand_structure_path"]))
        validation_status, validation_notes = validate_single_export(
            row=row,
            complex_path=complex_path,
            protein_path=protein_path,
            ligand_path=ligand_path,
        )
        row["validation_status"] = validation_status
        row["validation_notes"] = validation_notes
        validated_exported_count += 1

    output_stem = Path(workspace_cfg["interim_dir"]) / "master_manifest_validated"
    field_order = list(output_rows[0].keys()) if output_rows else []
    write_csv(output_stem.with_suffix(".csv"), output_rows, field_order)
    write_jsonl(output_stem.with_suffix(".jsonl"), output_rows)
    write_parquet(output_stem.with_suffix(".parquet"), output_rows, field_order)

    summary = build_validation_summary(output_rows, output_stem, validated_exported_count)
    write_json(Path(workspace_cfg["reports_dir"]) / "master_manifest_validated_summary.json", summary)
    write_markdown_report(Path(workspace_cfg["reports_dir"]) / "dataset_qc_report.md", summary)
    return summary


def validate_single_export(
    row: dict[str, object],
    complex_path: Path,
    protein_path: Path,
    ligand_path: Path,
) -> tuple[str, str | None]:
    """Validate a single exported complex."""

    if not complex_path.exists():
        return "missing_complex_file", "complex file missing"
    if not protein_path.exists():
        return "missing_protein_file", "protein file missing"
    if not ligand_path.exists():
        return "missing_ligand_file", "ligand file missing"

    complex_stats = analyze_exported_complex(complex_path)
    protein_stats = analyze_exported_complex(protein_path)
    ligand_stats = analyze_exported_complex(ligand_path)

    fail_reasons: list[str] = []
    if complex_stats["atom_chain_count"] != 1:
        fail_reasons.append("complex_atom_chain_count_not_1")
    if complex_stats["atom_count"] <= 0:
        fail_reasons.append("complex_missing_protein_atoms")
    if complex_stats["hetatm_count"] <= 0:
        fail_reasons.append("complex_missing_ligand_atoms")
    if not complex_stats["has_ter"]:
        fail_reasons.append("complex_missing_ter")
    if protein_stats["hetatm_count"] != 0:
        fail_reasons.append("protein_file_contains_hetatm")
    if ligand_stats["atom_count"] != 0:
        fail_reasons.append("ligand_file_contains_atom_records")
    if ligand_stats["hetatm_count"] != complex_stats["hetatm_count"]:
        fail_reasons.append("ligand_atom_count_mismatch")

    expected_bucket = row.get("final_dataset_bucket")
    if expected_bucket is None:
        fail_reasons.append("missing_final_bucket")

    if fail_reasons:
        return "failed", ";".join(fail_reasons)
    return "ok", None


def build_validation_summary(
    rows: list[dict[str, object]],
    output_stem: Path,
    validated_exported_count: int,
) -> dict[str, object]:
    """Build a summary of validation results."""

    validation_status_counts: dict[str, int] = defaultdict(int)
    bucket_counts: dict[str, int] = defaultdict(int)
    failed_reason_counts: dict[str, int] = defaultdict(int)

    for row in rows:
        status = str(row.get("validation_status"))
        validation_status_counts[status] += 1
        bucket = row.get("final_dataset_bucket")
        if bucket is not None and not (isinstance(bucket, float) and str(bucket) == "nan"):
            bucket_counts[str(bucket)] += 1
        notes = row.get("validation_notes")
        if isinstance(notes, str) and notes:
            for reason in notes.split(";"):
                failed_reason_counts[reason] += 1

    return {
        "row_count": len(rows),
        "validated_exported_count": validated_exported_count,
        "output_csv": str(output_stem.with_suffix(".csv")),
        "output_jsonl": str(output_stem.with_suffix(".jsonl")),
        "output_parquet": str(output_stem.with_suffix(".parquet")),
        "validation_status_counts": dict(sorted(validation_status_counts.items())),
        "final_dataset_bucket_counts": dict(sorted(bucket_counts.items())),
        "validation_failure_reason_counts": dict(sorted(failed_reason_counts.items())),
    }


def write_markdown_report(path: Path, summary: dict[str, object]) -> None:
    """Write a small human-readable QC report."""

    lines = [
        "# Dataset QC Report",
        "",
        f"- Total manifest rows: {summary['row_count']}",
        f"- Validated exported complexes: {summary['validated_exported_count']}",
        "",
        "## Validation Status Counts",
    ]
    for key, value in summary["validation_status_counts"].items():
        lines.append(f"- {key}: {value}")
    lines.extend(["", "## Final Dataset Bucket Counts"])
    for key, value in summary["final_dataset_bucket_counts"].items():
        lines.append(f"- {key}: {value}")
    failure_counts = summary["validation_failure_reason_counts"]
    if failure_counts:
        lines.extend(["", "## Validation Failure Reasons"])
        for key, value in failure_counts.items():
            lines.append(f"- {key}: {value}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
