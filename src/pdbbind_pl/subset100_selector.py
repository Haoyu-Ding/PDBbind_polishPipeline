"""Selection helpers for a 100-entry nonsugar subset under extra constraints."""

from __future__ import annotations

from pathlib import Path
import shutil
import tarfile
from typing import Any

from pdbbind_pl.subset_selector import pick_ligand_diverse_subset, representative_sort_key
from pdbbind_pl.utils_io import ensure_directory, load_simple_yaml, read_parquet_records, write_csv, write_json, write_jsonl, write_parquet

TARGET_COUNT = 100
EXCLUDED_PDB_IDS = {"1x07"}


def run_subset100_selection(project_root: Path, config_path: Path) -> dict[str, Any]:
    """Select 100 nonsugar entries with MW and protein-length constraints."""

    paths_cfg = load_simple_yaml(config_path)
    workspace_cfg = paths_cfg["workspace"]
    dataset_cfg = paths_cfg["dataset"]

    rows = read_parquet_records(Path(workspace_cfg["interim_dir"]) / "master_manifest_clustered_subset50.parquet")
    candidate_rows = [
        dict(row)
        for row in rows
        if row.get("validation_status") == "ok"
        and row.get("ligand_class") == "nonsugar"
        and row.get("ligand_mol_wt") is not None
        and float(row["ligand_mol_wt"]) > 150.0
        and isinstance(row.get("protein_sequence"), str)
        and len(str(row["protein_sequence"])) <= 300
        and str(row["pdb_id"]).lower() not in EXCLUDED_PDB_IDS
    ]

    with tarfile.open(Path(dataset_cfg["structure_archive"])) as archive:
        selected_rows = select_subset100_candidates(candidate_rows, archive)
    selected_ids = {str(row["pdb_id"]) for row in selected_rows}

    output_rows = []
    for row in rows:
        updated_row = dict(row)
        updated_row["subset100_nonsugar_selected"] = str(updated_row["pdb_id"]) in selected_ids
        output_rows.append(updated_row)

    subset_dir = Path(workspace_cfg["root_dir"]) / "data" / "subsets" / "subset100_nonsugar"
    copy_selected_files(selected_rows, subset_dir)

    output_stem = Path(workspace_cfg["interim_dir"]) / "subset100_nonsugar_selected"
    field_order = list(output_rows[0].keys()) if output_rows else []
    write_csv(output_stem.with_suffix(".csv"), output_rows, field_order)
    write_jsonl(output_stem.with_suffix(".jsonl"), output_rows)
    write_parquet(output_stem.with_suffix(".parquet"), output_rows, field_order)
    summary = build_summary(output_rows, selected_rows, output_stem, subset_dir)
    write_json(Path(workspace_cfg["reports_dir"]) / "subset100_nonsugar_selected_summary.json", summary)
    write_markdown_report(Path(workspace_cfg["reports_dir"]) / "subset100_nonsugar_selected_report.md", summary)
    return summary


def select_subset100_candidates(rows: list[dict[str, Any]], archive: tarfile.TarFile) -> list[dict[str, Any]]:
    """Select candidates with the same cluster-first and ligand-diversity logic as subset50."""

    cluster_groups: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        cluster_id = row.get("protein_cluster_0p7")
        cluster_key = str(cluster_id) if isinstance(cluster_id, str) and cluster_id else f"unclustered_{row['pdb_id']}"
        cluster_groups.setdefault(cluster_key, []).append(row)

    cluster_representatives = [
        sorted(group_rows, key=representative_sort_key)[0]
        for group_rows in cluster_groups.values()
    ]
    if len(cluster_representatives) >= TARGET_COUNT:
        return pick_ligand_diverse_subset(cluster_representatives, TARGET_COUNT, archive=archive)

    selected = list(cluster_representatives)
    selected_ids = {str(row["pdb_id"]) for row in selected}
    remaining = [row for row in sorted(rows, key=representative_sort_key) if str(row["pdb_id"]) not in selected_ids]
    selected.extend(
        pick_ligand_diverse_subset(
            remaining,
            TARGET_COUNT - len(selected),
            archive=archive,
            preselected=selected,
        )
    )
    return selected[:TARGET_COUNT]


def copy_selected_files(selected_rows: list[dict[str, Any]], subset_dir: Path) -> None:
    """Copy exported structures and a complex-only view for the new subset."""

    selected_structures_dir = subset_dir / "selected_structures"
    complex_only_dir = subset_dir / "complex_only"
    if selected_structures_dir.exists():
        shutil.rmtree(selected_structures_dir)
    if complex_only_dir.exists():
        shutil.rmtree(complex_only_dir)
    ensure_directory(selected_structures_dir)
    ensure_directory(complex_only_dir)

    for row in selected_rows:
        pdb_id = str(row["pdb_id"])
        target_dir = selected_structures_dir / pdb_id
        ensure_directory(target_dir)
        for field in ["final_complex_pdb_path", "final_protein_chain_pdb_path", "final_ligand_structure_path"]:
            source = row.get(field)
            if not isinstance(source, str) or not source:
                continue
            source_path = Path(source)
            shutil.copy2(source_path, target_dir / source_path.name)

        complex_source = row.get("final_complex_pdb_path")
        if isinstance(complex_source, str) and complex_source:
            source_path = Path(complex_source)
            shutil.copy2(source_path, complex_only_dir / f"{pdb_id}_complex.pdb")


def build_summary(
    output_rows: list[dict[str, Any]],
    selected_rows: list[dict[str, Any]],
    output_stem: Path,
    subset_dir: Path,
) -> dict[str, Any]:
    """Build summary statistics for the 100-entry nonsugar subset."""

    rigid_count = sum(row.get("nonsugar_flexibility_bucket") == "rigid" for row in selected_rows)
    flexible_count = sum(row.get("nonsugar_flexibility_bucket") == "flexible" for row in selected_rows)
    unique_clusters = len(
        {
            str(row["protein_cluster_0p7"])
            for row in selected_rows
            if isinstance(row.get("protein_cluster_0p7"), str) and row.get("protein_cluster_0p7")
        }
    )

    return {
        "row_count": len(output_rows),
        "selected_count": len(selected_rows),
        "selected_rigid_count": rigid_count,
        "selected_flexible_count": flexible_count,
        "unique_protein_cluster_count": unique_clusters,
        "output_csv": str(output_stem.with_suffix(".csv")),
        "output_jsonl": str(output_stem.with_suffix(".jsonl")),
        "output_parquet": str(output_stem.with_suffix(".parquet")),
        "selected_structures_dir": str(subset_dir / "selected_structures"),
        "complex_only_dir": str(subset_dir / "complex_only"),
    }


def write_markdown_report(path: Path, summary: dict[str, Any]) -> None:
    """Write a small report for the 100-entry nonsugar subset."""

    lines = [
        "# Subset100 Nonsugar Report",
        "",
        f"- Total manifest rows: {summary['row_count']}",
        f"- Selected count: {summary['selected_count']}",
        f"- Selected rigid count: {summary['selected_rigid_count']}",
        f"- Selected flexible count: {summary['selected_flexible_count']}",
        f"- Unique protein cluster count: {summary['unique_protein_cluster_count']}",
        "",
        f"- Selected structures dir: {summary['selected_structures_dir']}",
        f"- Complex-only dir: {summary['complex_only_dir']}",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
