"""Subset selection helpers for 50 representatives per final bucket."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
import shutil
import tarfile
from typing import Any

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit import RDLogger

from pdbbind_pl.utils_io import ensure_directory, load_simple_yaml, read_parquet_records, write_csv, write_json, write_jsonl, write_parquet

BUCKETS = ["sugar", "nonsugar_rigid", "nonsugar_flexible"]
TARGET_PER_BUCKET = 50


def run_subset_selection(project_root: Path, config_path: Path) -> dict[str, Any]:
    """Select 50 representatives from each final dataset bucket."""

    RDLogger.DisableLog("rdApp.warning")

    paths_cfg = load_simple_yaml(config_path)
    workspace_cfg = paths_cfg["workspace"]
    dataset_cfg = paths_cfg["dataset"]

    rows = read_parquet_records(Path(workspace_cfg["interim_dir"]) / "master_manifest_clustered_subset50.parquet")
    selected_ids: dict[str, set[str]] = {bucket: set() for bucket in BUCKETS}
    structure_archive = Path(dataset_cfg["structure_archive"])

    with tarfile.open(structure_archive) as archive:
        for bucket in BUCKETS:
            bucket_rows = [
                dict(row)
                for row in rows
                if row.get("final_dataset_bucket") == bucket and row.get("validation_status") == "ok"
            ]
            selected = select_bucket_rows(bucket_rows, bucket, archive)
            selected_ids[bucket] = {str(row["pdb_id"]) for row in selected}

    output_rows = []
    for row in rows:
        updated_row = dict(row)
        bucket = updated_row.get("final_dataset_bucket")
        updated_row["subset50_selected"] = False
        if isinstance(bucket, str) and str(updated_row["pdb_id"]) in selected_ids.get(bucket, set()):
            updated_row["subset50_selected"] = True
        output_rows.append(updated_row)

    subset_dir = Path(workspace_cfg["root_dir"]) / "data" / "subsets" / "subset50"
    copy_selected_files(output_rows, subset_dir)
    copy_complex_only_files(output_rows, subset_dir)

    output_stem = Path(workspace_cfg["interim_dir"]) / "subset50_selected"
    field_order = list(output_rows[0].keys()) if output_rows else []
    write_csv(output_stem.with_suffix(".csv"), output_rows, field_order)
    write_jsonl(output_stem.with_suffix(".jsonl"), output_rows)
    write_parquet(output_stem.with_suffix(".parquet"), output_rows, field_order)
    summary = build_subset_summary(output_rows, output_stem, subset_dir)
    write_json(Path(workspace_cfg["reports_dir"]) / "subset50_selected_summary.json", summary)
    write_markdown_report(Path(workspace_cfg["reports_dir"]) / "subset50_selected_report.md", summary)
    return summary


def select_bucket_rows(rows: list[dict[str, Any]], bucket: str, archive: tarfile.TarFile) -> list[dict[str, Any]]:
    """Select up to 50 rows from one bucket with protein cluster and ligand diversity constraints."""

    cluster_groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        cluster_id = row.get("protein_cluster_0p7")
        cluster_key = str(cluster_id) if isinstance(cluster_id, str) and cluster_id else f"unclustered_{row['pdb_id']}"
        cluster_groups[cluster_key].append(row)

    cluster_representatives = [
        sorted(group_rows, key=representative_sort_key)[0]
        for group_rows in cluster_groups.values()
    ]
    if len(cluster_representatives) >= TARGET_PER_BUCKET:
        return pick_ligand_diverse_subset(cluster_representatives, TARGET_PER_BUCKET, archive=archive)

    selected = list(cluster_representatives)
    selected_ids = {str(row["pdb_id"]) for row in selected}
    remaining = [row for row in sorted(rows, key=representative_sort_key) if str(row["pdb_id"]) not in selected_ids]
    selected.extend(
        pick_ligand_diverse_subset(
            remaining,
            TARGET_PER_BUCKET - len(selected),
            archive=archive,
            preselected=selected,
        )
    )
    return selected[:TARGET_PER_BUCKET]


def representative_sort_key(row: dict[str, Any]) -> tuple[float, str]:
    """Prefer better structural quality first."""

    resolution_value = row.get("resolution_value")
    resolution = float(resolution_value) if resolution_value is not None else float("inf")
    return (resolution, str(row["pdb_id"]))


def pick_ligand_diverse_subset(
    candidate_rows: list[dict[str, Any]],
    target_count: int,
    archive: tarfile.TarFile,
    preselected: list[dict[str, Any]] | None = None,
) -> list[dict[str, Any]]:
    """Greedy max-min selection by ligand fingerprint diversity with quality-aware seeding."""

    if target_count <= 0:
        return []

    selected = list(preselected or [])
    selected_ids = {str(row["pdb_id"]) for row in selected}
    candidates = [row for row in sorted(candidate_rows, key=representative_sort_key) if str(row["pdb_id"]) not in selected_ids]
    if not candidates:
        return [] if preselected is None else selected

    fingerprints = {str(row["pdb_id"]): load_ligand_fingerprint(row, archive) for row in candidates + selected}
    if preselected is None and candidates:
        selected.append(candidates.pop(0))
        selected_ids.add(str(selected[-1]["pdb_id"]))

    while len(selected) - len(preselected or []) < target_count and candidates:
        best_row = None
        best_score = None
        for row in candidates:
            row_fp = fingerprints.get(str(row["pdb_id"]))
            if row_fp is None:
                min_distance = 0.0
            else:
                distances = []
                for chosen in selected:
                    chosen_fp = fingerprints.get(str(chosen["pdb_id"]))
                    if chosen_fp is None:
                        continue
                    similarity = DataStructs.TanimotoSimilarity(row_fp, chosen_fp)
                    distances.append(1.0 - similarity)
                min_distance = min(distances) if distances else 1.0
            score = (min_distance, -representative_sort_key(row)[0], str(row["pdb_id"]))
            if best_score is None or score > best_score:
                best_score = score
                best_row = row
        if best_row is None:
            break
        selected.append(best_row)
        selected_ids.add(str(best_row["pdb_id"]))
        candidates = [row for row in candidates if str(row["pdb_id"]) != str(best_row["pdb_id"])]

    if preselected is None:
        return selected[:target_count]
    return selected


def load_ligand_fingerprint(row: dict[str, Any], archive: tarfile.TarFile):
    """Build a Morgan fingerprint from the original ligand structure files."""

    sdf_member = row.get("ligand_sdf_member")
    mol2_member = row.get("ligand_mol2_member")
    mol = load_ligand_mol_from_archive(
        archive,
        str(sdf_member) if isinstance(sdf_member, str) else "",
        str(mol2_member) if isinstance(mol2_member, str) else "",
    )
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)


def load_ligand_mol_from_archive(archive: tarfile.TarFile, sdf_member: str, mol2_member: str) -> Chem.Mol | None:
    """Load a ligand molecule from the original archive."""

    if sdf_member:
        handle = archive.extractfile(sdf_member)
        if handle is not None:
            block = handle.read().decode("utf-8", "ignore")
            mol = Chem.MolFromMolBlock(block, removeHs=False, sanitize=True)
            if mol is not None:
                return mol
    if mol2_member:
        handle = archive.extractfile(mol2_member)
        if handle is not None:
            block = handle.read().decode("utf-8", "ignore")
            mol = Chem.MolFromMol2Block(block, removeHs=False, sanitize=True)
            if mol is not None:
                return mol
    return None


def copy_selected_files(rows: list[dict[str, Any]], subset_dir: Path) -> None:
    """Copy selected exported structures into a dedicated subset directory."""

    final_dir = subset_dir / "selected_structures"
    complex_only_dir = subset_dir / "complex_only"
    if final_dir.exists():
        shutil.rmtree(final_dir)
    if complex_only_dir.exists():
        shutil.rmtree(complex_only_dir)
    ensure_directory(final_dir)
    for row in rows:
        if not row.get("subset50_selected"):
            continue
        bucket = str(row["final_dataset_bucket"])
        pdb_id = str(row["pdb_id"])
        target_dir = final_dir / bucket / pdb_id
        ensure_directory(target_dir)
        for field in ["final_complex_pdb_path", "final_protein_chain_pdb_path", "final_ligand_structure_path"]:
            source = row.get(field)
            if not isinstance(source, str) or not source:
                continue
            source_path = Path(source)
            shutil.copy2(source_path, target_dir / source_path.name)


def build_subset_summary(rows: list[dict[str, Any]], output_stem: Path, subset_dir: Path) -> dict[str, Any]:
    """Build a summary for the 50x3 subset selection."""

    counts: dict[str, int] = defaultdict(int)
    cluster_counts: dict[str, int] = defaultdict(int)
    for row in rows:
        if not row.get("subset50_selected"):
            continue
        bucket = str(row["final_dataset_bucket"])
        counts[bucket] += 1
        cluster_id = row.get("protein_cluster_0p7")
        if isinstance(cluster_id, str) and cluster_id:
            cluster_counts[f"{bucket}:{cluster_id}"] += 1

    unique_clusters_by_bucket: dict[str, int] = defaultdict(int)
    for key in cluster_counts:
        bucket, _cluster_id = key.split(":", 1)
        unique_clusters_by_bucket[bucket] += 1

    return {
        "row_count": len(rows),
        "output_csv": str(output_stem.with_suffix(".csv")),
        "output_jsonl": str(output_stem.with_suffix(".jsonl")),
        "output_parquet": str(output_stem.with_suffix(".parquet")),
        "selected_structures_dir": str(subset_dir / "selected_structures"),
        "complex_only_dir": str(subset_dir / "complex_only"),
        "selected_counts": dict(sorted(counts.items())),
        "unique_protein_clusters_by_bucket": dict(sorted(unique_clusters_by_bucket.items())),
    }


def write_markdown_report(path: Path, summary: dict[str, Any]) -> None:
    """Write a small report for the subset selection."""

    lines = [
        "# Subset50 Report",
        "",
        f"- Total manifest rows: {summary['row_count']}",
        "",
        "## Selected Counts",
    ]
    for key, value in summary["selected_counts"].items():
        lines.append(f"- {key}: {value}")
    lines.extend(["", "## Unique Protein Clusters By Bucket"])
    for key, value in summary["unique_protein_clusters_by_bucket"].items():
        lines.append(f"- {key}: {value}")
    lines.extend(
        [
            "",
            f"- Selected structures dir: {summary['selected_structures_dir']}",
            f"- Complex-only dir: {summary['complex_only_dir']}",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def copy_complex_only_files(rows: list[dict[str, Any]], subset_dir: Path) -> None:
    """Copy only complex PDB files into flat per-bucket directories for easy browsing."""

    complex_only_dir = subset_dir / "complex_only"
    ensure_directory(complex_only_dir)
    for bucket in BUCKETS:
        ensure_directory(complex_only_dir / bucket)

    for row in rows:
        if not row.get("subset50_selected"):
            continue
        bucket = row.get("final_dataset_bucket")
        source = row.get("final_complex_pdb_path")
        pdb_id = row.get("pdb_id")
        if not isinstance(bucket, str) or not isinstance(source, str) or not isinstance(pdb_id, str):
            continue
        source_path = Path(source)
        target_path = complex_only_dir / bucket / f"{pdb_id}_complex.pdb"
        shutil.copy2(source_path, target_path)
