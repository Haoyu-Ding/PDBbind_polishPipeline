"""Ligand and protein deduplication logic."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
import tarfile
from typing import Any

from rdkit import Chem, DataStructs
from rdkit import RDLogger
from rdkit.Chem import AllChem

from pdbbind_pl.utils_io import load_simple_yaml, read_parquet_records, write_csv, write_json, write_jsonl, write_parquet


def run_deduplication(project_root: Path, config_path: Path) -> dict[str, Any]:
    """Deduplicate nonsugar entries by ligand similarity and exact protein sequence."""

    RDLogger.DisableLog("rdApp.warning")

    paths_cfg = load_simple_yaml(config_path)
    workspace_cfg = paths_cfg["workspace"]
    dataset_cfg = paths_cfg["dataset"]
    filters_cfg = load_simple_yaml(project_root / "config" / "filters.yaml")
    profile = filters_cfg[filters_cfg["profiles"]["active"]]

    rows = read_parquet_records(Path(workspace_cfg["interim_dir"]) / "master_manifest_scored.parquet")
    threshold = float(profile["nonsugar_ligand_dedup"]["similarity_threshold_gt"])
    output_rows = [dict(row) for row in rows]

    with tarfile.open(Path(dataset_cfg["structure_archive"])) as archive:
        fingerprints = build_fingerprints(archive, output_rows)

    assign_ligand_clusters(output_rows, fingerprints, threshold)
    assign_protein_clusters(output_rows)
    assign_final_dataset_buckets(output_rows)

    output_stem = Path(workspace_cfg["interim_dir"]) / "master_manifest_deduped"
    field_order = list(output_rows[0].keys()) if output_rows else []
    write_csv(output_stem.with_suffix(".csv"), output_rows, field_order)
    write_jsonl(output_stem.with_suffix(".jsonl"), output_rows)
    write_parquet(output_stem.with_suffix(".parquet"), output_rows, field_order)
    summary = build_dedup_summary(output_rows, output_stem)
    write_json(Path(workspace_cfg["reports_dir"]) / "master_manifest_deduped_summary.json", summary)
    return summary


def build_fingerprints(archive: tarfile.TarFile, rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Build Morgan fingerprints for nonsugar rows."""

    fingerprints: dict[str, Any] = {}
    for row in rows:
        if row.get("ligand_class") != "nonsugar":
            continue
        pdb_id = str(row["pdb_id"])
        mol = load_ligand_mol(
            archive=archive,
            sdf_member=str(row.get("ligand_sdf_member") or ""),
            mol2_member=str(row.get("ligand_mol2_member") or ""),
        )
        if mol is None:
            continue
        fingerprints[pdb_id] = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
    return fingerprints


def load_ligand_mol(archive: tarfile.TarFile, sdf_member: str, mol2_member: str) -> Chem.Mol | None:
    """Load a ligand molecule from SDF first, then MOL2."""

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


def assign_ligand_clusters(rows: list[dict[str, Any]], fingerprints: dict[str, Any], threshold: float) -> None:
    """Assign ligand similarity clusters within rigid and flexible nonsugar buckets."""

    buckets: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        if row.get("ligand_class") != "nonsugar":
            row["ligand_dedup_cluster_id"] = None
            row["ligand_dedup_is_representative"] = None
            continue
        bucket = str(row.get("nonsugar_flexibility_bucket") or "unbucketed")
        buckets[bucket].append(row)

    for bucket, bucket_rows in buckets.items():
        representatives: list[dict[str, Any]] = []
        cluster_index = 0
        for row in sorted(bucket_rows, key=representative_sort_key):
            pdb_id = str(row["pdb_id"])
            fingerprint = fingerprints.get(pdb_id)
            if fingerprint is None:
                cluster_id = f"{bucket}_ligand_missingfp_{pdb_id}"
                row["ligand_dedup_cluster_id"] = cluster_id
                row["ligand_dedup_is_representative"] = True
                representatives.append(row)
                continue

            assigned = False
            for representative in representatives:
                rep_fp = fingerprints.get(str(representative["pdb_id"]))
                if rep_fp is None:
                    continue
                similarity = DataStructs.TanimotoSimilarity(fingerprint, rep_fp)
                if similarity > threshold:
                    row["ligand_dedup_cluster_id"] = representative["ligand_dedup_cluster_id"]
                    row["ligand_dedup_is_representative"] = False
                    assigned = True
                    break

            if assigned:
                continue

            cluster_index += 1
            cluster_id = f"{bucket}_ligand_{cluster_index:05d}"
            row["ligand_dedup_cluster_id"] = cluster_id
            row["ligand_dedup_is_representative"] = True
            representatives.append(row)


def assign_protein_clusters(rows: list[dict[str, Any]]) -> None:
    """Assign exact sequence protein clusters for ligand representatives."""

    groups: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        if row.get("ligand_class") != "nonsugar":
            row["protein_dedup_cluster_id"] = None
            row["protein_dedup_is_representative"] = None
            continue
        if not row.get("ligand_dedup_is_representative"):
            row["protein_dedup_cluster_id"] = None
            row["protein_dedup_is_representative"] = False
            continue
        bucket = str(row.get("nonsugar_flexibility_bucket") or "unbucketed")
        sequence_hash = str(row.get("protein_sequence_sha1") or f"missing_{row['pdb_id']}")
        groups[(bucket, sequence_hash)].append(row)

    for (bucket, sequence_hash), group_rows in groups.items():
        sorted_rows = sorted(group_rows, key=representative_sort_key)
        cluster_id = f"{bucket}_protein_{sequence_hash[:12]}"
        for index, row in enumerate(sorted_rows):
            row["protein_dedup_cluster_id"] = cluster_id
            row["protein_dedup_is_representative"] = index == 0


def assign_final_dataset_buckets(rows: list[dict[str, Any]]) -> None:
    """Assign final dataset membership buckets."""

    for row in rows:
        if not row.get("hard_filter_pass"):
            row["final_dataset_bucket"] = None
            continue
        if row.get("ligand_class") == "sugar":
            row["final_dataset_bucket"] = "sugar"
            continue
        if row.get("nonsugar_flexibility_bucket") == "intermediate":
            row["final_dataset_bucket"] = None
            continue
        if not row.get("ligand_dedup_is_representative"):
            row["final_dataset_bucket"] = None
            continue
        if not row.get("protein_dedup_is_representative"):
            row["final_dataset_bucket"] = None
            continue
        bucket = str(row.get("nonsugar_flexibility_bucket") or "")
        row["final_dataset_bucket"] = f"nonsugar_{bucket}" if bucket else None


def representative_sort_key(row: dict[str, Any]) -> tuple[float, str]:
    """Sort rows so the best representative is selected first."""

    resolution_value = row.get("resolution_value")
    numeric_resolution = float(resolution_value) if resolution_value is not None else float("inf")
    return (numeric_resolution, str(row["pdb_id"]))


def build_dedup_summary(rows: list[dict[str, Any]], output_stem: Path) -> dict[str, Any]:
    """Build a summary of deduplication results."""

    final_bucket_counts: dict[str, int] = {}
    ligand_rep_count = 0
    protein_rep_count = 0
    for row in rows:
        if row.get("ligand_dedup_is_representative") is True:
            ligand_rep_count += 1
        if row.get("protein_dedup_is_representative") is True:
            protein_rep_count += 1
        bucket = row.get("final_dataset_bucket")
        if bucket is not None:
            final_bucket_counts[str(bucket)] = final_bucket_counts.get(str(bucket), 0) + 1

    return {
        "row_count": len(rows),
        "output_csv": str(output_stem.with_suffix(".csv")),
        "output_jsonl": str(output_stem.with_suffix(".jsonl")),
        "output_parquet": str(output_stem.with_suffix(".parquet")),
        "ligand_representative_count": ligand_rep_count,
        "protein_representative_count": protein_rep_count,
        "final_dataset_bucket_counts": dict(sorted(final_bucket_counts.items())),
    }
