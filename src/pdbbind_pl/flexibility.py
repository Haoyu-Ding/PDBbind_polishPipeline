"""Nonsugar ligand flexibility scoring logic."""

from __future__ import annotations

from pathlib import Path
import tarfile
from typing import Any

from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import Lipinski, rdMolDescriptors

from pdbbind_pl.utils_io import load_simple_yaml, read_parquet_records, write_csv, write_json, write_jsonl, write_parquet


def run_flexibility_scoring(project_root: Path, config_path: Path) -> dict[str, Any]:
    """Score nonsugar ligands for flexibility and update the manifest."""

    RDLogger.DisableLog("rdApp.warning")

    paths_cfg = load_simple_yaml(config_path)
    workspace_cfg = paths_cfg["workspace"]
    dataset_cfg = paths_cfg["dataset"]
    filters_cfg = load_simple_yaml(project_root / "config" / "filters.yaml")
    profile = filters_cfg[filters_cfg["profiles"]["active"]]
    flex_rules = profile["nonsugar_flexibility"]

    rows = read_parquet_records(Path(workspace_cfg["interim_dir"]) / "master_manifest_classified.parquet")
    output_rows: list[dict[str, Any]] = []

    with tarfile.open(Path(dataset_cfg["structure_archive"])) as archive:
        for row in rows:
            updated_row = dict(row)
            if updated_row.get("ligand_class") != "nonsugar":
                output_rows.append(updated_row)
                continue

            mol = load_ligand_mol(
                archive=archive,
                sdf_member=str(updated_row.get("ligand_sdf_member") or ""),
                mol2_member=str(updated_row.get("ligand_mol2_member") or ""),
            )
            if mol is None:
                updated_row["ligand_parse_status"] = "failed_all"
                updated_row["ligand_rotatable_bonds"] = None
                updated_row["ligand_rotatable_bonds_ratio"] = None
                updated_row["ligand_formula"] = None
                updated_row["ligand_mol_wt"] = None
                updated_row["nonsugar_flexibility_bucket"] = None
                output_rows.append(updated_row)
                continue

            rotatable_bonds = int(Lipinski.NumRotatableBonds(mol))
            heavy_atom_count = int(updated_row.get("ligand_heavy_atom_count") or mol.GetNumHeavyAtoms() or 1)
            rotatable_ratio = rotatable_bonds / heavy_atom_count

            updated_row["ligand_rotatable_bonds"] = rotatable_bonds
            updated_row["ligand_rotatable_bonds_ratio"] = rotatable_ratio
            updated_row["ligand_formula"] = rdMolDescriptors.CalcMolFormula(mol)
            updated_row["ligand_mol_wt"] = float(rdMolDescriptors.CalcExactMolWt(mol))
            updated_row["nonsugar_flexibility_bucket"] = classify_flexibility_bucket(
                rotatable_bonds=rotatable_bonds,
                rotatable_ratio=rotatable_ratio,
                rules=flex_rules,
            )
            output_rows.append(updated_row)

    output_stem = Path(workspace_cfg["interim_dir"]) / "master_manifest_scored"
    field_order = list(output_rows[0].keys()) if output_rows else []
    write_csv(output_stem.with_suffix(".csv"), output_rows, field_order)
    write_jsonl(output_stem.with_suffix(".jsonl"), output_rows)
    write_parquet(output_stem.with_suffix(".parquet"), output_rows, field_order)
    summary = build_flexibility_summary(output_rows, output_stem)
    write_json(Path(workspace_cfg["reports_dir"]) / "master_manifest_scored_summary.json", summary)
    return summary


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


def classify_flexibility_bucket(rotatable_bonds: int, rotatable_ratio: float, rules: dict[str, Any]) -> str:
    """Assign a V2 flexibility bucket using absolute RB and RB ratio."""

    if (
        rotatable_bonds <= int(rules["rigid_max_rotatable_bonds"])
        and rotatable_ratio <= float(rules["rigid_max_rotatable_ratio"])
    ):
        return "rigid"
    if (
        rotatable_bonds >= int(rules["flexible_min_rotatable_bonds"])
        or rotatable_ratio >= float(rules["flexible_min_rotatable_ratio"])
    ):
        return "flexible"
    return "intermediate"


def build_flexibility_summary(rows: list[dict[str, Any]], output_stem: Path) -> dict[str, Any]:
    """Build a summary of nonsugar flexibility scoring."""

    bucket_counts: dict[str, int] = {}
    for row in rows:
        bucket = row.get("nonsugar_flexibility_bucket")
        if bucket is None:
            continue
        bucket_counts[str(bucket)] = bucket_counts.get(str(bucket), 0) + 1

    return {
        "row_count": len(rows),
        "output_csv": str(output_stem.with_suffix(".csv")),
        "output_jsonl": str(output_stem.with_suffix(".jsonl")),
        "output_parquet": str(output_stem.with_suffix(".parquet")),
        "nonsugar_flexibility_bucket_counts": dict(sorted(bucket_counts.items())),
    }
