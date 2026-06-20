"""Ligand sugar vs nonsugar classification logic."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import tarfile
from typing import Any

from rdkit import Chem
from rdkit import RDLogger

from pdbbind_pl.utils_io import load_simple_yaml, read_parquet_records, write_csv, write_json, write_jsonl, write_parquet


@dataclass
class SugarClassificationResult:
    """Final sugar classification payload."""

    ligand_class: str
    ligand_class_reason: str


def run_ligand_classification(project_root: Path, config_path: Path) -> dict[str, Any]:
    """Classify ligands into sugar or nonsugar under the V2 rules."""

    RDLogger.DisableLog("rdApp.warning")

    paths_cfg = load_simple_yaml(config_path)
    workspace_cfg = paths_cfg["workspace"]
    dataset_cfg = paths_cfg["dataset"]
    sugar_rules_cfg = load_simple_yaml(project_root / "config" / "sugar_rules.yaml")

    rows = read_parquet_records(Path(workspace_cfg["interim_dir"]) / "master_manifest_audited.parquet")
    output_rows: list[dict[str, Any]] = []

    with tarfile.open(Path(dataset_cfg["structure_archive"])) as archive:
        for row in rows:
            updated_row = dict(row)
            classification = classify_ligand_row(archive, updated_row, sugar_rules_cfg)
            updated_row["ligand_class"] = classification.ligand_class
            updated_row["ligand_class_reason"] = classification.ligand_class_reason
            output_rows.append(updated_row)

    output_stem = Path(workspace_cfg["interim_dir"]) / "master_manifest_classified"
    field_order = list(output_rows[0].keys()) if output_rows else []
    write_csv(output_stem.with_suffix(".csv"), output_rows, field_order)
    write_jsonl(output_stem.with_suffix(".jsonl"), output_rows)
    write_parquet(output_stem.with_suffix(".parquet"), output_rows, field_order)
    summary = build_classification_summary(output_rows, output_stem)
    write_json(Path(workspace_cfg["reports_dir"]) / "master_manifest_classified_summary.json", summary)
    return summary


def classify_ligand_row(
    archive: tarfile.TarFile,
    row: dict[str, Any],
    sugar_rules_cfg: dict[str, Any],
) -> SugarClassificationResult:
    """Classify a single ligand row."""

    if row.get("ligand_parse_status") != "ok":
        return SugarClassificationResult("nonsugar", "parse_failure_defaults_to_nonsugar")
    if row.get("ligand_component_count") != 1:
        return SugarClassificationResult("nonsugar", "multi_component_defaults_to_nonsugar")

    mol = load_ligand_mol(
        archive=archive,
        sdf_member=str(row.get("ligand_sdf_member") or ""),
        mol2_member=str(row.get("ligand_mol2_member") or ""),
    )
    if mol is None:
        return SugarClassificationResult("nonsugar", "rdkit_parse_failure")

    heuristics = sugar_rules_cfg["heuristics"]
    sugar_signal, reason = detect_sugar_signal(mol, heuristics)
    if sugar_signal:
        return SugarClassificationResult("sugar", reason)
    return SugarClassificationResult("nonsugar", reason)


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


def detect_sugar_signal(mol: Chem.Mol, heuristics: dict[str, Any]) -> tuple[bool, str]:
    """Return a sugar decision and machine-readable reason."""

    candidate_rings = []
    for ring in mol.GetRingInfo().AtomRings():
        ring_size = len(ring)
        if heuristics["require_5_or_6_member_ring"] and ring_size not in {5, 6}:
            continue
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        ring_oxygen_count = sum(atom.GetAtomicNum() == 8 for atom in ring_atoms)
        if heuristics["require_ring_oxygen_for_primary_sugar_call"] and ring_oxygen_count < int(
            heuristics["min_ring_oxygen_count"]
        ):
            continue
        candidate_rings.append((ring, ring_oxygen_count))

    if not candidate_rings:
        return False, "no_5_or_6_member_ring_with_ring_oxygen"

    oxygen_count = sum(atom.GetAtomicNum() == 8 for atom in mol.GetAtoms())
    carbon_count = sum(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
    if carbon_count == 0:
        return False, "no_carbon_backbone"
    oxygen_to_carbon_ratio = oxygen_count / carbon_count
    if oxygen_to_carbon_ratio < float(heuristics["min_oxygen_to_carbon_ratio"]):
        return False, "oxygen_to_carbon_ratio_too_low"

    exocyclic_oxygen_count = count_exocyclic_oxygen_atoms(mol)
    if exocyclic_oxygen_count < int(heuristics["min_exocyclic_oxygen_count"]):
        return False, "too_few_exocyclic_oxygens"

    ring_atom_indices = set()
    for ring, _ring_oxygen_count in candidate_rings:
        ring_atom_indices.update(ring)
    sugar_like_atom_indices = expand_ring_environment(mol, ring_atom_indices)
    heavy_atom_indices = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1}
    if not heavy_atom_indices:
        return False, "no_heavy_atoms"
    non_sugar_fraction = (len(heavy_atom_indices - sugar_like_atom_indices) / len(heavy_atom_indices))
    if non_sugar_fraction > float(heuristics["max_non_sugar_heavy_atom_fraction"]):
        return False, "non_sugar_scaffold_dominates"

    return True, "ring_oxygen_exocyclic_oxygen_ratio_match"


def count_exocyclic_oxygen_atoms(mol: Chem.Mol) -> int:
    """Count oxygen atoms directly attached to ring atoms while not being in a ring."""

    count = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8 or atom.IsInRing():
            continue
        if any(neighbor.IsInRing() for neighbor in atom.GetNeighbors()):
            count += 1
    return count


def expand_ring_environment(mol: Chem.Mol, ring_atom_indices: set[int]) -> set[int]:
    """Expand the sugar-like environment to ring atoms and their directly attached oxygens/carbons."""

    sugar_like = set(ring_atom_indices)
    for atom_idx in list(ring_atom_indices):
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() in {6, 8}:
                sugar_like.add(neighbor.GetIdx())
    return sugar_like


def build_classification_summary(rows: list[dict[str, Any]], output_stem: Path) -> dict[str, Any]:
    """Build a compact classification summary."""

    counts: dict[str, int] = {}
    reason_counts: dict[str, int] = {}
    for row in rows:
        ligand_class = str(row.get("ligand_class"))
        counts[ligand_class] = counts.get(ligand_class, 0) + 1
        reason = str(row.get("ligand_class_reason"))
        reason_counts[reason] = reason_counts.get(reason, 0) + 1

    return {
        "row_count": len(rows),
        "output_csv": str(output_stem.with_suffix(".csv")),
        "output_jsonl": str(output_stem.with_suffix(".jsonl")),
        "output_parquet": str(output_stem.with_suffix(".parquet")),
        "ligand_class_counts": dict(sorted(counts.items())),
        "ligand_class_reason_counts": dict(sorted(reason_counts.items())),
    }
