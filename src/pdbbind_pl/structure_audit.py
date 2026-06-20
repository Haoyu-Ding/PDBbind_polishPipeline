"""Protein and ligand structural audit logic."""

from __future__ import annotations

from collections import defaultdict
import hashlib
from math import dist
from pathlib import Path
import tarfile
from typing import Any

from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import rdMolDescriptors

from pdbbind_pl.utils_io import load_simple_yaml, read_parquet_records, write_csv, write_json, write_jsonl, write_parquet

AA3_TO_AA1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "MSE": "M",
    "SEC": "U",
    "PYL": "O",
}

METAL_CODES = {
    "LI", "NA", "K", "RB", "CS", "MG", "CA", "SR", "BA", "MN", "FE", "CO", "NI", "CU", "ZN",
    "CD", "HG", "AL", "GA", "IN", "TL", "PB", "AG", "AU", "PT", "PD", "CR", "V", "MO", "W", "YB",
}


def run_structure_audit(project_root: Path, config_path: Path) -> dict[str, Any]:
    """Audit protein and ligand structural fields and emit updated manifests."""

    RDLogger.DisableLog("rdApp.warning")

    paths_cfg = load_simple_yaml(config_path)
    workspace_cfg = paths_cfg["workspace"]
    dataset_cfg = paths_cfg["dataset"]
    filters_cfg = load_simple_yaml(project_root / "config" / "filters.yaml")
    active_profile_name = filters_cfg["profiles"]["active"]
    active_profile = filters_cfg[active_profile_name]

    rows = read_parquet_records(Path(workspace_cfg["interim_dir"]) / "master_manifest.parquet")
    audited_rows: list[dict[str, Any]] = []

    with tarfile.open(Path(dataset_cfg["structure_archive"])) as archive:
        for row in rows:
            audited_row = dict(row)
            ligand_result = audit_ligand_members(
                archive=archive,
                sdf_member=str(row.get("ligand_sdf_member") or ""),
                mol2_member=str(row.get("ligand_mol2_member") or ""),
            )
            audited_row.update(ligand_result)

            protein_result = audit_protein_member(
                archive=archive,
                member_name=str(row.get("protein_pdb_member") or ""),
                ligand_mol=ligand_result.get("_ligand_mol"),
                distance_cutoff=float(active_profile["contact_chain_rules"]["heavy_atom_distance_cutoff_angstrom"]),
            )
            audited_row.update(protein_result)
            audited_row.pop("_ligand_mol", None)

            fail_reasons = evaluate_hard_filter_failures(audited_row, active_profile)
            audited_row["hard_filter_fail_reasons"] = fail_reasons
            audited_row["hard_filter_pass"] = len(fail_reasons) == 0
            audited_rows.append(audited_row)

    output_stem = Path(workspace_cfg["interim_dir"]) / "master_manifest_audited"
    field_order = list(audited_rows[0].keys()) if audited_rows else []
    write_csv(output_stem.with_suffix(".csv"), audited_rows, field_order)
    write_jsonl(output_stem.with_suffix(".jsonl"), audited_rows)
    write_parquet(output_stem.with_suffix(".parquet"), audited_rows, field_order)
    summary = build_audit_summary(audited_rows, output_stem, active_profile_name)
    write_json(Path(workspace_cfg["reports_dir"]) / "master_manifest_audited_summary.json", summary)
    return summary


def audit_ligand_members(archive: tarfile.TarFile, sdf_member: str, mol2_member: str) -> dict[str, Any]:
    """Audit ligand parseability and compute core descriptors."""

    mol = load_ligand_mol(archive, sdf_member, mol2_member)
    if mol is None:
        return {
            "ligand_parse_status": "failed_all",
            "ligand_heavy_atom_count": None,
            "ligand_component_count": None,
            "ligand_component_rule": None,
            "ligand_primary_code": None,
            "ligand_formula": None,
            "ligand_mol_wt": None,
            "_ligand_mol": None,
        }

    primary_code = extract_primary_code_from_mol2(archive, mol2_member)
    component_count = len(Chem.GetMolFrags(mol, asMols=False))
    heavy_atom_count = int(mol.GetNumHeavyAtoms())
    return {
        "ligand_parse_status": "ok",
        "ligand_heavy_atom_count": heavy_atom_count,
        "ligand_component_count": component_count,
        "ligand_component_rule": "single_component" if component_count == 1 else "multi_component",
        "ligand_primary_code": primary_code,
        "ligand_formula": rdMolDescriptors.CalcMolFormula(mol),
        "ligand_mol_wt": float(rdMolDescriptors.CalcExactMolWt(mol)),
        "_ligand_mol": mol,
    }


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


def extract_primary_code_from_mol2(archive: tarfile.TarFile, mol2_member: str) -> str | None:
    """Extract a stable primary code from a MOL2 member when available."""

    if not mol2_member:
        return None
    handle = archive.extractfile(mol2_member)
    if handle is None:
        return None
    lines = [line.decode("utf-8", "ignore").rstrip("\n") for line in handle]
    if "@<TRIPOS>ATOM" not in lines or "@<TRIPOS>BOND" not in lines:
        return None
    atom_start = lines.index("@<TRIPOS>ATOM") + 1
    bond_start = lines.index("@<TRIPOS>BOND")
    codes: list[str] = []
    for line in lines[atom_start:bond_start]:
        if not line.strip():
            continue
        parts = line.split()
        if len(parts) >= 8:
            codes.append(parts[7])
    unique_codes = []
    seen = set()
    for code in codes:
        if code not in seen:
            seen.add(code)
            unique_codes.append(code)
    return unique_codes[0] if len(unique_codes) == 1 else None


def audit_protein_member(
    archive: tarfile.TarFile,
    member_name: str,
    ligand_mol: Chem.Mol | None,
    distance_cutoff: float,
) -> dict[str, Any]:
    """Audit chain statistics and determine the primary contact chain."""

    if not member_name:
        return empty_protein_result("missing_member")

    handle = archive.extractfile(member_name)
    if handle is None:
        return empty_protein_result("missing_member")

    chain_ids: list[str] = []
    seen_chain_ids: set[str] = set()
    sequence_by_chain: dict[str, list[str]] = defaultdict(list)
    seen_residues_by_chain: dict[str, set[tuple[str, str, str]]] = defaultdict(set)
    atom_count = 0
    chain_atom_coords: dict[str, list[tuple[float, float, float]]] = defaultdict(list)
    water_residues: set[tuple[str, str, str, str]] = set()
    metal_residues: set[tuple[str, str, str, str]] = set()
    other_hetero_residues: set[tuple[str, str, str, str]] = set()

    try:
        for raw_line in handle:
            line = raw_line.decode("utf-8", "ignore")
            record_name = line[:6]
            chain_id = (line[21].strip() if len(line) > 21 else "") or "_"
            residue_name = line[17:20].strip().upper()
            residue_seq = line[22:26].strip()
            insertion_code = line[26].strip() if len(line) > 26 else ""
            residue_key = (residue_name, residue_seq, insertion_code)
            hetero_key = (residue_name, chain_id, residue_seq, insertion_code)

            if record_name == "ATOM  ":
                atom_count += 1
                if chain_id not in seen_chain_ids:
                    seen_chain_ids.add(chain_id)
                    chain_ids.append(chain_id)
                if residue_key not in seen_residues_by_chain[chain_id]:
                    seen_residues_by_chain[chain_id].add(residue_key)
                    sequence_by_chain[chain_id].append(AA3_TO_AA1.get(residue_name, "X"))
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    chain_atom_coords[chain_id].append((x, y, z))
                except ValueError:
                    continue
            elif record_name == "HETATM":
                if residue_name == "HOH":
                    water_residues.add(hetero_key)
                elif residue_name in METAL_CODES:
                    metal_residues.add(hetero_key)
                else:
                    other_hetero_residues.add(hetero_key)
    except Exception:
        return empty_protein_result("parse_error")

    sequence_strings = {chain_id: "".join(sequence_by_chain[chain_id]) for chain_id in chain_ids}
    chain_lengths = {chain_id: len(sequence_strings[chain_id]) for chain_id in chain_ids}

    contact_counts = compute_contact_counts(chain_atom_coords, ligand_mol, distance_cutoff)
    primary_contact_chain_id, primary_contact_fraction = choose_primary_contact_chain(contact_counts)
    contact_chain_ids = [chain_id for chain_id, count in contact_counts.items() if count > 0]
    primary_sequence = sequence_strings.get(primary_contact_chain_id) if primary_contact_chain_id else None
    sequence_sha1 = hashlib.sha1(primary_sequence.encode("ascii")).hexdigest() if primary_sequence else None

    return {
        "protein_parse_status": "ok",
        "protein_chain_count": len(chain_ids),
        "protein_chain_ids": ",".join(chain_ids) if chain_ids else None,
        "protein_selected_chain_id": primary_contact_chain_id,
        "protein_sequence": primary_sequence,
        "protein_sequence_sha1": sequence_sha1,
        "protein_atom_count": atom_count,
        "protein_water_count": len(water_residues),
        "protein_metal_count": len(metal_residues),
        "protein_other_hetero_count": len(other_hetero_residues),
        "total_protein_chain_count": len(chain_ids),
        "contact_chain_ids": ",".join(contact_chain_ids) if contact_chain_ids else None,
        "contact_chain_count": len(contact_chain_ids),
        "primary_contact_chain_id": primary_contact_chain_id,
        "primary_contact_chain_length": chain_lengths.get(primary_contact_chain_id) if primary_contact_chain_id else None,
        "primary_contact_fraction": primary_contact_fraction,
    }


def empty_protein_result(parse_status: str) -> dict[str, Any]:
    """Return a standard empty protein result payload."""

    return {
        "protein_parse_status": parse_status,
        "protein_chain_count": None,
        "protein_chain_ids": None,
        "protein_selected_chain_id": None,
        "protein_sequence": None,
        "protein_sequence_sha1": None,
        "protein_atom_count": None,
        "protein_water_count": None,
        "protein_metal_count": None,
        "protein_other_hetero_count": None,
        "total_protein_chain_count": None,
        "contact_chain_ids": None,
        "contact_chain_count": None,
        "primary_contact_chain_id": None,
        "primary_contact_chain_length": None,
        "primary_contact_fraction": None,
    }


def compute_contact_counts(
    chain_atom_coords: dict[str, list[tuple[float, float, float]]],
    ligand_mol: Chem.Mol | None,
    distance_cutoff: float,
) -> dict[str, int]:
    """Compute heavy-atom contact counts between ligand atoms and each chain."""

    if ligand_mol is None or ligand_mol.GetNumConformers() == 0:
        return {chain_id: 0 for chain_id in chain_atom_coords}

    conformer = ligand_mol.GetConformer()
    ligand_coords: list[tuple[float, float, float]] = []
    for atom in ligand_mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        pos = conformer.GetAtomPosition(atom.GetIdx())
        ligand_coords.append((pos.x, pos.y, pos.z))

    contact_counts: dict[str, int] = {}
    for chain_id, coords in chain_atom_coords.items():
        count = 0
        for protein_coord in coords:
            for ligand_coord in ligand_coords:
                if dist(protein_coord, ligand_coord) <= distance_cutoff:
                    count += 1
                    break
        contact_counts[chain_id] = count
    return contact_counts


def choose_primary_contact_chain(contact_counts: dict[str, int]) -> tuple[str | None, float | None]:
    """Choose the primary contact chain and compute its contact fraction."""

    if not contact_counts:
        return None, None
    total_contacts = sum(contact_counts.values())
    if total_contacts <= 0:
        return None, None
    primary_chain = max(sorted(contact_counts), key=lambda chain_id: contact_counts[chain_id])
    primary_fraction = contact_counts[primary_chain] / total_contacts
    return primary_chain, primary_fraction


def evaluate_hard_filter_failures(row: dict[str, Any], profile: dict[str, Any]) -> list[str]:
    """Evaluate strict-v2 hard filter rules against an audited row."""

    fail_reasons: list[str] = []
    hard_exclusions = profile["hard_exclusions"]
    contact_rules = profile["contact_chain_rules"]
    ligand_size_rules = profile["ligand_size_rules"]
    protein_size_rules = profile["protein_size_rules"]

    if hard_exclusions["exclude_covalent_complex"] and row.get("flag_covalent_complex"):
        fail_reasons.append("covalent_complex")
    if hard_exclusions["exclude_incomplete_ligand"] and row.get("flag_incomplete_ligand"):
        fail_reasons.append("incomplete_ligand")
    if hard_exclusions["exclude_peptide_like_n_mer"] and row.get("flag_peptide_like_mer_annotation"):
        fail_reasons.append("peptide_like_n_mer")
    if hard_exclusions["exclude_isomer_annotations"] and row.get("flag_isomer_annotation"):
        fail_reasons.append("isomer_annotation")
    if hard_exclusions["exclude_unparseable_ligand"] and row.get("ligand_parse_status") != "ok":
        fail_reasons.append("unparseable_ligand")
    if hard_exclusions["exclude_multi_component_ligand"] and row.get("ligand_component_count") != 1:
        fail_reasons.append("multi_component_ligand")
    if hard_exclusions["exclude_no_primary_contact_chain"] and not row.get("primary_contact_chain_id"):
        fail_reasons.append("no_primary_contact_chain")

    primary_fraction = row.get("primary_contact_fraction")
    if (
        hard_exclusions["exclude_interface_dominated_binding"]
        and primary_fraction is not None
        and float(primary_fraction) < float(contact_rules["minimum_primary_contact_fraction"])
    ):
        fail_reasons.append("interface_dominated_binding")

    ligand_mol_wt = row.get("ligand_mol_wt")
    if (
        hard_exclusions["exclude_ligand_below_min_mw"]
        and ligand_mol_wt is not None
        and float(ligand_mol_wt) <= float(ligand_size_rules["min_mol_wt_gt"])
    ):
        fail_reasons.append("ligand_below_min_mw")

    primary_chain_length = row.get("primary_contact_chain_length")
    if (
        hard_exclusions["exclude_primary_chain_above_max_length"]
        and primary_chain_length is not None
        and int(primary_chain_length) > int(protein_size_rules["max_primary_chain_length_lte"])
    ):
        fail_reasons.append("primary_chain_above_max_length")

    if row.get("protein_parse_status") != "ok":
        fail_reasons.append("protein_parse_failure")
    return fail_reasons


def build_audit_summary(rows: list[dict[str, Any]], output_stem: Path, active_profile_name: str) -> dict[str, Any]:
    """Build a summary report for the audited manifest."""

    fail_reason_counts: dict[str, int] = defaultdict(int)
    for row in rows:
        for reason in row.get("hard_filter_fail_reasons") or []:
            fail_reason_counts[reason] += 1

    return {
        "active_profile": active_profile_name,
        "row_count": len(rows),
        "output_csv": str(output_stem.with_suffix(".csv")),
        "output_jsonl": str(output_stem.with_suffix(".jsonl")),
        "output_parquet": str(output_stem.with_suffix(".parquet")),
        "hard_filter_pass_count": sum(bool(row.get("hard_filter_pass")) for row in rows),
        "protein_parse_status_counts": count_by(rows, "protein_parse_status"),
        "ligand_parse_status_counts": count_by(rows, "ligand_parse_status"),
        "protein_chain_count_distribution": count_by(rows, "protein_chain_count"),
        "contact_chain_count_distribution": count_by(rows, "contact_chain_count"),
        "ligand_component_count_distribution": count_by(rows, "ligand_component_count"),
        "hard_filter_fail_reason_counts": dict(sorted(fail_reason_counts.items())),
    }


def count_by(rows: list[dict[str, Any]], field: str) -> dict[str, int]:
    """Count values for a field in a row collection."""

    counts: dict[str, int] = defaultdict(int)
    for row in rows:
        key = row.get(field)
        counts[str(key)] += 1
    return dict(sorted(counts.items(), key=lambda item: item[0]))
