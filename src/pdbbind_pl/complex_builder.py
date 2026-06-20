"""Final complex export logic."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import shutil
import tarfile
from typing import Any

from rdkit import Chem
from rdkit import RDLogger

from pdbbind_pl.utils_io import (
    ensure_directory,
    load_simple_yaml,
    read_parquet_records,
    write_csv,
    write_json,
    write_jsonl,
    write_parquet,
)


@dataclass
class ExportArtifacts:
    """File paths for one exported complex."""

    complex_pdb_path: Path
    protein_chain_pdb_path: Path
    ligand_structure_path: Path


def run_complex_export(project_root: Path, config_path: Path) -> dict[str, Any]:
    """Export cleaned complex PDB files for final dataset representatives."""

    RDLogger.DisableLog("rdApp.warning")

    paths_cfg = load_simple_yaml(config_path)
    workspace_cfg = paths_cfg["workspace"]
    dataset_cfg = paths_cfg["dataset"]

    rows = read_parquet_records(Path(workspace_cfg["interim_dir"]) / "master_manifest_deduped.parquet")
    output_rows = [dict(row) for row in rows]

    structure_archive = Path(dataset_cfg["structure_archive"])
    final_root = Path(workspace_cfg["final_dir"])
    if final_root.exists():
        shutil.rmtree(final_root)
    ensure_directory(final_root)

    exported_count = 0
    skipped_count = 0
    error_count = 0

    with tarfile.open(structure_archive) as archive:
        for row in output_rows:
            bucket = row.get("final_dataset_bucket")
            if bucket is None or (isinstance(bucket, float) and str(bucket) == "nan"):
                row["final_export_status"] = "not_selected"
                continue

            try:
                artifacts = export_single_complex(
                    archive=archive,
                    row=row,
                    final_root=final_root,
                )
                row["final_complex_pdb_path"] = str(artifacts.complex_pdb_path)
                row["final_protein_chain_pdb_path"] = str(artifacts.protein_chain_pdb_path)
                row["final_ligand_structure_path"] = str(artifacts.ligand_structure_path)
                row["final_export_status"] = "exported"
                exported_count += 1
            except SkipExport as exc:
                row["final_export_status"] = exc.reason
                skipped_count += 1
            except Exception as exc:  # pragma: no cover - batch robustness path
                row["final_export_status"] = f"export_error:{type(exc).__name__}"
                error_count += 1

    output_stem = Path(workspace_cfg["interim_dir"]) / "master_manifest_exported"
    field_order = list(output_rows[0].keys()) if output_rows else []
    write_csv(output_stem.with_suffix(".csv"), output_rows, field_order)
    write_jsonl(output_stem.with_suffix(".jsonl"), output_rows)
    write_parquet(output_stem.with_suffix(".parquet"), output_rows, field_order)

    summary = build_export_summary(
        rows=output_rows,
        output_stem=output_stem,
        exported_count=exported_count,
        skipped_count=skipped_count,
        error_count=error_count,
    )
    write_json(Path(workspace_cfg["reports_dir"]) / "master_manifest_exported_summary.json", summary)
    return summary


class SkipExport(Exception):
    """Raised when a row is intentionally skipped during export."""

    def __init__(self, reason: str) -> None:
        super().__init__(reason)
        self.reason = reason


def export_single_complex(archive: tarfile.TarFile, row: dict[str, Any], final_root: Path) -> ExportArtifacts:
    """Export one cleaned complex and its components."""

    bucket = str(row["final_dataset_bucket"])
    pdb_id = str(row["pdb_id"])
    chain_id = row.get("primary_contact_chain_id") or row.get("protein_selected_chain_id")
    if not chain_id:
        raise SkipExport("missing_selected_chain")

    protein_lines = extract_protein_chain_lines(archive, str(row.get("protein_pdb_member") or ""), str(chain_id))
    if not protein_lines:
        raise SkipExport("empty_protein_chain")

    ligand_mol = load_ligand_mol(
        archive=archive,
        sdf_member=str(row.get("ligand_sdf_member") or ""),
        mol2_member=str(row.get("ligand_mol2_member") or ""),
    )
    if ligand_mol is None:
        raise SkipExport("ligand_parse_failure")

    residue_name = infer_ligand_residue_name(row)
    ligand_lines = build_ligand_pdb_lines(ligand_mol, residue_name)
    if not ligand_lines:
        raise SkipExport("empty_ligand")

    target_dir = final_root / bucket / pdb_id
    ensure_directory(target_dir)
    protein_path = target_dir / f"{pdb_id}_protein_chain.pdb"
    ligand_path = target_dir / f"{pdb_id}_ligand.pdb"
    complex_path = target_dir / f"{pdb_id}_complex.pdb"

    protein_path.write_text("".join(protein_lines + ["END\n"]), encoding="utf-8")
    ligand_path.write_text("".join(ligand_lines + ["END\n"]), encoding="utf-8")
    complex_path.write_text("".join(protein_lines + ["TER\n"] + ligand_lines + ["END\n"]), encoding="utf-8")

    return ExportArtifacts(
        complex_pdb_path=complex_path,
        protein_chain_pdb_path=protein_path,
        ligand_structure_path=ligand_path,
    )


def extract_protein_chain_lines(archive: tarfile.TarFile, member_name: str, chain_id: str) -> list[str]:
    """Extract ATOM lines for a single protein chain."""

    handle = archive.extractfile(member_name)
    if handle is None:
        return []

    lines: list[str] = []
    target_chain = " " if chain_id == "_" else chain_id
    for raw_line in handle:
        line = raw_line.decode("utf-8", "ignore")
        if not line.startswith("ATOM  "):
            continue
        if len(line) <= 21:
            continue
        if line[21] != target_chain:
            continue
        lines.append(ensure_pdb_line(line))
    return lines


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


def infer_ligand_residue_name(row: dict[str, Any]) -> str:
    """Infer a residue name for exported ligand atoms."""

    primary_code = row.get("ligand_primary_code")
    if isinstance(primary_code, str) and primary_code:
        return primary_code[:3].upper()
    ligand_name = row.get("ligand_name_raw")
    if isinstance(ligand_name, str) and ligand_name:
        trimmed = "".join(ch for ch in ligand_name.upper() if ch.isalnum())
        return (trimmed[:3] or "LIG").ljust(3)
    return "LIG"


def build_ligand_pdb_lines(mol: Chem.Mol, residue_name: str) -> list[str]:
    """Build ligand HETATM lines from an RDKit molecule."""

    if mol.GetNumConformers() == 0:
        return []

    conformer = mol.GetConformer()
    lines: list[str] = []
    for atom_index, atom in enumerate(mol.GetAtoms(), start=1):
        position = conformer.GetAtomPosition(atom_index - 1)
        atom_name = format_atom_name(atom)
        element = atom.GetSymbol().upper()
        line = (
            f"HETATM{atom_index:5d} {atom_name:<4s} {residue_name:>3s} X{1:4d}    "
            f"{position.x:8.3f}{position.y:8.3f}{position.z:8.3f}"
            f"{1.00:6.2f}{0.00:6.2f}          {element:>2s}\n"
        )
        lines.append(line)
    return lines


def format_atom_name(atom: Chem.Atom) -> str:
    """Create a compact atom name for PDB export."""

    symbol = atom.GetSymbol().upper()
    serial = atom.GetIdx() + 1
    return f"{symbol}{serial}"[:4]


def ensure_pdb_line(line: str) -> str:
    """Ensure a PDB record line ends with a newline."""

    return line if line.endswith("\n") else f"{line}\n"


def build_export_summary(
    rows: list[dict[str, Any]],
    output_stem: Path,
    exported_count: int,
    skipped_count: int,
    error_count: int,
) -> dict[str, Any]:
    """Build a summary of the export stage."""

    final_status_counts: dict[str, int] = {}
    for row in rows:
        status = str(row.get("final_export_status"))
        final_status_counts[status] = final_status_counts.get(status, 0) + 1

    return {
        "row_count": len(rows),
        "output_csv": str(output_stem.with_suffix(".csv")),
        "output_jsonl": str(output_stem.with_suffix(".jsonl")),
        "output_parquet": str(output_stem.with_suffix(".parquet")),
        "exported_count": exported_count,
        "skipped_count": skipped_count,
        "error_count": error_count,
        "final_export_status_counts": dict(sorted(final_status_counts.items())),
    }


def analyze_exported_complex(path: Path) -> dict[str, Any]:
    """Analyze a simple exported complex PDB file."""

    atom_chain_ids: set[str] = set()
    atom_count = 0
    hetatm_count = 0
    has_ter = False
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        if raw_line.startswith("ATOM  "):
            atom_count += 1
            chain_id = (raw_line[21].strip() if len(raw_line) > 21 else "") or "_"
            atom_chain_ids.add(chain_id)
        elif raw_line.startswith("HETATM"):
            hetatm_count += 1
        elif raw_line == "TER":
            has_ter = True
    return {
        "atom_chain_count": len(atom_chain_ids),
        "atom_count": atom_count,
        "hetatm_count": hetatm_count,
        "has_ter": has_ter,
    }
