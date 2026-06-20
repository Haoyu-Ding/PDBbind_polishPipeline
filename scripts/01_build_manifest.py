#!/usr/bin/env python3
"""Stage 01: build the master manifest from index and archive paths."""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import sys
import tarfile

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from pdbbind_pl.index_parser import load_index_records
from pdbbind_pl.utils_io import load_simple_yaml, write_csv, write_json, write_jsonl, write_parquet

MANIFEST_FIELD_ORDER = [
    "pdb_id",
    "release_year",
    "structure_archive_member_dir",
    "protein_pdb_member",
    "ligand_sdf_member",
    "ligand_mol2_member",
    "source_index_member",
    "resolution_raw",
    "resolution_value",
    "structure_method_class",
    "binding_data_raw",
    "reference_raw",
    "ligand_name_raw",
    "index_comment_raw",
    "flag_covalent_complex",
    "flag_incomplete_ligand",
    "flag_isomer_annotation",
    "flag_redundant_annotation",
    "flag_peptide_like_mer_annotation",
    "peptide_like_mer_size",
    "ligand_parse_status",
    "ligand_formula",
    "ligand_heavy_atom_count",
    "ligand_mol_wt",
    "ligand_rotatable_bonds",
    "ligand_component_count",
    "ligand_component_rule",
    "ligand_class",
    "ligand_class_reason",
    "ligand_primary_code",
    "protein_parse_status",
    "protein_chain_count",
    "protein_chain_ids",
    "protein_selected_chain_id",
    "protein_sequence",
    "protein_sequence_sha1",
    "protein_atom_count",
    "protein_water_count",
    "protein_metal_count",
    "protein_other_hetero_count",
    "hard_filter_pass",
    "hard_filter_fail_reasons",
    "nonsugar_flexibility_bucket",
    "ligand_dedup_cluster_id",
    "ligand_dedup_is_representative",
    "protein_dedup_cluster_id",
    "protein_dedup_is_representative",
    "final_dataset_bucket",
    "final_export_status",
    "final_complex_pdb_path",
    "final_protein_chain_pdb_path",
    "final_ligand_structure_path",
    "validation_status",
    "validation_notes",
]

STRUCTURE_FILE_RE = re.compile(
    r"^P-L/(?P<year_bucket>[^/]+)/(?P<pdb_id>[0-9a-z]{4})/(?P<filename>[^/]+)$",
    re.IGNORECASE,
)


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

    index_records = load_index_records(
        index_archive=Path(dataset_cfg["index_archive"]),
        index_member=dataset_cfg["index_member"],
    )
    structure_index = collect_structure_members(Path(dataset_cfg["structure_archive"]))

    manifest_rows = [build_manifest_row(record, structure_index.get(record["pdb_id"])) for record in index_records]

    interim_dir = Path(workspace_cfg["interim_dir"])
    reports_dir = Path(workspace_cfg["reports_dir"])
    output_stem = interim_dir / "master_manifest"

    write_csv(output_stem.with_suffix(".csv"), manifest_rows, MANIFEST_FIELD_ORDER)
    write_jsonl(output_stem.with_suffix(".jsonl"), manifest_rows)
    write_parquet(output_stem.with_suffix(".parquet"), manifest_rows, MANIFEST_FIELD_ORDER)
    write_json(
        reports_dir / "master_manifest_summary.json",
        build_summary(
            manifest_rows=manifest_rows,
            structure_archive=Path(dataset_cfg["structure_archive"]),
            output_stem=output_stem,
        ),
    )

    print(f"Wrote {len(manifest_rows)} manifest rows")
    print(f"CSV: {output_stem.with_suffix('.csv')}")
    print(f"JSONL: {output_stem.with_suffix('.jsonl')}")
    print(f"Parquet: {output_stem.with_suffix('.parquet')}")


def collect_structure_members(structure_archive: Path) -> dict[str, dict[str, str]]:
    """Index structure archive members by PDB code."""

    structure_index: dict[str, dict[str, str]] = {}
    with tarfile.open(structure_archive) as archive:
        for member in archive:
            if not member.isfile():
                continue
            match = STRUCTURE_FILE_RE.match(member.name)
            if match is None:
                continue

            pdb_id = match.group("pdb_id").lower()
            filename = match.group("filename")
            entry = structure_index.setdefault(
                pdb_id,
                {
                    "structure_archive_member_dir": str(Path(member.name).parent),
                    "protein_pdb_member": "",
                    "ligand_sdf_member": "",
                    "ligand_mol2_member": "",
                },
            )

            lower_filename = filename.lower()
            if lower_filename.endswith("_protein.pdb"):
                entry["protein_pdb_member"] = member.name
            elif lower_filename.endswith("_ligand.sdf"):
                entry["ligand_sdf_member"] = member.name
            elif lower_filename.endswith("_ligand.mol2"):
                entry["ligand_mol2_member"] = member.name
    return structure_index


def build_manifest_row(index_record: dict[str, object], structure_record: dict[str, str] | None) -> dict[str, object]:
    """Combine index metadata with structure archive paths."""

    row = {
        "pdb_id": index_record["pdb_id"],
        "release_year": index_record["release_year"],
        "structure_archive_member_dir": "",
        "protein_pdb_member": "",
        "ligand_sdf_member": "",
        "ligand_mol2_member": "",
        "source_index_member": index_record["source_index_member"],
        "resolution_raw": index_record["resolution_raw"],
        "resolution_value": index_record["resolution_value"],
        "structure_method_class": index_record["structure_method_class"],
        "binding_data_raw": index_record["binding_data_raw"],
        "reference_raw": index_record["reference_raw"],
        "ligand_name_raw": index_record["ligand_name_raw"],
        "index_comment_raw": index_record["index_comment_raw"],
        "flag_covalent_complex": index_record["flag_covalent_complex"],
        "flag_incomplete_ligand": index_record["flag_incomplete_ligand"],
        "flag_isomer_annotation": index_record["flag_isomer_annotation"],
        "flag_redundant_annotation": index_record["flag_redundant_annotation"],
        "flag_peptide_like_mer_annotation": index_record["flag_peptide_like_mer_annotation"],
        "peptide_like_mer_size": index_record["peptide_like_mer_size"],
        "ligand_parse_status": "not_run",
        "ligand_formula": None,
        "ligand_heavy_atom_count": None,
        "ligand_mol_wt": None,
        "ligand_rotatable_bonds": None,
        "ligand_component_count": None,
        "ligand_component_rule": None,
        "ligand_class": None,
        "ligand_class_reason": None,
        "ligand_primary_code": None,
        "protein_parse_status": "not_run",
        "protein_chain_count": None,
        "protein_chain_ids": None,
        "protein_selected_chain_id": None,
        "protein_sequence": None,
        "protein_sequence_sha1": None,
        "protein_atom_count": None,
        "protein_water_count": None,
        "protein_metal_count": None,
        "protein_other_hetero_count": None,
        "hard_filter_pass": None,
        "hard_filter_fail_reasons": None,
        "nonsugar_flexibility_bucket": None,
        "ligand_dedup_cluster_id": None,
        "ligand_dedup_is_representative": None,
        "protein_dedup_cluster_id": None,
        "protein_dedup_is_representative": None,
        "final_dataset_bucket": None,
        "final_export_status": "pending",
        "final_complex_pdb_path": None,
        "final_protein_chain_pdb_path": None,
        "final_ligand_structure_path": None,
        "validation_status": None,
        "validation_notes": None,
    }

    if structure_record is not None:
        row["structure_archive_member_dir"] = structure_record["structure_archive_member_dir"]
        row["protein_pdb_member"] = structure_record["protein_pdb_member"]
        row["ligand_sdf_member"] = structure_record["ligand_sdf_member"]
        row["ligand_mol2_member"] = structure_record["ligand_mol2_member"]

    return row


def build_summary(
    manifest_rows: list[dict[str, object]],
    structure_archive: Path,
    output_stem: Path,
) -> dict[str, object]:
    """Build a compact manifest summary."""

    missing_protein = sum(not bool(row["protein_pdb_member"]) for row in manifest_rows)
    missing_sdf = sum(not bool(row["ligand_sdf_member"]) for row in manifest_rows)
    missing_mol2 = sum(not bool(row["ligand_mol2_member"]) for row in manifest_rows)

    return {
        "row_count": len(manifest_rows),
        "structure_archive": str(structure_archive),
        "output_csv": str(output_stem.with_suffix(".csv")),
        "output_jsonl": str(output_stem.with_suffix(".jsonl")),
        "output_parquet": str(output_stem.with_suffix(".parquet")),
        "missing_member_counts": {
            "protein_pdb_member": missing_protein,
            "ligand_sdf_member": missing_sdf,
            "ligand_mol2_member": missing_mol2,
        },
    }


if __name__ == "__main__":
    main()
