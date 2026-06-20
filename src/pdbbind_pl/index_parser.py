"""Utilities for parsing PDBbind index files."""

from __future__ import annotations

import re
import tarfile
from pathlib import Path
from typing import Any

INDEX_FIELD_ORDER = [
    "pdb_id",
    "release_year",
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
    "source_index_member",
]

INDEX_LINE_RE = re.compile(
    r"^(?P<pdb_id>\S+)\s+"
    r"(?P<resolution_raw>\S+)\s+"
    r"(?P<release_year>\d{4})\s+"
    r"(?P<binding_data_raw>\S+)\s*//\s*"
    r"(?P<index_comment_raw>.*)$"
)

MER_RE = re.compile(r"\((?P<mer_size>\d+)-mer\)", re.IGNORECASE)


def load_index_records(index_archive: Path, index_member: str) -> list[dict[str, Any]]:
    """Read and normalize the PDBbind general protein-ligand index."""

    records: list[dict[str, Any]] = []
    with tarfile.open(index_archive) as archive:
        member = archive.extractfile(index_member)
        if member is None:
            raise FileNotFoundError(f"Could not open {index_member} inside {index_archive}")
        for raw_line in member:
            line = raw_line.decode("utf-8", "ignore").rstrip("\n")
            if not line or line.startswith("#"):
                continue
            records.append(parse_index_line(line, index_member))
    return records


def parse_index_line(line: str, source_index_member: str) -> dict[str, Any]:
    """Parse one normalized PDBbind index line."""

    match = INDEX_LINE_RE.match(line.strip())
    if match is None:
        raise ValueError(f"Could not parse index line: {line}")

    row = match.groupdict()
    comment = row["index_comment_raw"].strip()
    reference_raw, ligand_name_raw = split_comment_fields(comment)
    mer_match = MER_RE.search(comment)

    resolution_raw = row["resolution_raw"]
    resolution_value = parse_resolution_value(resolution_raw)

    return {
        "pdb_id": row["pdb_id"].lower(),
        "release_year": int(row["release_year"]),
        "resolution_raw": resolution_raw,
        "resolution_value": resolution_value,
        "structure_method_class": classify_structure_method(resolution_raw),
        "binding_data_raw": row["binding_data_raw"],
        "reference_raw": reference_raw,
        "ligand_name_raw": ligand_name_raw,
        "index_comment_raw": comment,
        "flag_covalent_complex": "covalent complex" in comment.lower(),
        "flag_incomplete_ligand": "incomplete ligand" in comment.lower(),
        "flag_isomer_annotation": " isomer" in f" {comment.lower()}",
        "flag_redundant_annotation": "redundant to" in comment.lower(),
        "flag_peptide_like_mer_annotation": mer_match is not None,
        "peptide_like_mer_size": int(mer_match.group("mer_size")) if mer_match else None,
        "source_index_member": source_index_member,
    }


def split_comment_fields(comment: str) -> tuple[str | None, str | None]:
    """Split the raw comment into reference and the first ligand annotation."""

    parts = comment.split(maxsplit=1)
    reference_raw = parts[0] if parts else None
    ligand_name_match = re.search(r"\(([^()]*)\)", comment)
    ligand_name_raw = ligand_name_match.group(1).strip() if ligand_name_match else None
    return reference_raw, ligand_name_raw


def parse_resolution_value(resolution_raw: str) -> float | None:
    """Return a numeric resolution when available."""

    try:
        return float(resolution_raw)
    except ValueError:
        return None


def classify_structure_method(resolution_raw: str) -> str:
    """Map resolution tokens to a coarse structure method class."""

    if resolution_raw.upper() == "NMR":
        return "nmr"
    if parse_resolution_value(resolution_raw) is not None:
        return "xray"
    return "other"
