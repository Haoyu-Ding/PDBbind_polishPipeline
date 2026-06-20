"""Protein clustering helpers for subset selection."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
import subprocess
from typing import Any

from pdbbind_pl.utils_io import ensure_directory, load_simple_yaml, read_parquet_records, write_csv, write_json, write_jsonl, write_parquet


def run_protein_clustering(project_root: Path, config_path: Path) -> dict[str, Any]:
    """Cluster validated exported proteins with MMseqs2 at 0.7 identity."""

    paths_cfg = load_simple_yaml(config_path)
    workspace_cfg = paths_cfg["workspace"]

    rows = read_parquet_records(Path(workspace_cfg["interim_dir"]) / "master_manifest_validated.parquet")
    candidate_rows = [dict(row) for row in rows if row.get("validation_status") == "ok"]

    subset_dir = Path(workspace_cfg["root_dir"]) / "data" / "subsets" / "subset50"
    ensure_directory(subset_dir)
    fasta_path = subset_dir / "subset50_candidates.fasta"
    cluster_tsv_path = subset_dir / "subset50_protein_clusters.tsv"
    mmseqs_tmp_dir = subset_dir / "mmseqs_tmp"
    ensure_directory(mmseqs_tmp_dir)

    write_candidate_fasta(candidate_rows, fasta_path)
    run_mmseqs_easy_cluster(fasta_path, subset_dir / "subset50_mmseqs", mmseqs_tmp_dir)
    raw_cluster_tsv = subset_dir / "subset50_mmseqs_cluster.tsv"
    raw_cluster_map = load_cluster_assignments(raw_cluster_tsv)
    normalized_clusters = normalize_clusters(raw_cluster_map)
    write_cluster_tsv(cluster_tsv_path, normalized_clusters)

    updated_rows = []
    for row in rows:
        updated_row = dict(row)
        pdb_id = str(updated_row["pdb_id"])
        updated_row["protein_cluster_0p7"] = normalized_clusters.get(pdb_id)
        updated_rows.append(updated_row)

    output_stem = Path(workspace_cfg["interim_dir"]) / "master_manifest_clustered_subset50"
    field_order = list(updated_rows[0].keys()) if updated_rows else []
    write_csv(output_stem.with_suffix(".csv"), updated_rows, field_order)
    write_jsonl(output_stem.with_suffix(".jsonl"), updated_rows)
    write_parquet(output_stem.with_suffix(".parquet"), updated_rows, field_order)

    summary = build_cluster_summary(updated_rows, output_stem, cluster_tsv_path)
    write_json(Path(workspace_cfg["reports_dir"]) / "master_manifest_clustered_subset50_summary.json", summary)
    return summary


def write_candidate_fasta(rows: list[dict[str, Any]], fasta_path: Path) -> None:
    """Write one FASTA record per validated exported entry."""

    lines: list[str] = []
    for row in rows:
        sequence = row.get("protein_sequence")
        if not isinstance(sequence, str) or not sequence:
            continue
        pdb_id = str(row["pdb_id"])
        lines.append(f">{pdb_id}\n")
        lines.append(f"{sequence}\n")
    fasta_path.write_text("".join(lines), encoding="utf-8")


def run_mmseqs_easy_cluster(fasta_path: Path, result_prefix: Path, tmp_dir: Path) -> None:
    """Run MMseqs2 easy-cluster with 0.7 identity."""

    command = [
        "mmseqs",
        "easy-cluster",
        str(fasta_path),
        str(result_prefix),
        str(tmp_dir),
        "--min-seq-id",
        "0.7",
        "-c",
        "0.8",
        "--cov-mode",
        "0",
    ]
    subprocess.run(command, check=True)


def load_cluster_assignments(cluster_tsv_path: Path) -> dict[str, str]:
    """Load raw representative -> member assignments from MMseqs2 output."""

    assignments: dict[str, str] = {}
    for line in cluster_tsv_path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        representative, member = line.split("\t")
        assignments[member] = representative
    return assignments


def normalize_clusters(raw_assignments: dict[str, str]) -> dict[str, str]:
    """Rename raw representative ids to stable cluster IDs."""

    representative_to_cluster: dict[str, str] = {}
    cluster_members: dict[str, list[str]] = defaultdict(list)
    for member, representative in raw_assignments.items():
        cluster_members[representative].append(member)

    for index, representative in enumerate(sorted(cluster_members), start=1):
        representative_to_cluster[representative] = f"protein_cluster_0p7_{index:05d}"

    normalized: dict[str, str] = {}
    for member, representative in raw_assignments.items():
        normalized[member] = representative_to_cluster[representative]
    return normalized


def write_cluster_tsv(path: Path, assignments: dict[str, str]) -> None:
    """Write stable member -> cluster assignments."""

    lines = [f"{member}\t{cluster_id}\n" for member, cluster_id in sorted(assignments.items())]
    path.write_text("".join(lines), encoding="utf-8")


def build_cluster_summary(
    rows: list[dict[str, Any]],
    output_stem: Path,
    cluster_tsv_path: Path,
) -> dict[str, Any]:
    """Build summary statistics for protein clustering."""

    cluster_counts: dict[str, int] = defaultdict(int)
    assigned_count = 0
    for row in rows:
        cluster_id = row.get("protein_cluster_0p7")
        if isinstance(cluster_id, str) and cluster_id:
            assigned_count += 1
            cluster_counts[cluster_id] += 1

    return {
        "row_count": len(rows),
        "assigned_cluster_count": assigned_count,
        "unique_cluster_count": len(cluster_counts),
        "cluster_tsv": str(cluster_tsv_path),
        "output_csv": str(output_stem.with_suffix(".csv")),
        "output_jsonl": str(output_stem.with_suffix(".jsonl")),
        "output_parquet": str(output_stem.with_suffix(".parquet")),
    }
