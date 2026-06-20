"""Shared IO helpers for archive access and table writing."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml


def ensure_directory(path: Path) -> None:
    """Create a directory if it does not already exist."""

    path.mkdir(parents=True, exist_ok=True)


def ensure_parent_directory(path: Path) -> None:
    """Create the parent directory for a file path."""

    ensure_directory(path.parent)


def load_simple_yaml(path: Path) -> Any:
    """Load YAML using PyYAML."""

    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def write_json(path: Path, payload: Any) -> None:
    """Write indented JSON to disk."""

    ensure_parent_directory(path)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def write_jsonl(path: Path, records: list[dict[str, Any]]) -> None:
    """Write JSON Lines records."""

    ensure_parent_directory(path)
    with path.open("w", encoding="utf-8") as handle:
        for record in records:
            normalized = {key: normalize_json_value(value) for key, value in record.items()}
            handle.write(json.dumps(normalized, ensure_ascii=True, sort_keys=False))
            handle.write("\n")


def write_csv(path: Path, records: list[dict[str, Any]], fieldnames: list[str]) -> None:
    """Write a list of dictionaries as CSV."""

    ensure_parent_directory(path)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for record in records:
            row = {field: normalize_tabular_value(record.get(field)) for field in fieldnames}
            writer.writerow(row)


def write_parquet(path: Path, records: list[dict[str, Any]], fieldnames: list[str] | None = None) -> None:
    """Write Parquet using pandas and pyarrow."""

    ensure_parent_directory(path)
    if fieldnames is None:
        fieldnames = infer_fieldnames(records)
    frame = pd.DataFrame(records)
    if fieldnames:
        for field in fieldnames:
            if field not in frame.columns:
                frame[field] = None
        frame = frame[fieldnames]
    frame.to_parquet(path, index=False)


def read_parquet_records(path: Path) -> list[dict[str, Any]]:
    """Read Parquet into a list of dictionaries."""

    frame = pd.read_parquet(path)
    frame = frame.where(pd.notna(frame), None)
    return frame.to_dict(orient="records")


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    """Read a JSON Lines file into memory."""

    records: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                records.append(json.loads(line))
    return records


def normalize_tabular_value(value: Any) -> Any:
    """Convert nested values to stable flat text for CSV output."""

    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, np.ndarray):
        return json.dumps(value.tolist(), sort_keys=True)
    if isinstance(value, (list, dict)):
        return json.dumps(value, sort_keys=True)
    return value


def infer_fieldnames(records: list[dict[str, Any]]) -> list[str]:
    """Infer a stable field order from a list of dictionaries."""

    fieldnames: list[str] = []
    seen: set[str] = set()
    for record in records:
        for key in record:
            if key not in seen:
                seen.add(key)
                fieldnames.append(key)
    return fieldnames


def normalize_json_value(value: Any) -> Any:
    """Convert array-like values into JSON-serializable Python values."""

    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, list):
        return [normalize_json_value(item) for item in value]
    if isinstance(value, dict):
        return {key: normalize_json_value(item) for key, item in value.items()}
    return value
