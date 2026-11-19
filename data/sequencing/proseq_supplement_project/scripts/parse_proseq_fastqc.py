#!/usr/bin/env python3
"""Extract cleaning-relevant metrics from per-file FastQC archives."""

from __future__ import annotations

import argparse
import csv
import re
import zipfile
from pathlib import Path


def read_modules(zip_path: Path):
    with zipfile.ZipFile(zip_path) as archive:
        member = next(name for name in archive.namelist() if name.endswith("/fastqc_data.txt"))
        lines = archive.read(member).decode("utf-8", errors="replace").splitlines()

    modules = {}
    name = status = None
    headers = []
    rows = []
    for line in lines:
        if line.startswith(">>") and line != ">>END_MODULE":
            name, status = line[2:].split("\t", 1)
            headers, rows = [], []
        elif line == ">>END_MODULE":
            modules[name] = {"status": status, "headers": headers, "rows": rows}
            name = status = None
        elif name is not None:
            if line.startswith("#"):
                headers = line[1:].split("\t")
            elif not line.startswith("##"):
                rows.append(line.split("\t"))
    return modules


def module_dict(module):
    return {row[0]: row[1] for row in module.get("rows", []) if len(row) >= 2}


def base_endpoint(text: str) -> int:
    numbers = re.findall(r"\d+", text)
    return int(numbers[-1]) if numbers else 0


def fmt(value):
    if value is None:
        return "NA"
    if isinstance(value, float):
        return f"{value:.4f}"
    return str(value)


def parse_archive(zip_path: Path, metadata):
    modules = read_modules(zip_path)
    basic = module_dict(modules.get("Basic Statistics", {}))

    quality = modules.get("Per base sequence quality", {})
    quality_rows = []
    for row in quality.get("rows", []):
        if len(row) >= 2:
            quality_rows.append((base_endpoint(row[0]), float(row[1])))
    quality_rows.sort()

    n_module = modules.get("Per base N content", {})
    n_rows = []
    for row in n_module.get("rows", []):
        if len(row) >= 2:
            n_rows.append((base_endpoint(row[0]), float(row[1])))
    n_rows.sort()

    length_module = modules.get("Sequence Length Distribution", {})
    length_rows = []
    for row in length_module.get("rows", []):
        if len(row) >= 2 and float(row[1]) > 0:
            length_rows.append((row[0], base_endpoint(row[0]), float(row[1])))

    adapter_module = modules.get("Adapter Content", {})
    adapter_names = adapter_module.get("headers", [])[1:]
    adapter_values = {name: [] for name in adapter_names}
    for row in adapter_module.get("rows", []):
        for index, name in enumerate(adapter_names, start=1):
            if index < len(row):
                adapter_values[name].append(float(row[index]))

    adapter_summary = []
    for name, values in adapter_values.items():
        maximum = max(values, default=0.0)
        final = values[-1] if values else 0.0
        if maximum >= 0.1:
            adapter_summary.append(f"{name}:{maximum:.3f}")

    overrep = modules.get("Overrepresented sequences", {})
    overrep_rows = []
    for row in overrep.get("rows", [])[:10]:
        if len(row) >= 4:
            overrep_rows.append(f"{row[0]}|{row[2]}%|{row[3]}")

    filename = basic.get("Filename", zip_path.name.replace("_fastqc.zip", ".fastq.gz"))
    meta = metadata.get(filename, {})
    result = {
        "method": meta.get("method", "NA"),
        "group_id": meta.get("group_id", "NA"),
        "gse": meta.get("gse", "NA"),
        "gsm": meta.get("gsm", "NA"),
        "run": meta.get("run", "NA"),
        "layout": meta.get("layout", "NA"),
        "mate": meta.get("mate", "NA"),
        "filename": filename,
        "total_sequences": basic.get("Total Sequences", "NA"),
        "basic_length": basic.get("Sequence length", "NA"),
        "observed_length_min": min((item[1] for item in length_rows), default=None),
        "observed_length_max": max((item[1] for item in length_rows), default=None),
        "dominant_length_bin": max(length_rows, key=lambda item: item[2])[0] if length_rows else "NA",
        "gc_percent": basic.get("%GC", "NA"),
        "per_base_quality_status": quality.get("status", "NA"),
        "quality_min_mean": min((item[1] for item in quality_rows), default=None),
        "quality_last_mean": quality_rows[-1][1] if quality_rows else None,
        "n_content_status": n_module.get("status", "NA"),
        "n_max_percent": max((item[1] for item in n_rows), default=None),
        "n_last_percent": n_rows[-1][1] if n_rows else None,
        "adapter_status": adapter_module.get("status", "NA"),
        "adapter_signals_max_percent": ";".join(adapter_summary) or "none_ge_0.1pct",
        "overrepresented_status": overrep.get("status", "NA"),
        "top_overrepresented_sequences": ";".join(overrep_rows) or "none",
        "duplication_status": modules.get("Sequence Duplication Levels", {}).get("status", "NA"),
        "per_base_content_status": modules.get("Per base sequence content", {}).get("status", "NA"),
    }
    return {key: fmt(value) for key, value in result.items()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fastqc-dir", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    with args.manifest.open(encoding="utf-8-sig", newline="") as handle:
        metadata = {row["filename"]: row for row in csv.DictReader(handle, delimiter="\t")}

    records = [
        parse_archive(path, metadata)
        for path in sorted(args.fastqc_dir.glob("*_fastqc.zip"))
    ]
    if not records:
        raise SystemExit("No FastQC ZIP files found")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(records)
    print(f"records={len(records)} output={args.output}")


if __name__ == "__main__":
    main()
