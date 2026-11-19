#!/usr/bin/env python3
"""Audit ENA FASTQ availability against a run manifest."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--ena", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    manifest_rows = read_tsv(args.manifest)
    ena_rows = read_tsv(args.ena)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    ena_by_run = {row.get("run_accession", ""): row for row in ena_rows}
    ready: list[dict[str, str]] = []
    unavailable: list[dict[str, str]] = []
    issues: list[dict[str, str]] = []
    total_bytes = 0

    fields = [
        "run_accession",
        "provisional_sample_id",
        "expected_method",
        "expected_species",
        "expected_cell_type",
        "library_strategy",
        "library_layout",
        "scientific_name",
        "sample_accession",
        "sample_title",
        "experiment_title",
        "fastq_ftp",
        "fastq_md5",
        "fastq_bytes",
        "status",
        "issue",
    ]

    for manifest in manifest_rows:
        run = manifest["run_accession"]
        ena = ena_by_run.get(run)
        base = {
            "run_accession": run,
            "provisional_sample_id": manifest.get("provisional_sample_id", ""),
            "expected_method": manifest.get("method", ""),
            "expected_species": manifest.get("species", ""),
            "expected_cell_type": manifest.get("cell_type", ""),
        }

        if not ena:
            unavailable.append(
                {**base, **{f: "" for f in fields if f not in base},
                 "status": "unavailable", "issue": "No ENA read_run record returned"}
            )
            continue

        record = {
            **base,
            "library_strategy": ena.get("library_strategy", ""),
            "library_layout": ena.get("library_layout", ""),
            "scientific_name": ena.get("scientific_name", ""),
            "sample_accession": ena.get("sample_accession", ""),
            "sample_title": ena.get("sample_title", ""),
            "experiment_title": ena.get("experiment_title", ""),
            "fastq_ftp": ena.get("fastq_ftp", ""),
            "fastq_md5": ena.get("fastq_md5", ""),
            "fastq_bytes": ena.get("fastq_bytes", ""),
            "status": "ready",
            "issue": "",
        }

        problems: list[str] = []
        if not record["fastq_ftp"]:
            record["status"] = "unavailable"
            record["issue"] = "ENA record exists but fastq_ftp is empty"
            unavailable.append(record)
            continue
        if record["scientific_name"] and record["scientific_name"] != manifest.get(
            "species", ""
        ):
            problems.append("species mismatch")
        if record["library_layout"] not in {"SINGLE", "PAIRED"}:
            problems.append("unknown library layout")
        if record["library_strategy"] and record["library_strategy"] not in {
            "OTHER",
            "RNA-Seq",
        }:
            problems.append("unexpected ENA library strategy")

        byte_values = [x for x in record["fastq_bytes"].split(";") if x]
        try:
            total_bytes += sum(int(x) for x in byte_values)
        except ValueError:
            problems.append("invalid fastq_bytes")

        if problems:
            record["status"] = "review"
            record["issue"] = "; ".join(problems)
            issues.append(record)
        else:
            ready.append(record)

    write_tsv(args.output_dir / "ena_fastq_ready.tsv", ready, fields)
    write_tsv(args.output_dir / "ena_fastq_unavailable.tsv", unavailable, fields)
    write_tsv(args.output_dir / "ena_metadata_review.tsv", issues, fields)

    gib = total_bytes / 1024**3
    print(f"Manifest runs: {len(manifest_rows)}")
    print(f"Direct FASTQ ready: {len(ready)}")
    print(f"Unavailable: {len(unavailable)}")
    print(f"Needs metadata review: {len(issues)}")
    print(f"Reported FASTQ volume: {gib:.2f} GiB")


if __name__ == "__main__":
    main()
