#!/usr/bin/env python3
"""Build a run-level mNET-seq alignment manifest from audited metadata."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis-manifest", required=True, type=Path)
    parser.add_argument("--rrna-status", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def read_run_metadata(path: Path) -> dict[str, dict[str, str]]:
    runs: dict[str, dict[str, str]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            run = row["run"]
            current = {
                "run": run,
                "gse": row["gse"],
                "gsm": row["gsm"],
                "group_id": row["group_id"],
                "sample_label": row["sample_label"],
            }
            if run in runs and runs[run] != current:
                raise ValueError(f"Conflicting metadata for {run}")
            runs[run] = current
    return runs


def main() -> None:
    args = parse_args()
    metadata = read_run_metadata(args.analysis_manifest)
    output_rows: list[dict[str, str]] = []

    with args.rrna_status.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            run = row["run"]
            if row["status"] != "PASS":
                raise ValueError(f"rRNA filtering did not pass for {run}")
            if run not in metadata:
                raise ValueError(f"No audited metadata found for {run}")
            if row["layout"] != "PAIRED":
                raise ValueError(f"mNET-seq run {run} is not paired-end")
            output_rows.append(
                {
                    **metadata[run],
                    "input_r1": row["output_r1"],
                    "input_r2": row["output_r2"],
                    "rrna_percent": row["rrna_percent"],
                }
            )

    if len(output_rows) != len(metadata):
        seen = {row["run"] for row in output_rows}
        missing = sorted(set(metadata) - seen)
        raise ValueError(f"Missing rRNA-filtered runs: {', '.join(missing)}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "run",
        "gse",
        "gsm",
        "group_id",
        "sample_label",
        "input_r1",
        "input_r2",
        "rrna_percent",
    ]
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(sorted(output_rows, key=lambda row: (row["gse"], row["gsm"], row["run"])))

    print(f"Wrote {len(output_rows)} paired mNET-seq runs to {args.output}")


if __name__ == "__main__":
    main()
