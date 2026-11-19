#!/usr/bin/env python3
"""Expand one-row-per-Run cleaning metadata to one-row-per-clean-FASTQ."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cleaning-manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    with args.cleaning_manifest.open(encoding="utf-8-sig", newline="") as handle:
        runs = list(csv.DictReader(handle, delimiter="\t"))

    rows = []
    for run in runs:
        files = [(run["output_r1"], "1" if run["layout"] == "PAIRED" else "0")]
        if run["layout"] == "PAIRED":
            files.append((run["output_r2"], "2"))
        for path, mate in files:
            rows.append(
                {
                    "method": "mNET-seq_clean",
                    "group_id": run["group_id"],
                    "gse": run["gse"],
                    "gsm": run["gsm"],
                    "run": run["run"],
                    "layout": run["layout"],
                    "mate": mate,
                    "filename": Path(path).name,
                }
            )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"clean_fastq_rows={len(rows)} output={args.output}")


if __name__ == "__main__":
    main()
