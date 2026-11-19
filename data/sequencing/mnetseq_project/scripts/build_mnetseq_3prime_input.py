#!/usr/bin/env python3
"""Join mNET-seq metadata to Q20 BAM results for 3-prime extraction."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def formal_names(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for row in read_tsv(path):
        result[row["old_sample_id"]] = row["formal_sample_name"]
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cleaning-manifest", required=True, type=Path)
    parser.add_argument("--q20-summary", required=True, type=Path)
    parser.add_argument("--formal-names", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--exclude-runs",
        nargs="*",
        default=["SRR11521583", "SRR11521584"],
    )
    args = parser.parse_args()

    metadata = {row["run"]: row for row in read_tsv(args.cleaning_manifest)}
    names = formal_names(args.formal_names)
    excluded = set(args.exclude_runs)
    output_rows: list[dict[str, str]] = []

    for q20 in read_tsv(args.q20_summary):
        run = q20["run"]
        if run in excluded:
            continue
        if q20.get("status") != "PASS":
            continue
        if run not in metadata:
            raise ValueError(f"Q20 Run is absent from cleaning manifest: {run}")
        meta = metadata[run]
        if meta["layout"].upper() != "PAIRED":
            raise ValueError(f"mNET 3-prime extraction requires paired reads: {run}")

        q20_records = int(q20["q20_records"])
        pairs_text = q20.get("q20_pairs", "").strip()
        if pairs_text:
            retained_pairs = int(pairs_text)
        else:
            if q20_records % 2:
                raise ValueError(f"Odd Q20 record count for paired Run: {run}")
            retained_pairs = q20_records // 2

        old_sample_id = f"{meta['gse']}_{meta['gsm']}"
        output_rows.append(
            {
                "run": run,
                "gse": meta["gse"],
                "gsm": meta["gsm"],
                "sample_label": names.get(old_sample_id, old_sample_id),
                "retained_pairs": str(retained_pairs),
                "output_bam": q20["output_bam"],
                "status": "PASS",
                "note": "Q20_proper_pair_ready_for_R2_5prime_extraction",
            }
        )

    if not output_rows:
        raise ValueError("No PASS Q20 mNET-seq Runs remain after exclusions")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(output_rows[0]), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(sorted(output_rows, key=lambda row: row["run"]))

    sample_count = len({(row["gse"], row["gsm"]) for row in output_rows})
    print(
        f"3-prime inputs: {len(output_rows)} Runs, {sample_count} samples; "
        f"excluded={','.join(sorted(excluded)) or 'none'}"
    )


if __name__ == "__main__":
    main()
