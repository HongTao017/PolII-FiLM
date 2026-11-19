#!/usr/bin/env python3
"""Combine mNET-seq stage reports into run- and sample-level audit tables."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cleaning", required=True, type=Path)
    parser.add_argument("--rrna", required=True, type=Path)
    parser.add_argument("--alignment", required=True, type=Path)
    parser.add_argument("--three-prime", required=True, type=Path)
    parser.add_argument("--sample-merge", required=True, type=Path)
    parser.add_argument("--pi-summary", required=True, type=Path)
    parser.add_argument("--run-output", required=True, type=Path)
    parser.add_argument("--sample-output", required=True, type=Path)
    return parser.parse_args()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def index_rows(path: Path, key: str) -> dict[str, dict[str, str]]:
    return {row[key]: row for row in read_rows(path)}


def write_rows(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    cleaning = index_rows(args.cleaning, "run")
    rrna = index_rows(args.rrna, "run")
    alignment = index_rows(args.alignment, "run")
    three_prime = index_rows(args.three_prime, "run")

    run_ids = sorted(set(cleaning) | set(rrna) | set(alignment) | set(three_prime))
    run_rows: list[dict[str, str]] = []
    for run in run_ids:
        clean = cleaning.get(run, {})
        ribo = rrna.get(run, {})
        align = alignment.get(run, {})
        end = three_prime.get(run, {})
        stage_statuses = [
            ribo.get("status", "MISSING"),
            align.get("status", "MISSING"),
            end.get("status", "MISSING"),
        ]
        run_rows.append(
            {
                "run": run,
                "gse": align.get("gse", clean.get("gse", "NA")),
                "gsm": align.get("gsm", clean.get("gsm", "NA")),
                "sample_label": align.get("sample_label", "NA"),
                "layout": clean.get("layout", ribo.get("layout", "NA")),
                "adapter_r1": clean.get("adapter_r1", "NA"),
                "adapter_r2": clean.get("adapter_r2", "NA"),
                "rrna_input_fragments": ribo.get("input_fragments", "NA"),
                "non_rrna_fragments": ribo.get("retained_fragments", "NA"),
                "rrna_percent": ribo.get("rrna_percent", "NA"),
                "star_unique_percent": align.get("star_unique_percent", "NA"),
                "retained_unique_proper_pairs": align.get("retained_pairs", "NA"),
                "adapter_motif": end.get("adapter_motif", "NA"),
                "retained_3prime_ends": end.get("retained_3prime_ends", "NA"),
                "internal_priming_percent": end.get("internal_priming_percent", "NA"),
                "output_bam": align.get("output_bam", "NA"),
                "output_3prime_bed": end.get("output_bed_gz", "NA"),
                "pipeline_status": "PASS" if all(x == "PASS" for x in stage_statuses) else "CHECK",
            }
        )
    write_rows(args.run_output, run_rows)

    merged = index_rows(args.sample_merge, "sample_id")
    pi = index_rows(args.pi_summary, "sample_id")
    sample_ids = sorted(set(merged) | set(pi))
    sample_rows: list[dict[str, str]] = []
    for sample_id in sample_ids:
        sample = merged.get(sample_id, {})
        pi_row = pi.get(sample_id, {})
        statuses = [sample.get("status", "MISSING"), pi_row.get("status", "MISSING")]
        sample_rows.append(
            {
                "sample_id": sample_id,
                "gse": sample.get("gse", pi_row.get("gse", "NA")),
                "gsm": sample.get("gsm", pi_row.get("gsm", "NA")),
                "sample_label": sample.get("sample_label", "NA"),
                "technical_run_count": sample.get("technical_run_count", "NA"),
                "runs": sample.get("runs", "NA"),
                "total_3prime_ends": sample.get("total_3prime_ends", "NA"),
                "total_mane_transcripts": pi_row.get("total_transcripts", "NA"),
                "active_min10_transcripts": pi_row.get("active_min10_transcripts", "NA"),
                "output_3prime_bed": sample.get("output_bed_gz", "NA"),
                "output_pi": pi_row.get("output_pi", "NA"),
                "pipeline_status": "PASS" if all(x == "PASS" for x in statuses) else "CHECK",
            }
        )
    write_rows(args.sample_output, sample_rows)
    print(f"Run audit: {len(run_rows)} rows -> {args.run_output}")
    print(f"Sample audit: {len(sample_rows)} rows -> {args.sample_output}")


if __name__ == "__main__":
    main()
