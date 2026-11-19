#!/usr/bin/env python3
"""Build one conservative, FastQC-informed cleaning rule per PRO-seq run."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


SMALL_RNA_3 = "TGGAATTCTCGGGTGCCAAGG"
SMALL_RNA_5 = "GATCGTCGGACTGTAGAACTCTGAAC"
ILLUMINA_UNIVERSAL = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"

UNIVERSAL_GSE = {"GSE173783", "GSE233637", "GSE241764"}
NEXTSEQ_GSE = {"GSE179916", "GSE241764", "GSE296252"}


def read_tsv(path: Path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--metrics", required=True, type=Path)
    parser.add_argument("--project", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    manifest = read_tsv(args.manifest)
    metrics = {row["filename"]: row for row in read_tsv(args.metrics)}
    grouped = defaultdict(list)
    for row in manifest:
        grouped[row["run"]].append(row)

    output_rows = []
    for run, rows in sorted(grouped.items()):
        rows.sort(key=lambda row: int(row["mate"]))
        first = rows[0]
        gse = first["gse"]
        layout = first["layout"]
        by_mate = {row["mate"]: row for row in rows}

        if layout == "PAIRED":
            r1 = by_mate["1"]
            r2 = by_mate["2"]
            input_r1 = f"{args.project}/raw/{r1['filename']}"
            input_r2 = f"{args.project}/raw/{r2['filename']}"
            output_r1 = f"{args.project}/clean/{gse}/{run}_1.clean.fastq.gz"
            output_r2 = f"{args.project}/clean/{gse}/{run}_2.clean.fastq.gz"
            if gse in UNIVERSAL_GSE:
                adapter_r1 = adapter_r2 = ILLUMINA_UNIVERSAL
                adapter_type = "Illumina_universal_both_mates"
            else:
                adapter_r1, adapter_r2 = SMALL_RNA_3, SMALL_RNA_5
                adapter_type = "small_RNA_3prime_R1_and_5prime_R2"
            filenames = [r1["filename"], r2["filename"]]
        else:
            r1 = by_mate.get("0", rows[0])
            input_r1 = f"{args.project}/raw/{r1['filename']}"
            input_r2 = "NA"
            output_r1 = f"{args.project}/clean/{gse}/{run}.clean.fastq.gz"
            output_r2 = "NA"
            if gse in UNIVERSAL_GSE:
                adapter_r1 = ILLUMINA_UNIVERSAL
                adapter_type = "Illumina_universal"
            else:
                adapter_r1 = SMALL_RNA_3
                adapter_type = "small_RNA_3prime"
            adapter_r2 = "NA"
            filenames = [r1["filename"]]

        group_metrics = [metrics[name] for name in filenames]
        adapter_signals = ";".join(
            f"{item['filename']}={item['adapter_signals_max_percent']}"
            for item in group_metrics
        )
        quality_signals = ";".join(
            f"{item['filename']}:minQ={item['quality_min_mean']},lastQ={item['quality_last_mean']}"
            for item in group_metrics
        )
        n_signals = ";".join(
            f"{item['filename']}:maxN={item['n_max_percent']}%" for item in group_metrics
        )

        quality_mode = "nextseq_3prime" if gse in NEXTSEQ_GSE else "quality_3prime"
        output_rows.append(
            {
                "gse": gse,
                "group_id": first["group_id"],
                "gsm": first["gsm"],
                "run": run,
                "layout": layout,
                "input_r1": input_r1,
                "input_r2": input_r2,
                "output_r1": output_r1,
                "output_r2": output_r2,
                "adapter_type": adapter_type,
                "adapter_r1": adapter_r1,
                "adapter_r2": adapter_r2,
                "quality_mode": quality_mode,
                "quality_cutoff": "20",
                "trim_terminal_n": "yes",
                "max_n_fraction": "0.05",
                "minimum_length": "20",
                "poly_a_action": "retain_no_strong_signal",
                "umi_action": "retain_for_protocol_specific_extraction_before_alignment",
                "raw_adapter_evidence": adapter_signals,
                "raw_quality_evidence": quality_signals,
                "raw_n_evidence": n_signals,
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)
    print(f"runs={len(output_rows)} output={args.output}")


if __name__ == "__main__":
    main()
