#!/usr/bin/env python3
"""Build protocol-aware, FastQC-informed cleaning rules for mNET-seq Runs."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


SMALL_RNA_R1 = "TGGAATTCTCGGGTGCCAAGG"
SMALL_RNA_R2 = "GATCGTCGGACTGTAGAACTCTGAAC"
ILLUMINA_UNIVERSAL_R1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"


PROTOCOLS = {
    "GSE148433": {
        "adapter_r1": "TGGAATTCTCGG",
        "adapter_r2": "GATCGTCGGACT",
        "overlap": "12",
        "minimum_length": "25",
        "error_rate": "0.10",
        "match_read_wildcards": "no",
        "source_rule": "published_Cutadapt_O12_m25",
    },

    "GSE81662": {
        "adapter_r1": SMALL_RNA_R1,
        "adapter_r2": SMALL_RNA_R2,
        "overlap": "3",
        "minimum_length": "10",
        "error_rate": "0.10",
        "match_read_wildcards": "no",
        "source_rule": "published_full_smallRNA_adapters_m10",
    },
    "GSE138356": {
        "adapter_r1": ILLUMINA_UNIVERSAL_R1,
        "adapter_r2": SMALL_RNA_R2,
        "overlap": "5",
        "minimum_length": "10",
        "error_rate": "0.05",
        "match_read_wildcards": "yes",
        "source_rule": "FastQC_R1_universal_R2_smallRNA5_general_mNET_pipeline",
    },
    "GSE115290": {
        "adapter_r1": ILLUMINA_UNIVERSAL_R1,
        "adapter_r2": SMALL_RNA_R2,
        "overlap": "5",
        "minimum_length": "10",
        "error_rate": "0.05",
        "match_read_wildcards": "yes",
        "source_rule": "FastQC_TruSeq_R1_smallRNA5_R2_corrected_v2",
        "output_version": "v2",
    },
}

DEFAULT_PROTOCOL = {
    "adapter_r1": SMALL_RNA_R1,
    "adapter_r2": SMALL_RNA_R2,
    "overlap": "5",
    "minimum_length": "10",
    "error_rate": "0.05",
    "match_read_wildcards": "yes",
    "source_rule": "NEBNext_smallRNA_protocol_FastQC_verified",
}


NON_MNET_RUNS = {"SRR14777732", "SRR14777733", "SRR14777736", "SRR14777737"}


def read_tsv(path: Path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def number(row: dict[str, str], key: str, default: float = 99.0) -> float:
    try:
        return float(row.get(key, "NA"))
    except (TypeError, ValueError):
        return default


def quality_rule(group_metrics: list[dict[str, str]]) -> tuple[str, str]:
    """Only trim clearly poor 3-prime tails; preserve mNET read-2 5-prime ends."""
    worst_last = min(number(item, "quality_last_mean") for item in group_metrics)
    worst_mean = min(number(item, "quality_min_mean") for item in group_metrics)
    if min(worst_last, worst_mean) < 15:
        return "quality_3prime", "15"
    if min(worst_last, worst_mean) < 20:
        return "quality_3prime", "10"
    return "none", "0"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--metrics", required=True, type=Path)
    parser.add_argument("--project", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    manifest = read_tsv(args.manifest)
    metrics = {row["filename"]: row for row in read_tsv(args.metrics)}
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    excluded_runs: set[str] = set()
    for row in manifest:


        if row["run"] in NON_MNET_RUNS:
            excluded_runs.add(row["run"])
            continue
        grouped[row["run"]].append(row)

    output_rows = []
    for run, rows in sorted(grouped.items()):
        rows.sort(key=lambda row: int(row["mate"]))
        first = rows[0]
        gse = first["gse"]
        layout = first["layout"].upper()
        protocol = PROTOCOLS.get(gse, DEFAULT_PROTOCOL)
        by_mate = {row["mate"]: row for row in rows}

        if layout == "PAIRED":
            r1, r2 = by_mate["1"], by_mate["2"]
            filenames = [r1["filename"], r2["filename"]]
            input_r1 = f"{args.project}/raw/{r1['filename']}"
            input_r2 = f"{args.project}/raw/{r2['filename']}"
            version = f".{protocol['output_version']}" if protocol.get("output_version") else ""
            output_r1 = f"{args.project}/clean/{gse}/{run}_1{version}.clean.fastq.gz"
            output_r2 = f"{args.project}/clean/{gse}/{run}_2{version}.clean.fastq.gz"
            adapter_r2 = protocol["adapter_r2"]
        else:
            r1 = by_mate.get("0", rows[0])
            filenames = [r1["filename"]]
            input_r1 = f"{args.project}/raw/{r1['filename']}"
            input_r2 = "NA"
            output_r1 = f"{args.project}/clean/{gse}/{run}.clean.fastq.gz"
            output_r2 = "NA"
            adapter_r2 = "NA"

        missing_metrics = [name for name in filenames if name not in metrics]
        if missing_metrics:
            raise ValueError(f"Missing FastQC metrics for {run}: {missing_metrics}")
        group_metrics = [metrics[name] for name in filenames]
        quality_mode, quality_cutoff = quality_rule(group_metrics)

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
                "adapter_r1": protocol["adapter_r1"],
                "adapter_r2": adapter_r2,
                "adapter_overlap": protocol["overlap"],
                "error_rate": protocol["error_rate"],
                "match_read_wildcards": protocol["match_read_wildcards"],
                "quality_mode": quality_mode,
                "quality_cutoff": quality_cutoff,
                "trim_terminal_n": "yes",
                "max_n_fraction": "0.05",
                "minimum_length": protocol["minimum_length"],
                "umi_action": "retain_no_inline_UMI_confirmed",
                "mnet_coordinate_note": "paired_R2_5prime_base_with_R1_strand",
                "source_rule": protocol["source_rule"],
                "raw_adapter_evidence": ";".join(
                    f"{item['filename']}={item['adapter_signals_max_percent']}"
                    for item in group_metrics
                ),
                "raw_quality_evidence": ";".join(
                    f"{item['filename']}:minQ={item['quality_min_mean']},lastQ={item['quality_last_mean']}"
                    for item in group_metrics
                ),
                "raw_n_evidence": ";".join(
                    f"{item['filename']}:maxN={item['n_max_percent']}%"
                    for item in group_metrics
                ),
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)
    print(
        f"mnet_runs={len(output_rows)} excluded_non_mnet_runs={len(excluded_runs)} "
        f"excluded={','.join(sorted(excluded_runs)) or 'none'} output={args.output}"
    )


if __name__ == "__main__":
    main()
