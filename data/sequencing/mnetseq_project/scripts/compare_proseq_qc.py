#!/usr/bin/env python3
"""Compare raw and cleaned FastQC metrics one file at a time."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def number(value, default=0.0):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def adapter_max(text):
    maximum = 0.0
    if not text or text.startswith("none_"):
        return maximum
    for item in text.split(";"):
        if ":" not in item:
            continue
        name, value = item.rsplit(":", 1)
        if name in {"PolyA", "PolyG"}:
            continue
        maximum = max(maximum, number(value))
    return maximum


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-metrics", required=True, type=Path)
    parser.add_argument("--clean-metrics", required=True, type=Path)
    parser.add_argument("--cleaning-manifest", required=True, type=Path)
    parser.add_argument("--cleaning-status", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    raw = {row["filename"]: row for row in read_tsv(args.raw_metrics)}
    clean = {row["filename"]: row for row in read_tsv(args.clean_metrics)}
    status = {row["run"]: row for row in read_tsv(args.cleaning_status)}

    records = []
    for run_row in read_tsv(args.cleaning_manifest):
        run = run_row["run"]
        run_status = status.get(run, {})
        pairs = [(run_row["input_r1"], run_row["output_r1"], "1" if run_row["layout"] == "PAIRED" else "0")]
        if run_row["layout"] == "PAIRED":
            pairs.append((run_row["input_r2"], run_row["output_r2"], "2"))

        for raw_path, clean_path, mate in pairs:
            raw_name = Path(raw_path).name
            clean_name = Path(clean_path).name
            before = raw.get(raw_name)
            after = clean.get(clean_name)
            reasons = []
            if before is None or after is None:
                verdict = "FAIL_MISSING_QC"
            else:
                after_adapter = adapter_max(after["adapter_signals_max_percent"])
                after_last_q = number(after["quality_last_mean"])
                after_n = number(after["n_max_percent"])
                retention = number(run_status.get("retention_percent"), 100.0)
                if after_adapter > 5:
                    reasons.append(f"residual_adapter={after_adapter:.2f}%")
                if after_last_q < 20:
                    reasons.append(f"last_mean_Q={after_last_q:.2f}")
                if after_n > 5:
                    reasons.append(f"max_N={after_n:.2f}%")
                if retention < 30:
                    reasons.append(f"retention={retention:.2f}%")
                verdict = "PASS" if not reasons else "REVIEW"

            records.append(
                {
                    "gse": run_row["gse"],
                    "run": run,
                    "layout": run_row["layout"],
                    "mate": mate,
                    "raw_filename": raw_name,
                    "clean_filename": clean_name,
                    "retention_percent": run_status.get("retention_percent", "NA"),
                    "raw_length": before["basic_length"] if before else "NA",
                    "clean_length": after["basic_length"] if after else "NA",
                    "raw_last_mean_q": before["quality_last_mean"] if before else "NA",
                    "clean_last_mean_q": after["quality_last_mean"] if after else "NA",
                    "raw_max_n_percent": before["n_max_percent"] if before else "NA",
                    "clean_max_n_percent": after["n_max_percent"] if after else "NA",
                    "raw_adapter_max_percent": f"{adapter_max(before['adapter_signals_max_percent']):.4f}" if before else "NA",
                    "clean_adapter_max_percent": f"{adapter_max(after['adapter_signals_max_percent']):.4f}" if after else "NA",
                    "raw_adapter_status": before["adapter_status"] if before else "NA",
                    "clean_adapter_status": after["adapter_status"] if after else "NA",
                    "verdict": verdict,
                    "review_reason": ";".join(reasons) if reasons else "none",
                }
            )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(records)
    print(f"files={len(records)} output={args.output}")


if __name__ == "__main__":
    main()
