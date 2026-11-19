#!/usr/bin/env python3
"""Append an R2-only mNET-seq sample to the main PI matrix by MANE ID."""

from __future__ import annotations

import argparse
import csv
import os
from collections import OrderedDict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-matrix", required=True, type=Path)
    parser.add_argument("--extra-matrix", required=True, type=Path)
    parser.add_argument("--base-summary", required=True, type=Path)
    parser.add_argument("--extra-summary", required=True, type=Path)
    parser.add_argument("--output-matrix", required=True, type=Path)
    parser.add_argument("--output-summary", required=True, type=Path)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def read_matrix(path: Path) -> tuple[str, list[str], OrderedDict[str, list[str]]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError(f"Empty matrix: {path}") from exc
        if len(header) < 2:
            raise ValueError(f"Matrix has no sample columns: {path}")
        key_name, samples = header[0], header[1:]
        if len(samples) != len(set(samples)):
            raise ValueError(f"Duplicate sample columns in {path}")
        rows: OrderedDict[str, list[str]] = OrderedDict()
        for line_number, row in enumerate(reader, 2):
            if not row:
                continue
            if len(row) != len(header):
                raise ValueError(
                    f"Wrong column count in {path}, line {line_number}: "
                    f"expected {len(header)}, found {len(row)}"
                )
            key = row[0]
            if not key:
                raise ValueError(f"Empty MANE key in {path}, line {line_number}")
            if key in rows:
                raise ValueError(f"Duplicate MANE key in {path}: {key}")
            rows[key] = [value if value else "NA" for value in row[1:]]
    return key_name, samples, rows


def read_summary(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Summary has no header: {path}")
        return list(reader.fieldnames), list(reader)


def atomic_tsv(path: Path, rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(f"{path.name}.incomplete.{os.getpid()}")
    with temp.open("w", newline="", encoding="utf-8") as handle:
        csv.writer(handle, delimiter="\t", lineterminator="\n").writerows(rows)
    temp.replace(path)


def main() -> None:
    args = parse_args()
    for path in (
        args.base_matrix,
        args.extra_matrix,
        args.base_summary,
        args.extra_summary,
    ):
        if not path.is_file():
            raise FileNotFoundError(path)
    for path in (args.output_matrix, args.output_summary):
        if path.exists() and not args.force:
            raise FileExistsError(f"Output exists; use --force to replace it: {path}")

    base_key, base_samples, base_rows = read_matrix(args.base_matrix)
    extra_key, extra_samples, extra_rows = read_matrix(args.extra_matrix)
    if base_key != extra_key:
        raise ValueError(f"Matrix key headers differ: {base_key!r} vs {extra_key!r}")
    overlap = sorted(set(base_samples) & set(extra_samples))
    if overlap:
        raise ValueError(f"Sample columns already exist in both matrices: {overlap}")

    union_keys = list(base_rows)
    union_keys.extend(key for key in extra_rows if key not in base_rows)
    matrix_output = [[base_key, *base_samples, *extra_samples]]
    for key in union_keys:
        matrix_output.append(
            [
                key,
                *(base_rows.get(key) or ["NA"] * len(base_samples)),
                *(extra_rows.get(key) or ["NA"] * len(extra_samples)),
            ]
        )
    atomic_tsv(args.output_matrix, matrix_output)

    base_fields, base_summary = read_summary(args.base_summary)
    extra_fields, extra_summary = read_summary(args.extra_summary)
    summary_fields = list(base_fields)
    summary_fields.extend(field for field in extra_fields if field not in summary_fields)
    sample_field = "sample_id" if "sample_id" in summary_fields else summary_fields[0]
    summary_rows = base_summary + extra_summary
    sample_ids = [row.get(sample_field, "") for row in summary_rows]
    if "" in sample_ids or len(sample_ids) != len(set(sample_ids)):
        raise ValueError(f"Missing or duplicate {sample_field} values in combined summary")
    args.output_summary.parent.mkdir(parents=True, exist_ok=True)
    summary_temp = args.output_summary.with_name(
        f"{args.output_summary.name}.incomplete.{os.getpid()}"
    )
    with summary_temp.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=summary_fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in summary_rows:
            writer.writerow({field: row.get(field, "NA") or "NA" for field in summary_fields})
    summary_temp.replace(args.output_summary)

    base_only = sum(key not in extra_rows for key in base_rows)
    extra_only = sum(key not in base_rows for key in extra_rows)
    shared = sum(key in extra_rows for key in base_rows)
    print(f"Base samples: {len(base_samples)}")
    print(f"Added samples: {len(extra_samples)} ({', '.join(extra_samples)})")
    print(f"Combined samples: {len(base_samples) + len(extra_samples)}")
    print(f"Base transcripts: {len(base_rows)}")
    print(f"Extra transcripts: {len(extra_rows)}")
    print(f"Shared transcripts: {shared}")
    print(f"Base-only transcripts: {base_only}")
    print(f"Extra-only transcripts: {extra_only}")
    print(f"Combined transcript union: {len(union_keys)}")
    print(f"Combined matrix: {args.output_matrix}")
    print(f"Combined summary: {args.output_summary}")


if __name__ == "__main__":
    main()
