#!/usr/bin/env python3
"""Combine per-sample model-ready PI tables into analysis matrices."""

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import OrderedDict
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, MutableMapping, Sequence, Tuple


ID_COLUMNS = ["gene", "transcript", "chrom", "strand", "tss_bed0", "tes_bed0"]
VALUE_COLUMNS = ["pi_raw", "log2_pi_pseudocount"]
REQUIRED_COLUMNS = set(ID_COLUMNS + VALUE_COLUMNS)
INPUT_SUFFIX = ".PI.model_ready.tsv"


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Merge model-ready PI tables by MANE transcript."
    )
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--missing-value",
        default="NA",
        help="Text written for a transcript absent from a sample (default: NA)",
    )
    return parser.parse_args(argv)


def chromosome_key(chrom: str) -> Tuple[int, str]:
    match = re.fullmatch(r"chr(\d+|X|Y)", chrom)
    if not match:
        return (1000, chrom)
    token = match.group(1)
    if token.isdigit():
        return (int(token), "")
    return ({"X": 23, "Y": 24}[token], "")


def finite_number(text: str, column: str, path: Path, line_number: int) -> str:
    try:
        value = float(text)
    except ValueError as exc:
        raise ValueError(
            f"{path}:{line_number}: invalid {column} value: {text!r}"
        ) from exc
    if not math.isfinite(value):
        raise ValueError(
            f"{path}:{line_number}: non-finite {column} value: {text!r}"
        )
    return text


def read_sample(
    path: Path,
) -> Tuple[Dict[str, Tuple[str, ...]], Dict[str, Dict[str, str]]]:
    metadata: Dict[str, Tuple[str, ...]] = {}
    values: Dict[str, Dict[str, str]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path}: missing header")
        missing = REQUIRED_COLUMNS.difference(reader.fieldnames)
        if missing:
            raise ValueError(f"{path}: missing columns: {', '.join(sorted(missing))}")
        for line_number, row in enumerate(reader, 2):
            transcript = row["transcript"]
            if not transcript:
                raise ValueError(f"{path}:{line_number}: empty transcript identifier")
            if transcript in values:
                raise ValueError(f"{path}:{line_number}: duplicate transcript {transcript}")
            metadata[transcript] = tuple(row[column] for column in ID_COLUMNS)
            values[transcript] = {
                column: finite_number(row[column], column, path, line_number)
                for column in VALUE_COLUMNS
            }
    return metadata, values


def write_matrix(
    path: Path,
    transcripts: Iterable[str],
    metadata: Mapping[str, Tuple[str, ...]],
    sample_names: Sequence[str],
    sample_values: Mapping[str, Mapping[str, Mapping[str, str]]],
    value_column: str,
    missing_value: str,
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(ID_COLUMNS + list(sample_names))
        for transcript in transcripts:
            row = list(metadata[transcript])
            row.extend(
                sample_values[sample].get(transcript, {}).get(value_column, missing_value)
                for sample in sample_names
            )
            writer.writerow(row)


def write_coverage(
    path: Path,
    transcripts: Iterable[str],
    metadata: Mapping[str, Tuple[str, ...]],
    sample_names: Sequence[str],
    sample_values: Mapping[str, Mapping[str, Mapping[str, str]]],
) -> None:
    sample_total = len(sample_names)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(ID_COLUMNS + ["samples_with_model_ready_pi", "sample_coverage_pct"])
        for transcript in transcripts:
            present = sum(transcript in sample_values[sample] for sample in sample_names)
            coverage = 100 * present / sample_total if sample_total else 0
            writer.writerow(list(metadata[transcript]) + [present, f"{coverage:.2f}"])


def write_long_table(
    path: Path,
    transcripts: Iterable[str],
    metadata: Mapping[str, Tuple[str, ...]],
    sample_names: Sequence[str],
    sample_values: Mapping[str, Mapping[str, Mapping[str, str]]],
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample"] + ID_COLUMNS + VALUE_COLUMNS)
        for sample in sample_names:
            values = sample_values[sample]
            for transcript in transcripts:
                if transcript not in values:
                    continue
                writer.writerow(
                    [sample]
                    + list(metadata[transcript])
                    + [values[transcript][column] for column in VALUE_COLUMNS]
                )


def build_matrices(input_dir: Path, output_dir: Path, missing_value: str) -> None:
    paths = sorted(input_dir.glob(f"*{INPUT_SUFFIX}"))
    if not paths:
        raise ValueError(f"No *{INPUT_SUFFIX} files found in {input_dir}")

    metadata: MutableMapping[str, Tuple[str, ...]] = OrderedDict()
    sample_values: Dict[str, Dict[str, Dict[str, str]]] = {}
    sample_names: List[str] = []

    for path in paths:
        sample = path.name[: -len(INPUT_SUFFIX)]
        if sample in sample_values:
            raise ValueError(f"Duplicate sample name: {sample}")
        sample_metadata, values = read_sample(path)
        for transcript, current in sample_metadata.items():
            previous = metadata.get(transcript)
            if previous is not None and previous != current:
                raise ValueError(
                    f"Inconsistent annotation for transcript {transcript}: "
                    f"{previous!r} versus {current!r} in {path}"
                )
            metadata.setdefault(transcript, current)
        sample_names.append(sample)
        sample_values[sample] = values

    transcripts = sorted(
        metadata,
        key=lambda tx: (
            chromosome_key(metadata[tx][2]),
            int(metadata[tx][4]),
            metadata[tx][0],
            tx,
        ),
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    write_matrix(
        output_dir / "human_transcript_log2_pi_matrix.tsv",
        transcripts,
        metadata,
        sample_names,
        sample_values,
        "log2_pi_pseudocount",
        missing_value,
    )
    write_matrix(
        output_dir / "human_transcript_pi_raw_matrix.tsv",
        transcripts,
        metadata,
        sample_names,
        sample_values,
        "pi_raw",
        missing_value,
    )
    write_coverage(
        output_dir / "human_transcript_pi_coverage.tsv",
        transcripts,
        metadata,
        sample_names,
        sample_values,
    )
    write_long_table(
        output_dir / "human_transcript_pi_long.tsv",
        transcripts,
        metadata,
        sample_names,
        sample_values,
    )
    print(f"Samples: {len(sample_names)}")
    print(f"Union of transcripts: {len(transcripts)}")
    print(f"Output directory: {output_dir}")


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        build_matrices(args.input_dir, args.output_dir, args.missing_value)
        return 0
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
