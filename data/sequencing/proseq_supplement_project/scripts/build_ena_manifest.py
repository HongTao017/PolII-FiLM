#!/usr/bin/env python3
"""Build a verified ENA FASTQ download manifest for supplemental PRO-seq runs."""

from __future__ import annotations

import argparse
import csv
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path


OUTPUT_COLUMNS = [
    "method", "group_id", "cell_types", "gse", "reference", "sample_label",
    "gsm", "run", "layout", "reported_layout", "mate", "filename", "url",
    "md5", "bytes", "scientific_name", "sample_accession", "experiment_accession",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {
        "group_id", "method", "species", "cell_type", "reference", "gse",
        "run", "formal_sample",
    }
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"Sample table is missing columns: {sorted(required)}")
    return rows


def ena_record(run: str) -> dict[str, str]:
    fields = [
        "run_accession", "experiment_accession", "sample_accession",
        "scientific_name", "library_layout", "fastq_ftp", "fastq_md5",
        "fastq_bytes",
    ]
    query = urllib.parse.urlencode(
        {
            "accession": run,
            "result": "read_run",
            "fields": ",".join(fields),
            "format": "json",
        }
    )
    url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?{query}"
    request = urllib.request.Request(url, headers={"User-Agent": "PROseqSupplement/1.0"})
    last_error: Exception | None = None
    for attempt in range(5):
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                records = json.load(response)
            if len(records) != 1:
                raise ValueError(f"ENA returned {len(records)} records for {run}")
            return records[0]
        except Exception as exc:  # network failures are retried with backoff
            last_error = exc
            time.sleep(2 ** attempt)
    raise RuntimeError(f"ENA query failed for {run}: {last_error}")


def split_field(record: dict[str, str], name: str) -> list[str]:
    value = record.get(name, "")
    return value.split(";") if value else []


def main() -> None:
    args = parse_args()
    output_rows: list[dict[str, str]] = []
    for sample in read_tsv(args.samples):
        run = sample["run"]
        record = ena_record(run)
        layout = record["library_layout"].upper()
        urls = split_field(record, "fastq_ftp")
        md5s = split_field(record, "fastq_md5")
        sizes = split_field(record, "fastq_bytes")
        expected = 2 if layout == "PAIRED" else 1
        if not (len(urls) == len(md5s) == len(sizes) == expected):
            raise ValueError(
                f"{run}: layout={layout}, FASTQ={len(urls)}, MD5={len(md5s)}, bytes={len(sizes)}"
            )
        for index, (ftp, md5, size) in enumerate(zip(urls, md5s, sizes), 1):
            filename = ftp.rsplit("/", 1)[-1]
            output_rows.append(
                {
                    "method": sample["method"],
                    "group_id": sample["group_id"],
                    "cell_types": sample["cell_type"],
                    "gse": sample["gse"],
                    "reference": sample["reference"],
                    "sample_label": sample["formal_sample"],
                    "gsm": "NA",
                    "run": run,
                    "layout": layout,
                    "reported_layout": layout,
                    "mate": str(index if layout == "PAIRED" else 0),
                    "filename": filename,
                    "url": "https://" + ftp.removeprefix("ftp://").removeprefix("https://"),
                    "md5": md5,
                    "bytes": size,
                    "scientific_name": record["scientific_name"],
                    "sample_accession": record["sample_accession"],
                    "experiment_accession": record["experiment_accession"],
                }
            )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(output_rows)
    temporary.replace(args.output)
    print(f"Runs: {len({row['run'] for row in output_rows})}")
    print(f"FASTQ files: {len(output_rows)}")
    print(f"Manifest: {args.output}")


if __name__ == "__main__":
    main()
