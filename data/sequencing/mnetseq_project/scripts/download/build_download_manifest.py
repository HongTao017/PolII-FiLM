#!/usr/bin/env python3
"""Build a nascent-RNA FASTQ download manifest from the sample table.

The input is the tab-separated 93-sample table used in this project. PRO-seq
GSM accessions are resolved to SRA runs through NCBI RunInfo, then each run is
resolved to the public ENA FASTQ files and checksums.
"""

from __future__ import annotations

import argparse
import csv
import io
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path


USER_AGENT = "PROseqManifestBuilder/1.0 (public-data-retrieval)"
GSM_RE = re.compile(r"GSM\d+")
GSE_RE = re.compile(r"GSE\d+")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-table", required=True, type=Path)
    parser.add_argument("--method", default="PRO-seq")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--timeout", type=int, default=60)
    return parser.parse_args()


def fetch_text(url: str, timeout: int, attempts: int = 5) -> str:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    last_error: Exception | None = None
    for attempt in range(1, attempts + 1):
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return response.read().decode("utf-8-sig")
        except Exception as exc:
            last_error = exc
            if attempt < attempts:
                time.sleep(min(2**attempt, 15))
    raise RuntimeError(f"Failed to retrieve {url}: {last_error}")


def read_assay_samples(path: Path, method: str) -> list[dict[str, str]]:
    samples: list[dict[str, str]] = []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("Method", "").strip() != method:
                continue
            labels = [item.strip() for item in row["Sample IDs"].split(";")]
            for label in labels:
                match = GSM_RE.search(label)
                if not match:
                    raise ValueError(f"No GSM accession in {method} label: {label}")
                samples.append(
                    {
                        "group_id": row["GroupID"].strip(),
                        "cell_types": row["Cell types"].strip(),
                        "gse": (
                            GSE_RE.search(label).group(0)
                            if GSE_RE.search(label)
                            else row["GSE"].strip()
                        ),
                        "reference": row["Reference"].strip(),
                        "sample_label": label,
                        "gsm": match.group(0),
                    }
                )
    if not samples:
        raise ValueError(f"No {method} samples were found")
    return samples


def gsm_runs(gsm: str, timeout: int) -> list[dict[str, str]]:
    url = (
        "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc="
        + urllib.parse.quote(gsm)
    )
    text = fetch_text(url, timeout)
    rows = list(csv.DictReader(io.StringIO(text)))
    rows = [row for row in rows if row.get("Run", "").startswith("SRR")]
    if not rows:
        raise ValueError(f"No SRR run was returned for {gsm}")
    return rows


def ena_fastqs(run: str, timeout: int) -> list[dict[str, str]]:
    fields = ",".join(
        [
            "run_accession",
            "library_layout",
            "fastq_ftp",
            "fastq_md5",
            "fastq_bytes",
            "scientific_name",
            "sample_accession",
            "experiment_accession",
        ]
    )
    url = (
        "https://www.ebi.ac.uk/ena/portal/api/filereport?"
        + urllib.parse.urlencode(
            {
                "accession": run,
                "result": "read_run",
                "fields": fields,
                "format": "tsv",
            }
        )
    )
    text = fetch_text(url, timeout)
    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    if len(rows) != 1 or not rows[0].get("fastq_ftp"):
        raise ValueError(f"ENA did not return FASTQ files for {run}")
    row = rows[0]
    urls = row["fastq_ftp"].split(";")
    md5s = row["fastq_md5"].split(";")
    sizes = row["fastq_bytes"].split(";")
    if not (len(urls) == len(md5s) == len(sizes)):
        raise ValueError(f"Inconsistent ENA FASTQ metadata for {run}")

    inferred_layout = row["library_layout"]
    if (
        len(urls) == 2
        and urls[0].endswith("_1.fastq.gz")
        and urls[1].endswith("_2.fastq.gz")
    ):
        inferred_layout = "PAIRED"

    result: list[dict[str, str]] = []
    for index, (ftp_path, md5, size) in enumerate(zip(urls, md5s, sizes), 1):
        url_value = ftp_path if "://" in ftp_path else "https://" + ftp_path
        filename = url_value.rsplit("/", 1)[-1]
        if inferred_layout == "PAIRED":
            mate = "1" if filename.endswith("_1.fastq.gz") else "2"
        else:
            mate = "0"
        result.append(
            {
                "run": run,
                "layout": inferred_layout,
                "reported_layout": row["library_layout"],
                "mate": mate,
                "filename": filename,
                "url": url_value,
                "md5": md5,
                "bytes": size,
                "scientific_name": row.get("scientific_name", ""),
                "sample_accession": row.get("sample_accession", ""),
                "experiment_accession": row.get("experiment_accession", ""),
            }
        )
    return result


def main() -> int:
    args = parse_args()
    samples = read_assay_samples(args.sample_table, args.method)
    output_rows: list[dict[str, str]] = []
    seen_files: set[str] = set()

    for number, sample in enumerate(samples, 1):
        print(f"[{number}/{len(samples)}] resolving {sample['gsm']}", file=sys.stderr)
        for run_row in gsm_runs(sample["gsm"], args.timeout):
            run = run_row["Run"]
            for fastq in ena_fastqs(run, args.timeout):
                if fastq["filename"] in seen_files:
                    continue
                seen_files.add(fastq["filename"])
                output_rows.append(
                    {
                        "method": args.method,
                        **sample,
                        **fastq,
                    }
                )
        time.sleep(0.12)

    fieldnames = [
        "method",
        "group_id",
        "cell_types",
        "gse",
        "reference",
        "sample_label",
        "gsm",
        "run",
        "layout",
        "reported_layout",
        "mate",
        "filename",
        "url",
        "md5",
        "bytes",
        "scientific_name",
        "sample_accession",
        "experiment_accession",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)

    runs = {row["run"] for row in output_rows}
    layouts = {row["run"]: row["layout"] for row in output_rows}
    paired = sum(value == "PAIRED" for value in layouts.values())
    total_bytes = sum(int(row["bytes"]) for row in output_rows)
    print(
        f"Wrote {len(output_rows)} FASTQ records for {len(runs)} runs "
        f"({paired} paired runs; {total_bytes / 1024**3:.2f} GiB)",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
