#!/usr/bin/env python3
"""Run one FastQC-informed cutadapt command per PRO-seq run."""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


SUMMARY_PATTERNS = {
    "total_reads": re.compile(r"Total read pairs processed:\s+([\d,]+)|Total reads processed:\s+([\d,]+)"),
    "written_reads": re.compile(r"Pairs written \(passing filters\):\s+([\d,]+)|Reads written \(passing filters\):\s+([\d,]+)"),
}


def read_tsv(path: Path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def gzip_ok(path: Path) -> bool:
    if not path.is_file() or path.stat().st_size == 0:
        return False
    try:
        with gzip.open(path, "rb") as handle:
            while handle.read(1024 * 1024):
                pass
        return True
    except (OSError, EOFError):
        return False


def parse_count(pattern, text):
    match = pattern.search(text)
    if not match:
        return "NA"
    value = next(item for item in match.groups() if item is not None)
    return value.replace(",", "")


def clean_one(row, log_dir: Path, threads: int):
    run = row["run"]
    output_r1 = Path(row["output_r1"])
    output_r2 = None if row["output_r2"] == "NA" else Path(row["output_r2"])
    output_r1.parent.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / f"{run}.cutadapt.log"

    expected_outputs = [output_r1] + ([output_r2] if output_r2 else [])
    if all(gzip_ok(path) for path in expected_outputs):
        return {
            "run": run,
            "gse": row["gse"],
            "layout": row["layout"],
            "status": "SKIP_VALID_EXISTING",
            "seconds": "0",
            "total_reads": "NA",
            "written_reads": "NA",
            "retention_percent": "NA",
            "log": str(log_path),
        }

    for path in expected_outputs:
        if path.exists():
            path.unlink()

    command = [
        "cutadapt",
        "-j", str(threads),
        "-e", "0.10",
        "-O", "5",
        "--trim-n",
        "--max-n", row["max_n_fraction"],
        "-m", row["minimum_length"],
    ]
    if row["quality_mode"] == "nextseq_3prime":
        command.extend(["--nextseq-trim", row["quality_cutoff"]])
    else:
        command.extend(["-q", f"0,{row['quality_cutoff']}"])

    command.extend(["-a", row["adapter_r1"]])
    if row["layout"] == "PAIRED":
        command.extend(
            [
                "-A", row["adapter_r2"],
                "--pair-filter", "any",
                "-o", row["output_r1"],
                "-p", row["output_r2"],
                row["input_r1"],
                row["input_r2"],
            ]
        )
    else:
        command.extend(["-o", row["output_r1"], row["input_r1"]])

    started = time.time()
    with log_path.open("w", encoding="utf-8") as log_handle:
        log_handle.write("COMMAND\t" + " ".join(command) + "\n")
        completed = subprocess.run(
            command,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
    seconds = int(time.time() - started)
    text = log_path.read_text(encoding="utf-8", errors="replace")
    total = parse_count(SUMMARY_PATTERNS["total_reads"], text)
    written = parse_count(SUMMARY_PATTERNS["written_reads"], text)
    retention = "NA"
    if total != "NA" and written != "NA" and int(total):
        retention = f"{100 * int(written) / int(total):.2f}"

    outputs_valid = completed.returncode == 0 and all(
        path.is_file() and path.stat().st_size > 0 for path in expected_outputs
    )
    return {
        "run": run,
        "gse": row["gse"],
        "layout": row["layout"],
        "status": "PASS" if outputs_valid else f"FAIL_EXIT_{completed.returncode}",
        "seconds": str(seconds),
        "total_reads": total,
        "written_reads": written,
        "retention_percent": retention,
        "log": str(log_path),
    }


def write_status(path: Path, records):
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "run", "gse", "layout", "status", "seconds", "total_reads",
        "written_reads", "retention_percent", "log",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(sorted(records, key=lambda item: item["run"]))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--log-dir", required=True, type=Path)
    parser.add_argument("--status", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--threads-per-run", type=int, default=2)
    args = parser.parse_args()

    rows = read_tsv(args.manifest)
    results = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(clean_one, row, args.log_dir, args.threads_per_run): row
            for row in rows
        }
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            write_status(args.status, results)
            print(
                f"{result['status']} {result['run']} "
                f"retention={result['retention_percent']}% completed={len(results)}/{len(rows)}",
                flush=True,
            )

    failed = [item for item in results if not item["status"].startswith(("PASS", "SKIP"))]
    print(f"FINISH runs={len(results)} failed={len(failed)} status={args.status}")
    raise SystemExit(1 if failed else 0)


if __name__ == "__main__":
    main()
