#!/usr/bin/env python3
"""Clean mNET-seq FASTQs without breaking paired-end relationships."""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


TOTAL = re.compile(r"Total read pairs processed:\s+([\d,]+)|Total reads processed:\s+([\d,]+)")
WRITTEN = re.compile(r"Pairs written \(passing filters\):\s+([\d,]+)|Reads written \(passing filters\):\s+([\d,]+)")


def read_tsv(path: Path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def gzip_ok(path: Path) -> bool:
    if not path.is_file() or path.stat().st_size == 0:
        return False
    result = subprocess.run(
        ["gzip", "-t", str(path)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
    )
    return result.returncode == 0


def parsed_count(pattern: re.Pattern[str], text: str) -> str:
    match = pattern.search(text)
    if not match:
        return "NA"
    return next(value for value in match.groups() if value).replace(",", "")


def clean_one(row: dict[str, str], log_dir: Path, threads: int) -> dict[str, str]:
    run = row["run"]
    outputs = [Path(row["output_r1"])]
    if row["output_r2"] != "NA":
        outputs.append(Path(row["output_r2"]))
    for path in outputs:
        path.parent.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    log = log_dir / f"{run}.cutadapt.log"

    if all(gzip_ok(path) for path in outputs):
        return {"run": run, "gse": row["gse"], "layout": row["layout"],
                "status": "SKIP_VALID_EXISTING", "seconds": "0", "total_fragments": "NA",
                "written_fragments": "NA", "retention_percent": "NA", "log": str(log)}

    for path in outputs:
        path.unlink(missing_ok=True)

    command = ["cutadapt", "-j", str(threads), "-e", row["error_rate"],
               "-O", row["adapter_overlap"], "-n", "1", "--trim-n",
               "--max-n", row["max_n_fraction"], "-m", row["minimum_length"]]
    if row["match_read_wildcards"] == "yes":
        command.append("--match-read-wildcards")
    if row["quality_mode"] == "quality_3prime":
        command.extend(["-q", f"0,{row['quality_cutoff']}"])
    command.extend(["-a", row["adapter_r1"]])

    if row["layout"] == "PAIRED":
        command.extend(["-A", row["adapter_r2"], "--pair-filter", "any",
                        "-o", row["output_r1"], "-p", row["output_r2"],
                        row["input_r1"], row["input_r2"]])
    else:
        command.extend(["-o", row["output_r1"], row["input_r1"]])

    started = time.time()
    with log.open("w", encoding="utf-8") as handle:
        handle.write("COMMAND\t" + " ".join(command) + "\n")
        completed = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT)
    elapsed = int(time.time() - started)
    text = log.read_text(encoding="utf-8", errors="replace")
    total, written = parsed_count(TOTAL, text), parsed_count(WRITTEN, text)
    retention = "NA"
    if total != "NA" and written != "NA" and int(total):
        retention = f"{100 * int(written) / int(total):.2f}"
    valid = completed.returncode == 0 and all(path.is_file() and path.stat().st_size for path in outputs)
    return {"run": run, "gse": row["gse"], "layout": row["layout"],
            "status": "PASS" if valid else f"FAIL_EXIT_{completed.returncode}",
            "seconds": str(elapsed), "total_fragments": total,
            "written_fragments": written, "retention_percent": retention, "log": str(log)}


def write_status(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["run", "gse", "layout", "status", "seconds", "total_fragments",
              "written_fragments", "retention_percent", "log"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda item: item["run"]))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--log-dir", required=True, type=Path)
    parser.add_argument("--status", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--threads-per-run", type=int, default=2)
    args = parser.parse_args()
    rows, results = read_tsv(args.manifest), []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(clean_one, row, args.log_dir, args.threads_per_run): row for row in rows}
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            write_status(args.status, results)
            print(f"{result['status']} {result['run']} retention={result['retention_percent']}% "
                  f"completed={len(results)}/{len(rows)}", flush=True)
    failed = [row for row in results if not row["status"].startswith(("PASS", "SKIP"))]
    print(f"FINISH runs={len(results)} failed={len(failed)} status={args.status}")
    raise SystemExit(1 if failed else 0)


if __name__ == "__main__":
    main()
