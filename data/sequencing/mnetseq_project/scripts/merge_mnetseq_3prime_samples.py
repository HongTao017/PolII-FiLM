#!/usr/bin/env python3
"""Merge technical mNET-seq Runs by GEO sample (GSE + GSM)."""

from __future__ import annotations

import argparse
import csv
import os
import shlex
import shutil
import subprocess
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-status", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--temp-dir", required=True, type=Path)
    parser.add_argument("--status", required=True, type=Path)
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--threads-per-job", type=int, default=8)
    return parser.parse_args()


def quote(value: str | Path) -> str:
    return shlex.quote(str(value))


def write_status(path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames = [
        "sample_id",
        "gse",
        "gsm",
        "sample_label",
        "technical_run_count",
        "runs",
        "total_3prime_ends",
        "output_bed_gz",
        "status",
        "note",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda row: row["sample_id"]))


def merge_one(rows: list[dict[str, str]], args: argparse.Namespace) -> dict[str, str]:
    first = rows[0]
    sample_id = f"{first['gse']}_{first['gsm']}"
    output = args.output_dir / f"{sample_id}.mNET.3prime.bed.gz"
    output_index = Path(f"{output}.tbi")
    temp_output = args.output_dir / f"{sample_id}.mNET.3prime.tmp.bed.gz"
    sample_temp = args.temp_dir / sample_id
    count_file = sample_temp / "total_count.txt"
    result = {
        "sample_id": sample_id,
        "gse": first["gse"],
        "gsm": first["gsm"],
        "sample_label": first["sample_label"],
        "technical_run_count": str(len(rows)),
        "runs": ";".join(sorted(row["run"] for row in rows)),
        "total_3prime_ends": "NA",
        "output_bed_gz": str(output),
        "status": "FAILED",
        "note": "",
    }
    try:
        inputs = [Path(row["output_bed_gz"]) for row in rows]
        for row, path in zip(rows, inputs):
            if row["status"] != "PASS" or not path.is_file():
                raise ValueError(f"Run is unavailable: {row['run']}")
        sample_temp.mkdir(parents=True, exist_ok=True)
        args.output_dir.mkdir(parents=True, exist_ok=True)
        count_path = str(count_file).replace("\\", "\\\\").replace('"', '\\"')
        if not (output.is_file() and output_index.is_file() and count_file.is_file()):
            expected_total = sum(int(row["retained_3prime_ends"]) for row in rows)
            if len(inputs) == 1:
                output.unlink(missing_ok=True)
                output_index.unlink(missing_ok=True)
                os.link(inputs[0], output)
                os.link(Path(f"{inputs[0]}.tbi"), output_index)
                count_file.write_text(f"{expected_total}\n", encoding="utf-8")
            else:
                input_list = " ".join(quote(path) for path in inputs)
                final_awk = (
                    "BEGIN{FS=OFS=\"\\t\";sum=0} "
                    "{sum+=$5;print $1,$2,$3,\".\",$5,$4} "
                    f'END{{print sum > "{count_path}"}}'
                )
                command = " | ".join(
                    [
                        f"zcat {input_list}",
                        (
                            f"LC_ALL=C sort --parallel={args.threads_per_job} -S 8G "
                            f"-T {quote(sample_temp)} -k1,1 -k2,2n -k3,3n -k6,6"
                        ),
                        "bedtools groupby -g 1,2,3,6 -c 5 -o sum",
                        f"awk {quote(final_awk)}",
                        f"bgzip -@ {args.threads_per_job} -c > {quote(temp_output)}",
                    ]
                )
                subprocess.run(["bash", "-c", f"set -o pipefail; {command}"], check=True)
                temp_output.replace(output)
                subprocess.run(["tabix", "-f", "-p", "bed", str(output)], check=True)
        result.update(
            {
                "total_3prime_ends": count_file.read_text(encoding="utf-8").strip(),
                "status": "PASS",
                "note": "technical_runs_summed_by_coordinate_and_strand",
            }
        )
    except Exception as exc:
        result["note"] = f"{type(exc).__name__}: {exc}"
    return result


def main() -> None:
    args = parse_args()
    for tool in ("bedtools", "bgzip", "tabix"):
        if not shutil.which(tool):
            raise RuntimeError(f"Required tool not found: {tool}")
    groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    with args.run_status.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            groups[(row["gse"], row["gsm"])].append(row)

    results: list[dict[str, str]] = []
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = [executor.submit(merge_one, rows, args) for rows in groups.values()]
        for index, future in enumerate(as_completed(futures), 1):
            result = future.result()
            results.append(result)
            write_status(args.status, results)
            print(
                f"[{index}/{len(groups)}] {result['sample_id']}: {result['status']}, "
                f"runs={result['technical_run_count']}, ends={result['total_3prime_ends']}",
                flush=True,
            )
    passed = sum(row["status"] == "PASS" for row in results)
    print(f"Completed: {passed} PASS, {len(results) - passed} FAILED")
    if passed != len(results):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
