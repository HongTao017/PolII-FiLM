#!/usr/bin/env python3
"""Align paired mNET-seq reads and retain high-confidence human pairs."""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


STANDARD_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--star-index", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--log-dir", required=True, type=Path)
    parser.add_argument("--status", required=True, type=Path)
    parser.add_argument("--jobs", type=int, default=3)
    parser.add_argument("--threads-per-job", type=int, default=8)
    parser.add_argument("--mapq", type=int, default=7)
    return parser.parse_args()


def run_command(command: list[str], stdout: Path | None = None) -> None:
    if stdout is None:
        subprocess.run(command, check=True)
        return
    with stdout.open("w", encoding="utf-8") as handle:
        subprocess.run(command, check=True, stdout=handle, stderr=subprocess.STDOUT)


def parse_star_log(path: Path) -> tuple[int, float]:
    text = path.read_text(encoding="utf-8", errors="replace")
    input_match = re.search(r"Number of input reads\s*\|\s*(\d+)", text)
    unique_match = re.search(r"Uniquely mapped reads %\s*\|\s*([0-9.]+)%", text)
    if not input_match or not unique_match:
        raise ValueError(f"Could not parse STAR summary: {path}")
    return int(input_match.group(1)), float(unique_match.group(1))


def count_records(bam: Path, threads: int) -> int:
    result = subprocess.run(
        ["samtools", "view", "-@", str(threads), "-c", str(bam)],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    )
    return int(result.stdout.strip())


def align_one(
    row: dict[str, str],
    args: argparse.Namespace,
) -> dict[str, str]:
    run = row["run"]
    final_bam = args.output_dir / f"{run}.GRCh38.standard.unique.proper.bam"
    final_bai = Path(f"{final_bam}.bai")
    run_log_dir = args.log_dir / run
    star_prefix = run_log_dir / "STAR."
    star_bam = run_log_dir / "STAR.Aligned.sortedByCoord.out.bam"
    command_log = run_log_dir / "command.txt"

    run_log_dir.mkdir(parents=True, exist_ok=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    result = {
        **row,
        "star_input_pairs": "NA",
        "star_unique_percent": "NA",
        "retained_pairs": "NA",
        "output_bam": str(final_bam),
        "status": "FAILED",
        "note": "",
    }

    try:
        for key in ("input_r1", "input_r2"):
            path = Path(row[key])
            if not path.is_file() or path.stat().st_size == 0:
                raise FileNotFoundError(path)

        if final_bam.is_file() and final_bai.is_file():
            subprocess.run(["samtools", "quickcheck", str(final_bam)], check=True)
        else:
            star_command = [
                "STAR",
                "--runThreadN",
                str(args.threads_per_job),
                "--genomeDir",
                str(args.star_index),
                "--readFilesIn",
                row["input_r1"],
                row["input_r2"],
                "--readFilesCommand",
                "zcat",
                "--outFileNamePrefix",
                str(star_prefix),
                "--outSAMtype",
                "BAM",
                "SortedByCoordinate",
                "--outSAMattributes",
                "NH",
                "HI",
                "AS",
                "nM",
                "MD",
                "--outFilterMultimapNmax",
                "1",
                "--outFilterMismatchNoverLmax",
                "0.10",
                "--alignEndsType",
                "EndToEnd",
                "--alignIntronMin",
                "11",
                "--alignIntronMax",
                "1000000",
                "--outSAMunmapped",
                "None",
                "--outSAMstrandField",
                "intronMotif",
                "--limitBAMsortRAM",
                "12000000000",
            ]
            command_log.write_text(" ".join(star_command) + "\n", encoding="utf-8")
            run_command(star_command, run_log_dir / "STAR.stdout.log")
            run_command(
                [
                    "samtools",
                    "index",
                    "-@",
                    str(args.threads_per_job),
                    str(star_bam),
                ]
            )

            filtered_tmp = Path(f"{final_bam}.tmp")
            filter_command = [
                "samtools",
                "view",
                "-@",
                str(args.threads_per_job),
                "-b",
                "-q",
                str(args.mapq),
                "-f",
                "2",
                "-F",
                "2304",
                "-o",
                str(filtered_tmp),
                str(star_bam),
                *STANDARD_CHROMS,
            ]
            with command_log.open("a", encoding="utf-8") as handle:
                handle.write(" ".join(filter_command) + "\n")
            run_command(filter_command)
            filtered_tmp.replace(final_bam)
            run_command(["samtools", "index", "-@", str(args.threads_per_job), str(final_bam)])
            subprocess.run(["samtools", "quickcheck", str(final_bam)], check=True)
            star_bam.unlink(missing_ok=True)
            Path(f"{star_bam}.bai").unlink(missing_ok=True)

        input_pairs, unique_percent = parse_star_log(run_log_dir / "STAR.Log.final.out")
        records = count_records(final_bam, args.threads_per_job)
        if records % 2:
            raise ValueError(f"Odd number of records in paired BAM: {records}")
        result.update(
            {
                "star_input_pairs": str(input_pairs),
                "star_unique_percent": f"{unique_percent:.2f}",
                "retained_pairs": str(records // 2),
                "status": "PASS",
                "note": "unique_proper_standard_chromosomes",
            }
        )
    except Exception as exc:  # Keep the batch running and report the failing Run.
        result["note"] = f"{type(exc).__name__}: {exc}"
    return result


def write_status(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "run",
        "gse",
        "gsm",
        "group_id",
        "sample_label",
        "input_r1",
        "input_r2",
        "rrna_percent",
        "star_input_pairs",
        "star_unique_percent",
        "retained_pairs",
        "output_bam",
        "status",
        "note",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda row: (row["gse"], row["gsm"], row["run"])))


def main() -> None:
    args = parse_args()
    for tool in ("STAR", "samtools"):
        if not shutil.which(tool):
            raise RuntimeError(f"Required tool not found: {tool}")
    if not (args.star_index / "Genome").is_file():
        raise FileNotFoundError(f"Incomplete STAR index: {args.star_index}")

    with args.manifest.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    results: list[dict[str, str]] = []
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {executor.submit(align_one, row, args): row["run"] for row in rows}
        for index, future in enumerate(as_completed(futures), 1):
            result = future.result()
            results.append(result)
            write_status(args.status, results)
            print(
                f"[{index}/{len(rows)}] {result['run']}: {result['status']}, "
                f"unique={result['star_unique_percent']}%, pairs={result['retained_pairs']}",
                flush=True,
            )

    passed = sum(row["status"] == "PASS" for row in results)
    print(f"Completed: {passed} PASS, {len(results) - passed} FAILED")
    print(f"Status: {args.status}")
    if passed != len(results):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
