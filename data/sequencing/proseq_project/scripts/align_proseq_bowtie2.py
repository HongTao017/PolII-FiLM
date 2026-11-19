#!/usr/bin/env python3
"""Align rRNA-filtered PRO-seq FASTQs to GRCh38 with Bowtie2."""

from __future__ import annotations

import argparse
import csv
import re
import shlex
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


STANDARD_CHROMS = [f"chr{number}" for number in range(1, 23)] + ["chrX", "chrY"]
ALIGNMENT_RATE_RE = re.compile(r"([0-9]+(?:\.[0-9]+)?)% overall alignment rate")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rrna-status", required=True, type=Path)
    parser.add_argument("--bowtie2-index", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--log-dir", required=True, type=Path)
    parser.add_argument("--status", required=True, type=Path)
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--threads-per-job", type=int, default=8)
    parser.add_argument("--mapq", type=int, default=20)
    return parser.parse_args()


def quote(value: str | Path) -> str:
    return shlex.quote(str(value))


def read_jobs(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {
        "gse",
        "group_id",
        "run",
        "layout",
        "output_r1",
        "output_r2",
        "retained_fragments",
        "status",
    }
    if not rows or not required.issubset(rows[0]):
        raise ValueError("The rRNA status table is missing required columns")
    failed = [row["run"] for row in rows if row["status"] != "PASS"]
    if failed:
        raise ValueError("rRNA filtering did not pass for: " + ", ".join(failed))
    if len({row["run"] for row in rows}) != len(rows):
        raise ValueError("Duplicate Run accessions were found")
    return rows


def check_index(prefix: Path) -> None:
    small = [Path(f"{prefix}.{suffix}.bt2") for suffix in (1, 2, 3, 4)] + [
        Path(f"{prefix}.rev.{suffix}.bt2") for suffix in (1, 2)
    ]
    large = [Path(f"{prefix}.{suffix}.bt2l") for suffix in (1, 2, 3, 4)] + [
        Path(f"{prefix}.rev.{suffix}.bt2l") for suffix in (1, 2)
    ]
    if not all(path.is_file() for path in small) and not all(path.is_file() for path in large):
        raise FileNotFoundError(f"Bowtie2 index is incomplete: {prefix}")


def count_records(bam: Path, threads: int) -> int:
    result = subprocess.run(
        ["samtools", "view", "-@", str(threads), "-c", str(bam)],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    )
    return int(result.stdout.strip())


def bam_standard_chroms(bam: Path) -> list[str]:
    result = subprocess.run(
        ["samtools", "view", "-H", str(bam)],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    )
    present = {
        token.split(":", 1)[1]
        for line in result.stdout.splitlines()
        if line.startswith("@SQ\t")
        for token in line.split("\t")
        if token.startswith("SN:")
    }
    chroms = [chrom for chrom in STANDARD_CHROMS if chrom in present]
    if len(chroms) != 24:
        missing = sorted(set(STANDARD_CHROMS) - set(chroms))
        raise ValueError("BAM lacks standard chromosomes: " + ", ".join(missing))
    return chroms


def parse_alignment_rate(log: Path) -> str:
    match = ALIGNMENT_RATE_RE.search(log.read_text(encoding="utf-8", errors="replace"))
    return match.group(1) if match else "NA"


def quickcheck(bam: Path) -> None:
    subprocess.run(["samtools", "quickcheck", "-v", str(bam)], check=True)


def align_one(
    row: dict[str, str],
    args: argparse.Namespace,
) -> dict[str, str]:
    run = row["run"]
    layout = row["layout"].upper()
    output = args.output_dir / f"{run}.GRCh38.standard.primary.q{args.mapq}.bam"
    index = Path(f"{output}.bai")
    log = args.log_dir / f"{run}.bowtie2.log"
    temp_bam = args.output_dir / f".{run}.GRCh38.full.sorted.tmp.bam"
    temp_index = Path(f"{temp_bam}.bai")
    started = time.time()
    result = {
        "run": run,
        "gse": row["gse"],
        "group_id": row["group_id"],
        "layout": layout,
        "input_fragments": row["retained_fragments"],
        "bowtie2_alignment_percent": "NA",
        "q20_records": "NA",
        "q20_fragments": "NA",
        "q20_retention_percent": "NA",
        "output_bam": str(output),
        "output_bai": str(index),
        "seconds": "NA",
        "status": "FAILED",
        "note": "",
    }
    try:
        r1 = Path(row["output_r1"])
        r2 = Path(row["output_r2"]) if row["output_r2"] != "NA" else None
        if not r1.is_file():
            raise FileNotFoundError(r1)
        if layout == "PAIRED" and (r2 is None or not r2.is_file()):
            raise FileNotFoundError(r2 or "missing R2")
        if layout not in {"SINGLE", "PAIRED"}:
            raise ValueError(f"Unsupported layout: {layout}")

        if output.is_file() and index.is_file():
            quickcheck(output)
        else:
            output.unlink(missing_ok=True)
            index.unlink(missing_ok=True)
            temp_bam.unlink(missing_ok=True)
            temp_index.unlink(missing_ok=True)
            input_args = (
                f"-U {quote(r1)}"
                if layout == "SINGLE"
                else f"-1 {quote(r1)} -2 {quote(r2)}"
            )
            align_command = (
                f"bowtie2 --very-sensitive --threads {args.threads_per_job} "
                f"-x {quote(args.bowtie2_index)} {input_args} 2> {quote(log)} | "
                f"samtools view -@ {args.threads_per_job} -b - | "
                f"samtools sort -@ {args.threads_per_job} -o {quote(temp_bam)} -"
            )
            subprocess.run(
                ["bash", "-c", f"set -o pipefail; {align_command}"], check=True
            )
            subprocess.run(
                ["samtools", "index", "-@", str(args.threads_per_job), str(temp_bam)],
                check=True,
            )
            chroms = bam_standard_chroms(temp_bam)
            filter_command = [
                "samtools",
                "view",
                "-@",
                str(args.threads_per_job),
                "-b",
                "-q",
                str(args.mapq),
                "-F",
                "2308",
            ]
            if layout == "PAIRED":
                filter_command.extend(["-f", "2"])
            filter_command.extend(["-o", str(output), str(temp_bam), *chroms])
            subprocess.run(filter_command, check=True)
            subprocess.run(
                ["samtools", "index", "-@", str(args.threads_per_job), str(output)],
                check=True,
            )
            quickcheck(output)
            temp_bam.unlink(missing_ok=True)
            temp_index.unlink(missing_ok=True)

        records = count_records(output, args.threads_per_job)
        fragments = records // 2 if layout == "PAIRED" else records
        if layout == "PAIRED" and records % 2:
            raise ValueError(f"Proper-pair BAM has an odd record count: {records}")
        input_fragments = int(row["retained_fragments"])
        retention = 100 * fragments / input_fragments if input_fragments else 0.0
        result.update(
            {
                "bowtie2_alignment_percent": parse_alignment_rate(log)
                if log.is_file()
                else "NA",
                "q20_records": str(records),
                "q20_fragments": str(fragments),
                "q20_retention_percent": f"{retention:.2f}",
                "seconds": str(round(time.time() - started)),
                "status": "PASS",
                "note": (
                    "standard_chromosomes_primary_mapq20_proper_pairs"
                    if layout == "PAIRED"
                    else "standard_chromosomes_primary_mapq20"
                ),
            }
        )
    except Exception as exc:
        temp_bam.unlink(missing_ok=True)
        temp_index.unlink(missing_ok=True)
        result["seconds"] = str(round(time.time() - started))
        result["note"] = f"{type(exc).__name__}: {exc}"
    return result


def write_status(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda row: row["run"]))


def main() -> int:
    args = parse_args()
    for tool in ("bowtie2", "samtools"):
        if not shutil.which(tool):
            raise RuntimeError(f"Required program is unavailable: {tool}")
    check_index(args.bowtie2_index)
    jobs = read_jobs(args.rrna_status)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.log_dir.mkdir(parents=True, exist_ok=True)

    results: list[dict[str, str]] = []
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {executor.submit(align_one, row, args): row["run"] for row in jobs}
        for number, future in enumerate(as_completed(futures), 1):
            result = future.result()
            results.append(result)
            write_status(args.status, results)
            print(
                f"[{number}/{len(jobs)}] {result['run']}: {result['status']}, "
                f"alignment={result['bowtie2_alignment_percent']}%, "
                f"q20_fragments={result['q20_fragments']}",
                flush=True,
            )

    passed = sum(row["status"] == "PASS" for row in results)
    print(f"Completed: {passed} PASS, {len(results) - passed} FAILED")
    return 0 if passed == len(results) else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {type(exc).__name__}: {exc}", file=sys.stderr)
        raise SystemExit(1)
