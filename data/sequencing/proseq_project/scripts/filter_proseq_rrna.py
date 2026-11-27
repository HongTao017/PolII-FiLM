#!/usr/bin/env python3
"""Remove human rRNA reads from cleaned PRO-seq FASTQ files."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import shlex
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path


MISSING = {"", "NA", "N/A", "."}


@dataclass(frozen=True)
class Job:
    gse: str
    group_id: str
    run: str
    layout: str
    input_r1: Path
    input_r2: Path | None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Filter cleaned PRO-seq FASTQs against a Bowtie2 rRNA index."
    )
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--rrna-index", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--log-dir", required=True, type=Path)
    parser.add_argument("--status", required=True, type=Path)
    parser.add_argument("--threads-per-job", type=int, default=2)
    parser.add_argument("--jobs", type=int, default=3)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def read_jobs(path: Path) -> list[Job]:
    jobs: list[Job] = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"gse", "group_id", "run", "layout", "output_r1", "output_r2"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Manifest lacks columns: {', '.join(sorted(missing))}")

        for row in reader:
            input_r1_text = row["output_r1"].strip()
            input_r2_text = row["output_r2"].strip()
            input_r2 = None if input_r2_text in MISSING else Path(input_r2_text)
            jobs.append(
                Job(
                    gse=row["gse"].strip(),
                    group_id=row["group_id"].strip(),
                    run=row["run"].strip(),
                    layout=row["layout"].strip().upper(),
                    input_r1=Path(input_r1_text),
                    input_r2=input_r2,
                )
            )
    return jobs


def require_inputs(job: Job) -> None:
    if not job.input_r1.is_file():
        raise FileNotFoundError(f"Missing R1: {job.input_r1}")
    if job.layout == "PAIRED":
        if job.input_r2 is None or not job.input_r2.is_file():
            raise FileNotFoundError(f"Missing R2 for paired Run {job.run}")
    elif job.layout != "SINGLE":
        raise ValueError(f"Unsupported layout for {job.run}: {job.layout}")


def output_paths(job: Job, output_dir: Path) -> tuple[Path, Path | None]:
    if job.layout == "PAIRED":
        return (
            output_dir / f"{job.run}_1.non_rRNA.fastq.gz",
            output_dir / f"{job.run}_2.non_rRNA.fastq.gz",
        )
    return output_dir / f"{job.run}.non_rRNA.fastq.gz", None


def temporary_fastq_path(path: Path) -> Path:
    """Return a temporary name that still ends in .gz for samtools fastq."""
    suffix = ".fastq.gz"
    if not path.name.endswith(suffix):
        raise ValueError(f"Expected a .fastq.gz output path: {path}")
    return path.with_name(path.name[: -len(suffix)] + ".tmp" + suffix)


def gzip_ok(path: Path) -> bool:
    result = subprocess.run(
        ["gzip", "-t", str(path)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
    )
    return result.returncode == 0


def count_fastq_records(path: Path) -> int:
    lines = 0
    with gzip.open(path, "rb") as handle:
        for _ in handle:
            lines += 1
    if lines % 4:
        raise ValueError(f"FASTQ line count is not divisible by four: {path}")
    return lines // 4


def parse_counts(log_path: Path, layout: str, output_r1: Path) -> tuple[int, int]:
    text = log_path.read_text(encoding="utf-8", errors="replace")
    total_match = re.search(r"(?m)^(\d+) reads; of these:", text)
    if not total_match:
        raise ValueError("Could not parse the Bowtie2 input read count")
    total = int(total_match.group(1))

    if layout == "SINGLE":
        retained_match = re.search(
            r"(?m)^\s*(\d+) \([^\n]+\) aligned 0 times\s*$", text
        )
        retained = int(retained_match.group(1)) if retained_match else count_fastq_records(output_r1)
    else:
        processed = re.findall(r"processed (\d+) reads", text)
        retained = int(processed[-1]) // 2 if processed else count_fastq_records(output_r1)
    return total, retained


def run_command(command: str, log_path: Path) -> None:

    result = subprocess.run(["bash", "-c", command])
    if result.returncode:
        raise subprocess.CalledProcessError(result.returncode, command)
    if not log_path.is_file():
        raise RuntimeError(f"Expected log was not created: {log_path}")


def filter_job(
    job: Job,
    index: Path,
    output_dir: Path,
    log_dir: Path,
    threads: int,
    overwrite: bool,
) -> dict[str, str | int | float]:
    output_r1, output_r2 = output_paths(job, output_dir)
    log_path = log_dir / f"{job.run}.rRNA.bowtie2.log"
    require_inputs(job)

    outputs = [output_r1] + ([output_r2] if output_r2 else [])
    already_done = all(path and path.is_file() and gzip_ok(path) for path in outputs)

    if overwrite or not already_done:
        for path in outputs:
            if path:
                path.unlink(missing_ok=True)

        q = shlex.quote
        if job.layout == "SINGLE":
            temp_r1 = temporary_fastq_path(output_r1)
            temp_r1.unlink(missing_ok=True)
            command = (
                "set -o pipefail; "
                f"bowtie2 --very-sensitive-local --threads {threads} "
                f"-x {q(str(index))} -U {q(str(job.input_r1))} "
                f"--un-gz {q(str(temp_r1))} -S /dev/null "
                f"2> {q(str(log_path))}"
            )
            run_command(command, log_path)
            temp_r1.replace(output_r1)
        else:
            assert job.input_r2 is not None and output_r2 is not None
            temp_r1 = temporary_fastq_path(output_r1)
            temp_r2 = temporary_fastq_path(output_r2)
            temp_r1.unlink(missing_ok=True)
            temp_r2.unlink(missing_ok=True)
            command = (
                "set -o pipefail; "
                f"bowtie2 --very-sensitive-local --threads {threads} "
                f"-x {q(str(index))} -1 {q(str(job.input_r1))} -2 {q(str(job.input_r2))} "
                f"2> {q(str(log_path))} | "
                f"samtools view -@ {threads} -b -f 12 -F 2304 - | "
                f"samtools sort -n -@ {threads} -m 768M - | "
                f"samtools fastq -@ {threads} -n "
                f"-1 {q(str(temp_r1))} -2 {q(str(temp_r2))} "
                f"-0 /dev/null -s /dev/null - 2>> {q(str(log_path))}"
            )
            run_command(command, log_path)
            temp_r1.replace(output_r1)
            temp_r2.replace(output_r2)

    for path in outputs:
        if path is None or not path.is_file() or not gzip_ok(path):
            raise RuntimeError(f"Invalid filtered FASTQ: {path}")

    total, retained = parse_counts(log_path, job.layout, output_r1)
    if job.layout == "PAIRED" and output_r2 is not None:
        retained_r2 = count_fastq_records(output_r2)
        if retained_r2 != retained:
            raise ValueError(
                f"Paired output count mismatch for {job.run}: R1={retained}, R2={retained_r2}"
            )
    removed = total - retained
    percent = 100.0 * removed / total if total else 0.0
    return {
        "gse": job.gse,
        "group_id": job.group_id,
        "run": job.run,
        "layout": job.layout,
        "input_r1": str(job.input_r1),
        "input_r2": str(job.input_r2) if job.input_r2 else "NA",
        "output_r1": str(output_r1),
        "output_r2": str(output_r2) if output_r2 else "NA",
        "input_fragments": total,
        "retained_fragments": retained,
        "removed_fragments": removed,
        "rrna_percent": f"{percent:.2f}",
        "status": "PASS",
        "note": "existing_valid_outputs" if already_done and not overwrite else "filtered_and_verified",
    }


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.log_dir.mkdir(parents=True, exist_ok=True)
    args.status.parent.mkdir(parents=True, exist_ok=True)

    jobs = read_jobs(args.manifest)
    if not jobs:
        raise ValueError("No runs found in the cleaning manifest")

    results: list[dict[str, str | int | float]] = []
    failures: list[str] = []
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {
            pool.submit(
                filter_job,
                job,
                args.rrna_index,
                args.output_dir,
                args.log_dir,
                args.threads_per_job,
                args.overwrite,
            ): job
            for job in jobs
        }
        completed = 0
        for future in as_completed(futures):
            job = futures[future]
            completed += 1
            try:
                result = future.result()
                results.append(result)
                print(
                    f"[{completed}/{len(jobs)}] {job.run}: PASS, "
                    f"rRNA={result['rrna_percent']}%",
                    flush=True,
                )
            except Exception as exc:
                failures.append(f"{job.run}: {exc}")
                print(f"[{completed}/{len(jobs)}] {job.run}: FAILED: {exc}", flush=True)

    fields = [
        "gse",
        "group_id",
        "run",
        "layout",
        "input_r1",
        "input_r2",
        "output_r1",
        "output_r2",
        "input_fragments",
        "retained_fragments",
        "removed_fragments",
        "rrna_percent",
        "status",
        "note",
    ]
    with tempfile.NamedTemporaryFile(
        "w", delete=False, dir=args.status.parent, newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(sorted(results, key=lambda row: str(row["run"])))
        temporary_status = Path(handle.name)
    os.replace(temporary_status, args.status)

    print(f"Completed: {len(results)} PASS, {len(failures)} FAILED", flush=True)
    print(f"Status: {args.status}", flush=True)
    if failures:
        raise SystemExit("\n".join(failures))


if __name__ == "__main__":
    main()
