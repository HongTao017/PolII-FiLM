#!/usr/bin/env python3
"""Extract strand-aware mNET-seq RNA 3' ends and remove internal priming."""

from __future__ import annotations

import argparse
import csv
import shlex
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment-status", required=True, type=Path)
    parser.add_argument("--cleaning-manifest", required=True, type=Path)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--temp-dir", required=True, type=Path)
    parser.add_argument("--log-dir", required=True, type=Path)
    parser.add_argument("--status", required=True, type=Path)
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--threads-per-job", type=int, default=8)
    return parser.parse_args()


def quote(value: str | Path) -> str:
    return shlex.quote(str(value))


def read_adapter_motifs(path: Path) -> dict[str, str]:
    motifs: dict[str, str] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            adapter = row["adapter_r1"].upper()
            if len(adapter) < 3:
                raise ValueError(f"R1 adapter is too short for {row['run']}: {adapter}")
            motifs[row["run"]] = adapter[:3]
    return motifs


def write_status(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "run",
        "gse",
        "gsm",
        "sample_label",
        "input_pairs",
        "adapter_motif",
        "retained_3prime_ends",
        "internal_priming_removed",
        "internal_priming_percent",
        "output_bed_gz",
        "status",
        "note",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda row: (row["gse"], row["gsm"], row["run"])))


def process_one(
    row: dict[str, str], motif: str, args: argparse.Namespace
) -> dict[str, str]:
    run = row["run"]
    bam = Path(row["output_bam"])
    output = args.output_dir / f"{run}.mNET.3prime.bed.gz"
    output_index = Path(f"{output}.tbi")
    temp_output = args.output_dir / f"{run}.mNET.3prime.tmp.bed.gz"
    run_temp = args.temp_dir / run
    run_log = args.log_dir / f"{run}.extract_3prime.log"
    count_file = run_temp / "retained_count.txt"
    input_pairs = int(row["retained_pairs"])

    result = {
        "run": run,
        "gse": row["gse"],
        "gsm": row["gsm"],
        "sample_label": row["sample_label"],
        "input_pairs": str(input_pairs),
        "adapter_motif": motif,
        "retained_3prime_ends": "NA",
        "internal_priming_removed": "NA",
        "internal_priming_percent": "NA",
        "output_bed_gz": str(output),
        "status": "FAILED",
        "note": "",
    }

    try:
        if row["status"] != "PASS":
            raise ValueError("Alignment status is not PASS")
        subprocess.run(["samtools", "quickcheck", str(bam)], check=True)
        run_temp.mkdir(parents=True, exist_ok=True)
        args.output_dir.mkdir(parents=True, exist_ok=True)
        args.log_dir.mkdir(parents=True, exist_ok=True)

        if not (output.is_file() and output_index.is_file() and count_file.is_file()):
            motif_length = len(motif)
            plus_r2_transform = (
                "BEGIN{OFS=\"\\t\"} "
                f"{{site=$2;tx=\"-\";s=site-{motif_length};e=site;"
                "if(s>=0)print $1,s,e,$1\"|\"site\"|\"tx\"|\"$3,0,tx}"
            )
            minus_r2_transform = (
                "BEGIN{OFS=\"\\t\"} "
                f"{{site=$2;tx=\"+\";s=site+1;e=s+{motif_length};"
                "print $1,s,e,$1\"|\"site\"|\"tx\"|\"$3,0,tx}"
            )
            fasta_filter = (
                "BEGIN{FS=OFS=\"\\t\"} "
                f"{{split($1,b,\"[(]\");id=b[1];"
                f"if(toupper($2)!=\"{motif}\"){{split(id,a,\"[|]\");"
                "print a[1],a[2],a[2]+1,\".\",a[4],a[3]}}"
            )
            count_path = str(count_file).replace("\\", "\\\\").replace('"', '\\"')
            final_format = (
                "BEGIN{FS=OFS=\"\\t\";sum=0} "
                "{sum+=$5;print $1,$2,$3,$4,$5,$6} "
                f'END{{print sum > "{count_path}"}}'
            )
            r2_5prime_counts = (
                "{ "
                f"samtools view -@ {args.threads_per_job} -bh -f 128 -F 2304 {quote(bam)} | "
                "bedtools genomecov -ibam stdin -5 -dz -strand + | "
                f"awk {quote(plus_r2_transform)}; "
                f"samtools view -@ {args.threads_per_job} -bh -f 128 -F 2304 {quote(bam)} | "
                "bedtools genomecov -ibam stdin -5 -dz -strand - | "
                f"awk {quote(minus_r2_transform)}; "
                "}"
            )
            command = " | ".join(
                [
                    r2_5prime_counts,
                    (
                        f"bedtools getfasta -fi {quote(args.fasta)} -bed stdin "
                        "-s -nameOnly -tab"
                    ),
                    f"awk {quote(fasta_filter)}",
                    (
                        f"LC_ALL=C sort --parallel={args.threads_per_job} -S 8G "
                        f"-T {quote(run_temp)} -k1,1 -k2,2n -k3,3n -k6,6"
                    ),
                    f"awk {quote(final_format)}",
                    f"bgzip -@ {args.threads_per_job} -c > {quote(temp_output)}",
                ]
            )
            with run_log.open("w", encoding="utf-8") as handle:
                handle.write(command + "\n")
                completed = subprocess.run(
                    ["bash", "-c", f"set -o pipefail; {command}"],
                    stdout=handle,
                    stderr=subprocess.STDOUT,
                )
            if completed.returncode:
                raise subprocess.CalledProcessError(completed.returncode, command)
            temp_output.replace(output)
            subprocess.run(["tabix", "-f", "-p", "bed", str(output)], check=True)

        retained = int(count_file.read_text(encoding="utf-8").strip())
        if retained > input_pairs:
            raise ValueError(f"Retained ends exceed input pairs: {retained} > {input_pairs}")
        removed = input_pairs - retained
        percent = 100.0 * removed / input_pairs if input_pairs else 0.0
        result.update(
            {
                "retained_3prime_ends": str(retained),
                "internal_priming_removed": str(removed),
                "internal_priming_percent": f"{percent:.4f}",
                "status": "PASS",
                "note": "R2_5prime_coordinate_R1_direction_internal_priming_filtered",
            }
        )
    except Exception as exc:
        result["note"] = f"{type(exc).__name__}: {exc}"
    return result


def main() -> None:
    args = parse_args()
    for tool in ("samtools", "bedtools", "bgzip", "tabix"):
        if not shutil.which(tool):
            raise RuntimeError(f"Required tool not found: {tool}")
    if not args.fasta.is_file() or not Path(f"{args.fasta}.fai").is_file():
        raise FileNotFoundError("FASTA and FASTA index are required")

    motifs = read_adapter_motifs(args.cleaning_manifest)
    with args.alignment_status.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    results: list[dict[str, str]] = []
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {
            executor.submit(process_one, row, motifs[row["run"]], args): row["run"]
            for row in rows
        }
        for index, future in enumerate(as_completed(futures), 1):
            result = future.result()
            results.append(result)
            write_status(args.status, results)
            print(
                f"[{index}/{len(rows)}] {result['run']}: {result['status']}, "
                f"retained={result['retained_3prime_ends']}, "
                f"internal={result['internal_priming_percent']}%",
                flush=True,
            )

    passed = sum(row["status"] == "PASS" for row in results)
    print(f"Completed: {passed} PASS, {len(results) - passed} FAILED")
    if passed != len(results):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
