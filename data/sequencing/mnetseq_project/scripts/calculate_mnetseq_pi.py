#!/usr/bin/env python3
"""Calculate strand-specific MANE pausing indices from masked mNET-seq 3' ends."""

from __future__ import annotations

import argparse
import csv
import math
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-status", required=True, type=Path)
    parser.add_argument("--formal-names", required=True, type=Path)
    parser.add_argument("--tss-tes", required=True, type=Path)
    parser.add_argument("--fai", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--matrix", required=True, type=Path)
    parser.add_argument("--pause-upstream", type=int, default=50)
    parser.add_argument("--pause-downstream", type=int, default=250)
    parser.add_argument("--min-transcript-length", type=int, default=2000)
    parser.add_argument("--minimum-total-count", type=int, default=10)
    parser.add_argument("--jobs", type=int, default=4)
    return parser.parse_args()


def read_chrom_sizes(path: Path) -> dict[str, int]:
    sizes: dict[str, int] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            chrom, length, *_ = line.rstrip("\n").split("\t")
            sizes[chrom] = int(length)
    if not sizes:
        raise ValueError(f"No chromosome sizes found in {path}")
    return sizes


def build_regions(
    tss_tes: Path,
    chrom_sizes: dict[str, int],
    regions_bed: Path,
    metadata_tsv: Path,
    upstream: int,
    downstream: int,
    min_length: int,
) -> list[dict[str, str]]:
    transcripts: list[dict[str, str]] = []
    region_rows: list[tuple[str, int, int, str, int, str, int, str]] = []

    with tss_tes.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"gene_transcript", "chrom", "strand", "tss_bed0", "tes_bed0"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"TSS/TES table lacks required columns: {sorted(required)}")

        for row in reader:
            chrom, strand = row["chrom"], row["strand"]
            if chrom not in chrom_sizes or strand not in {"+", "-"}:
                continue
            tss, tes = int(row["tss_bed0"]), int(row["tes_bed0"])
            transcript_length = abs(tes - tss) + 1
            if transcript_length < min_length:
                continue

            if strand == "+":
                pause_start = max(0, tss - upstream)
                pause_end = min(chrom_sizes[chrom], tss + downstream)
                body_start = pause_end
                body_end = min(chrom_sizes[chrom], tes + 1)
            else:
                pause_start = max(0, tss - downstream + 1)
                pause_end = min(chrom_sizes[chrom], tss + upstream + 1)
                body_start = max(0, tes)
                body_end = pause_start

            if pause_end <= pause_start or body_end <= body_start:
                continue

            name = row["gene_transcript"]
            gene, transcript = name.split("|", 1) if "|" in name else (name, name)
            index = len(transcripts)
            transcripts.append(
                {
                    "tx_index": str(index),
                    "gene": gene,
                    "transcript": transcript,
                    "chrom": chrom,
                    "strand": strand,
                    "tss_bed0": str(tss),
                    "tes_bed0": str(tes),
                    "transcript_length": str(transcript_length),
                    "pause_length": str(pause_end - pause_start),
                    "body_length": str(body_end - body_start),
                }
            )
            region_rows.extend(
                [
                    (chrom, pause_start, pause_end, f"{index}|pause", 0, strand, index, "pause"),
                    (chrom, body_start, body_end, f"{index}|body", 0, strand, index, "body"),
                ]
            )

    if not transcripts:
        raise ValueError("No MANE transcripts remained after region filtering")

    region_rows.sort(key=lambda row: (row[0], row[1], row[2], row[5], row[7]))
    with regions_bed.open("w", encoding="utf-8", newline="") as handle:
        for row in region_rows:
            handle.write("\t".join(map(str, row)) + "\n")
    with metadata_tsv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(transcripts[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(transcripts)
    return transcripts


def read_formal_names(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"old_sample_id", "formal_sample_name"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"Formal-name table lacks required columns: {sorted(required)}")
    mapping = {row["old_sample_id"]: row["formal_sample_name"] for row in rows}
    if len(mapping) != len(rows) or len(set(mapping.values())) != len(mapping):
        raise ValueError("Duplicate old or formal sample names were found")
    return mapping


def normalize_sample(row: dict[str, str], names: dict[str, str]) -> dict[str, str]:
    original_id = row.get("sample_id") or row.get("sample") or ""
    if not original_id:
        raise ValueError("Sample status row has no sample/sample_id value")
    if original_id not in names:
        raise ValueError(f"No formal sample name for {original_id}")
    parts = original_id.split("_", 1)
    return {
        "sample_id": names[original_id],
        "original_sample_id": original_id,
        "gse": row.get("gse") or parts[0],
        "gsm": row.get("gsm") or (parts[1] if len(parts) == 2 else "NA"),
        "technical_run_count": row.get("technical_run_count", "NA"),
        "signal": row.get("output_bed_gz", ""),
        "status": row.get("status", ""),
        "removed_percent": row.get("removed_percent", "NA"),
    }


def calculate_one(
    sample: dict[str, str],
    transcripts: list[dict[str, str]],
    regions_bed: Path,
    output_dir: Path,
    active_dir: Path,
    model_dir: Path,
    minimum_total: int,
) -> tuple[dict[str, str], dict[str, float | None]]:
    sample_id = sample["sample_id"]
    output = output_dir / f"{sample_id}.mNET.PI.tsv"
    active_output = active_dir / f"{sample_id}.mNET.PI.active_min10.tsv"
    model_output = model_dir / f"{sample_id}.mNET.PI.model_ready.tsv"
    mapped = output_dir / f"{sample_id}.mapped_regions.tmp.tsv"
    signal = Path(sample["signal"])
    result = {
        "sample_id": sample_id,
        "original_sample_id": sample["original_sample_id"],
        "gse": sample["gse"],
        "gsm": sample["gsm"],
        "technical_run_count": sample["technical_run_count"],
        "processing_removed_percent": sample["removed_percent"],
        "total_transcripts": str(len(transcripts)),
        "active_min10_transcripts": "NA",
        "model_ready_transcripts": "NA",
        "output_pi": str(output),
        "status": "FAILED",
        "note": "",
    }
    matrix_values: dict[str, float | None] = {}
    fieldnames = [
        "gene", "transcript", "chrom", "strand", "tss_bed0", "tes_bed0",
        "transcript_length", "pause_length", "pause_count", "body_length",
        "body_count", "total_count", "pause_density", "body_density", "pi_raw",
        "log2_pi_pseudocount", "active_min10", "model_ready",
    ]

    try:
        if sample["status"] != "PASS" or not signal.is_file():
            raise ValueError(f"Masked sample signal is unavailable: {signal}")
        with mapped.open("w", encoding="utf-8") as handle:
            subprocess.run(
                ["bedtools", "map", "-a", str(regions_bed), "-b", str(signal),
                 "-s", "-c", "5", "-o", "sum", "-null", "0"],
                check=True,
                stdout=handle,
                stderr=subprocess.PIPE,
                text=True,
            )

        counts = [{"pause": 0, "body": 0} for _ in transcripts]
        with mapped.open(encoding="utf-8") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                counts[int(fields[6])][fields[7]] = int(float(fields[8]))

        active_rows: list[dict[str, object]] = []
        model_rows: list[dict[str, object]] = []
        active = model_ready = 0
        with output.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for transcript, count in zip(transcripts, counts):
                pause_count, body_count = count["pause"], count["body"]
                pause_length = int(transcript["pause_length"])
                body_length = int(transcript["body_length"])
                total = pause_count + body_count
                is_active = total >= minimum_total
                is_model_ready = is_active and body_count > 0
                active += int(is_active)
                model_ready += int(is_model_ready)

                pause_density = pause_count / pause_length
                body_density = body_count / body_length
                pi_raw = pause_density / body_density if body_count > 0 else None
                log2_pi = math.log2(
                    ((pause_count + 1) / pause_length)
                    / ((body_count + 1) / body_length)
                )
                row: dict[str, object] = {
                    "gene": transcript["gene"],
                    "transcript": transcript["transcript"],
                    "chrom": transcript["chrom"],
                    "strand": transcript["strand"],
                    "tss_bed0": transcript["tss_bed0"],
                    "tes_bed0": transcript["tes_bed0"],
                    "transcript_length": transcript["transcript_length"],
                    "pause_length": pause_length,
                    "pause_count": pause_count,
                    "body_length": body_length,
                    "body_count": body_count,
                    "total_count": total,
                    "pause_density": f"{pause_density:.9g}",
                    "body_density": f"{body_density:.9g}",
                    "pi_raw": "NA" if pi_raw is None else f"{pi_raw:.9g}",
                    "log2_pi_pseudocount": f"{log2_pi:.9g}",
                    "active_min10": "yes" if is_active else "no",
                    "model_ready": "yes" if is_model_ready else "no",
                }
                writer.writerow(row)
                if is_active:
                    active_rows.append(row)
                if is_model_ready:
                    model_rows.append(row)
                key = f"{transcript['gene']}|{transcript['transcript']}"
                matrix_values[key] = log2_pi if is_model_ready else None

        for path, rows in ((active_output, active_rows), (model_output, model_rows)):
            with path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)

        mapped.unlink(missing_ok=True)
        result.update(
            {
                "active_min10_transcripts": str(active),
                "model_ready_transcripts": str(model_ready),
                "status": "PASS",
                "note": "masked_3prime_TSS_minus50_plus250_body_to_TES_strand_specific",
            }
        )
    except Exception as exc:
        result["note"] = f"{type(exc).__name__}: {exc}"
    return result, matrix_values


def main() -> int:
    args = parse_args()
    try:
        if not shutil.which("bedtools"):
            raise RuntimeError("bedtools is required")
        for path in (args.sample_status, args.formal_names, args.tss_tes, args.fai):
            if not path.is_file():
                raise FileNotFoundError(path)
        if args.pause_upstream < 0 or args.pause_downstream <= 0:
            raise ValueError("Invalid pause-region offsets")

        args.output_dir.mkdir(parents=True, exist_ok=True)
        args.summary.parent.mkdir(parents=True, exist_ok=True)
        args.matrix.parent.mkdir(parents=True, exist_ok=True)
        active_dir = args.output_dir.parent / "active_min10"
        model_dir = args.output_dir.parent / "model_ready"
        active_dir.mkdir(parents=True, exist_ok=True)
        model_dir.mkdir(parents=True, exist_ok=True)

        regions_bed = args.output_dir.parent / "MANE_Select.mNET.PI_regions.bed"
        metadata_tsv = args.output_dir.parent / "MANE_Select.mNET.PI_regions.tsv"
        transcripts = build_regions(
            args.tss_tes, read_chrom_sizes(args.fai), regions_bed, metadata_tsv,
            args.pause_upstream, args.pause_downstream, args.min_transcript_length,
        )
        with args.sample_status.open(newline="", encoding="utf-8-sig") as handle:
            names = read_formal_names(args.formal_names)
            samples = [
                normalize_sample(row, names)
                for row in csv.DictReader(handle, delimiter="\t")
            ]
        if not samples:
            raise ValueError("No samples found in processing-mask status table")

        summaries: list[dict[str, str]] = []
        matrix: dict[str, dict[str, float | None]] = {}
        with ThreadPoolExecutor(max_workers=args.jobs) as executor:
            futures = [
                executor.submit(
                    calculate_one, sample, transcripts, regions_bed, args.output_dir,
                    active_dir, model_dir, args.minimum_total_count,
                )
                for sample in samples
            ]
            for index, future in enumerate(as_completed(futures), 1):
                summary, values = future.result()
                summaries.append(summary)
                matrix[summary["sample_id"]] = values
                print(
                    f"[{index}/{len(samples)}] {summary['sample_id']}: {summary['status']}, "
                    f"active={summary['active_min10_transcripts']}, "
                    f"model={summary['model_ready_transcripts']}",
                    flush=True,
                )

        with args.summary.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(summaries[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(sorted(summaries, key=lambda row: row["sample_id"]))

        sample_ids = sorted(matrix)
        all_keys = [f"{row['gene']}|{row['transcript']}" for row in transcripts]
        # Match the GRO/PRO matrix rule: retain the union of transcripts with
        # model-ready PI in at least one sample, rather than all MANE rows.
        keys = [
            key
            for key in all_keys
            if any(matrix[sample].get(key) is not None for sample in sample_ids)
        ]
        with args.matrix.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["gene_transcript", *sample_ids])
            for key in keys:
                writer.writerow(
                    [key]
                    + [
                        "NA" if matrix[sample].get(key) is None
                        else f"{matrix[sample][key]:.9g}"
                        for sample in sample_ids
                    ]
                )

        failed = [row for row in summaries if row["status"] != "PASS"]
        print(f"MANE PI transcripts: {len(transcripts)}")
        print(f"Model-ready matrix transcripts: {len(keys)}")
        print(f"Completed: {len(summaries) - len(failed)} PASS, {len(failed)} FAILED")
        return 1 if failed else 0
    except Exception as exc:
        print(f"ERROR: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
