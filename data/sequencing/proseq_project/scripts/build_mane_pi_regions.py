#!/usr/bin/env python3
"""Build strand-aware MANE TSS/TES, pause and gene-body regions."""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bed", required=True, type=Path)
    parser.add_argument("--bam", required=True, type=Path)
    parser.add_argument("--tss-tes", required=True, type=Path)
    parser.add_argument("--regions-bed", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--pause-upstream", type=int, default=50)
    parser.add_argument("--pause-downstream", type=int, default=250)
    parser.add_argument("--body-gap", type=int, default=0)
    parser.add_argument("--min-transcript-length", type=int, default=2000)
    return parser.parse_args()


def bam_contigs(bam: Path) -> dict[str, int]:
    result = subprocess.run(
        ["samtools", "view", "-H", str(bam)],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    contigs: dict[str, int] = {}
    for line in result.stdout.splitlines():
        if not line.startswith("@SQ\t"):
            continue
        values: dict[str, str] = {}
        for token in line.split("\t")[1:]:
            if ":" in token:
                key, value = token.split(":", 1)
                values[key] = value
        if "SN" in values and "LN" in values:
            contigs[values["SN"]] = int(values["LN"])
    if not contigs:
        raise ValueError("No @SQ contigs found in BAM header")
    return contigs


def split_name(name: str) -> tuple[str, str]:
    if "|" in name:
        return tuple(name.split("|", 1))  # type: ignore[return-value]
    return name, ""


def main() -> int:
    args = parse_args()
    try:
        if not args.bed.is_file():
            raise ValueError(f"BED does not exist: {args.bed}")
        if not args.bam.is_file():
            raise ValueError(f"BAM does not exist: {args.bam}")
        if args.pause_upstream < 0 or args.pause_downstream <= 0:
            raise ValueError("Invalid pause-region offsets")
        if args.body_gap < 0 or args.min_transcript_length < 0:
            raise ValueError("Negative body gap or transcript-length cutoff")

        contigs = bam_contigs(args.bam)
        args.tss_tes.parent.mkdir(parents=True, exist_ok=True)
        args.regions_bed.parent.mkdir(parents=True, exist_ok=True)
        args.metadata.parent.mkdir(parents=True, exist_ok=True)

        total = retained = skipped_short = skipped_contig = skipped_body = 0

        with (
            args.bed.open(encoding="utf-8") as source,
            args.tss_tes.open("w", encoding="utf-8", newline="") as tss_handle,
            args.regions_bed.open("w", encoding="utf-8", newline="") as bed_handle,
            args.metadata.open("w", encoding="utf-8", newline="") as metadata_handle,
        ):
            tss_writer = csv.writer(tss_handle, delimiter="\t", lineterminator="\n")
            bed_writer = csv.writer(bed_handle, delimiter="\t", lineterminator="\n")
            metadata_writer = csv.writer(
                metadata_handle, delimiter="\t", lineterminator="\n"
            )

            tss_writer.writerow(
                [
                    "gene",
                    "transcript",
                    "chrom",
                    "strand",
                    "transcript_start_bed0",
                    "transcript_end_bed0",
                    "tss_bed0",
                    "tes_bed0",
                    "transcript_length",
                    "pi_region_status",
                ]
            )
            metadata_writer.writerow(
                [
                    "gene",
                    "transcript",
                    "chrom",
                    "strand",
                    "tss_bed0",
                    "tes_bed0",
                    "pause_start",
                    "pause_end",
                    "pause_length",
                    "body_start",
                    "body_end",
                    "body_length",
                ]
            )

            for line_number, raw_line in enumerate(source, 1):
                line = raw_line.rstrip("\n")
                if not line or line.startswith(("#", "track")):
                    continue
                fields = line.split("\t")
                if len(fields) < 6:
                    raise ValueError(f"BED line {line_number} has fewer than 6 fields")

                chrom, start_text, end_text, name, score, strand = fields[:6]
                start, end = int(start_text), int(end_text)
                if start < 0 or end <= start or strand not in {"+", "-"}:
                    raise ValueError(f"Invalid BED record on line {line_number}")

                total += 1
                gene, transcript = split_name(name)
                length = end - start
                tss = start if strand == "+" else end - 1
                tes = end - 1 if strand == "+" else start
                status = "retained"

                if length < args.min_transcript_length:
                    status = "skip_short"
                    skipped_short += 1
                elif chrom not in contigs:
                    status = "skip_contig"
                    skipped_contig += 1
                else:
                    chrom_length = contigs[chrom]
                    if strand == "+":
                        pause_start = max(0, tss - args.pause_upstream)
                        pause_end = min(chrom_length, tss + args.pause_downstream)
                        body_start = tss + args.pause_downstream + args.body_gap
                        body_end = end
                    else:
                        pause_start = max(0, tss - args.pause_downstream + 1)
                        pause_end = min(
                            chrom_length, tss + args.pause_upstream + 1
                        )
                        body_start = start
                        body_end = tss - args.pause_downstream - args.body_gap + 1

                    if pause_end <= pause_start or body_end <= body_start:
                        status = "skip_empty_body"
                        skipped_body += 1
                    else:
                        retained += 1
                        pause_name = f"{gene}|{transcript}|pause"
                        body_name = f"{gene}|{transcript}|body"
                        bed_writer.writerow(
                            [
                                chrom,
                                pause_start,
                                pause_end,
                                pause_name,
                                score,
                                strand,
                            ]
                        )
                        bed_writer.writerow(
                            [chrom, body_start, body_end, body_name, score, strand]
                        )
                        metadata_writer.writerow(
                            [
                                gene,
                                transcript,
                                chrom,
                                strand,
                                tss,
                                tes,
                                pause_start,
                                pause_end,
                                pause_end - pause_start,
                                body_start,
                                body_end,
                                body_end - body_start,
                            ]
                        )

                tss_writer.writerow(
                    [
                        gene,
                        transcript,
                        chrom,
                        strand,
                        start,
                        end,
                        tss,
                        tes,
                        length,
                        status,
                    ]
                )

        print(f"Input MANE transcripts: {total}")
        print(f"Retained PI transcripts: {retained}")
        print(f"Skipped short transcripts: {skipped_short}")
        print(f"Skipped absent contigs: {skipped_contig}")
        print(f"Skipped empty gene bodies: {skipped_body}")
        print(f"TSS/TES: {args.tss_tes}")
        print(f"PI regions BED: {args.regions_bed}")
        print(f"PI region metadata: {args.metadata}")
        return 0
    except (OSError, ValueError, subprocess.CalledProcessError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
