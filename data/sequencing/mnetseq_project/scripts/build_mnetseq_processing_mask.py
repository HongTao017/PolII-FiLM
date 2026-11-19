#!/usr/bin/env python3
"""Build strand-aware mNET-seq splice and 3-prime cleavage masks from MANE GTF."""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path


ATTRIBUTE_RE = re.compile(r'([A-Za-z0-9_.-]+)\s+"([^"]*)"')


@dataclass
class Transcript:
    chrom: str
    strand: str
    start0: int
    end0: int
    gene: str
    transcript_id: str
    exons: list[tuple[int, int]] = field(default_factory=list)


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def attributes(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for key, value in ATTRIBUTE_RE.findall(text):
        if key in result:
            result[key] = f"{result[key]};{value}"
        else:
            result[key] = value
    return result


def is_mane_select(attrs: dict[str, str]) -> bool:
    return "MANE Select" in attrs.get("tag", "")


def chromosome_key(chrom: str) -> tuple[int, str]:
    token = chrom.removeprefix("chr")
    if token.isdigit():
        return int(token), ""
    if token == "X":
        return 23, ""
    if token == "Y":
        return 24, ""
    return 1000, chrom


def load_mane_transcripts(path: Path) -> dict[str, Transcript]:
    transcripts: dict[str, Transcript] = {}
    with open_text(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "transcript":
                continue
            attrs = attributes(fields[8])
            if not is_mane_select(attrs):
                continue
            transcript_id = attrs.get("transcript_id", "")
            if not transcript_id:
                continue
            gene = (
                attrs.get("gene_name")
                or attrs.get("gene")
                or attrs.get("gene_id")
                or "NA"
            )
            transcripts[transcript_id] = Transcript(
                chrom=fields[0],
                strand=fields[6],
                start0=int(fields[3]) - 1,
                end0=int(fields[4]),
                gene=gene,
                transcript_id=transcript_id,
            )

    with open_text(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "exon":
                continue
            attrs = attributes(fields[8])
            transcript_id = attrs.get("transcript_id", "")
            if transcript_id not in transcripts:
                continue
            transcripts[transcript_id].exons.append(
                (int(fields[3]) - 1, int(fields[4]))
            )

    empty = [tx_id for tx_id, tx in transcripts.items() if not tx.exons]
    if empty:
        raise ValueError(f"MANE transcripts without exons: {len(empty)}")
    if not transcripts:
        raise ValueError("No MANE Select transcripts were found in the GTF")
    return transcripts


def interval(position: int, flank: int) -> tuple[int, int]:
    return max(0, position - flank), position + flank + 1


def build_sites(
    transcripts: dict[str, Transcript], flank: int
) -> tuple[list[tuple[str, int, int, str, str, str]], list[tuple[str, int, int, str, str, str]]]:
    splice_sites: set[tuple[str, int, int, str, str, str]] = set()
    cleavage_sites: set[tuple[str, int, int, str, str, str]] = set()

    for tx in transcripts.values():
        exons = sorted(set(tx.exons))
        tes_position = tx.end0 - 1 if tx.strand == "+" else tx.start0
        start, end = interval(tes_position, flank)
        cleavage_sites.add(
            (
                tx.chrom,
                start,
                end,
                f"{tx.gene}|{tx.transcript_id}|three_prime_cleavage",
                "three_prime_cleavage",
                tx.strand,
            )
        )

        for left, right in zip(exons, exons[1:]):
            intron_start = left[1]
            intron_end = right[0]
            if intron_end <= intron_start:
                continue
            if tx.strand == "+":
                exon_3prime = intron_start - 1
                intron_3prime = intron_end - 1
            else:
                exon_3prime = intron_end
                intron_3prime = intron_start

            for position, site_type in (
                (exon_3prime, "splice_exon_3prime"),
                (intron_3prime, "splice_intron_3prime"),
            ):
                start, end = interval(position, flank)
                splice_sites.add(
                    (
                        tx.chrom,
                        start,
                        end,
                        f"{tx.gene}|{tx.transcript_id}|{site_type}",
                        site_type,
                        tx.strand,
                    )
                )

    key = lambda row: (chromosome_key(row[0]), row[1], row[2], row[5], row[4], row[3])
    return sorted(splice_sites, key=key), sorted(cleavage_sites, key=key)


def write_bed(path: Path, rows: list[tuple[str, int, int, str, str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for chrom, start, end, name, _site_type, strand in rows:
            writer.writerow((chrom, start, end, name, 0, strand))


def write_combined(
    path: Path,
    rows: list[tuple[str, int, int, str, str, str]],
) -> int:
    grouped: dict[tuple[str, int, int, str], dict[str, set[str]]] = defaultdict(
        lambda: {"names": set(), "types": set()}
    )
    for chrom, start, end, name, site_type, strand in rows:
        group = grouped[(chrom, start, end, strand)]
        group["names"].add(name)
        group["types"].add(site_type)

    output_rows = []
    for (chrom, start, end, strand), values in grouped.items():
        site_types = ",".join(sorted(values["types"]))
        names = ",".join(sorted(values["names"]))
        output_rows.append((chrom, start, end, f"{site_types}|{names}", 0, strand))
    output_rows.sort(
        key=lambda row: (chromosome_key(row[0]), row[1], row[2], row[5], row[3])
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows(output_rows)
    return len(output_rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--flank",
        type=int,
        default=0,
        help="Bases added on each side; default 0 masks only the exact nucleotide",
    )
    args = parser.parse_args()
    if args.flank < 0:
        raise ValueError("--flank cannot be negative")

    transcripts = load_mane_transcripts(args.gtf)
    splice, cleavage = build_sites(transcripts, args.flank)
    splice_path = args.output_dir / "MANE_Select.splice_processing_sites.bed"
    cleavage_path = args.output_dir / "MANE_Select.three_prime_cleavage_sites.bed"
    combined_path = args.output_dir / "MANE_Select.mNET_processing_mask.bed"
    summary_path = args.output_dir / "MANE_Select.mNET_processing_mask.summary.tsv"

    write_bed(splice_path, splice)
    write_bed(cleavage_path, cleavage)
    combined_count = write_combined(combined_path, splice + cleavage)

    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("metric", "value"))
        writer.writerow(("mane_select_transcripts", len(transcripts)))
        writer.writerow(("splice_processing_records", len(splice)))
        writer.writerow(("three_prime_cleavage_records", len(cleavage)))
        writer.writerow(("combined_unique_intervals", combined_count))
        writer.writerow(("flank_each_side_bp", args.flank))

    print(f"MANE Select transcripts: {len(transcripts)}")
    print(f"Splice processing sites: {len(splice)}")
    print(f"3-prime cleavage sites: {len(cleavage)}")
    print(f"Combined unique mask intervals: {combined_count}")
    print(f"Output: {combined_path}")


if __name__ == "__main__":
    main()
