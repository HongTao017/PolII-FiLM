#!/usr/bin/env python3
"""Build a strand-aware mask for mNET-seq RNA-processing intermediates."""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path


TRANSCRIPT_RE = re.compile(r'transcript_id "([^"]+)"')


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--tss-tes", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def main() -> None:
    args = parse_args()
    exons: dict[tuple[str, str, str], list[tuple[int, int]]] = defaultdict(list)
    with open_text(args.gtf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            match = TRANSCRIPT_RE.search(fields[8])
            if not match:
                continue
            chrom, strand = fields[0], fields[6]
            start0, end0 = int(fields[3]) - 1, int(fields[4])
            exons[(match.group(1), chrom, strand)].append((start0, end0))

    mask: set[tuple[str, int, str, str]] = set()
    for (transcript, chrom, strand), intervals in exons.items():
        ordered = sorted(intervals, reverse=(strand == "-"))
        for upstream, downstream in zip(ordered, ordered[1:]):
            if strand == "+":
                exon_last = upstream[1] - 1
                intron_last = downstream[0] - 1
            else:
                exon_last = upstream[0]
                intron_last = downstream[1]
            if exon_last >= 0:
                mask.add((chrom, exon_last, strand, "exon_5ss_cleavage"))
            if intron_last >= 0:
                mask.add((chrom, intron_last, strand, "intron_3ss_cleavage"))

    with args.tss_tes.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            mask.add((row["chrom"], int(row["tes_bed0"]), row["strand"], "annotated_TES"))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        for chrom, position, strand, reason in sorted(mask, key=lambda item: (item[0], item[1], item[2], item[3])):
            handle.write(f"{chrom}\t{position}\t{position + 1}\t{reason}\t0\t{strand}\n")
    print(f"Wrote {len(mask)} strand-aware processing positions to {args.output}")


if __name__ == "__main__":
    main()
