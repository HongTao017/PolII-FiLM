#!/usr/bin/env python3
"""Calculate strand-specific promoter-proximal pausing indices.

The default regions reproduce the definition used by Fisher et al.:

  pause region: TSS - 50 bp to TSS + 250 bp (300 bp, BED half-open)
  gene body:    immediately after the pause region to TES

Input annotations must be BED6 transcript intervals. The fourth column is used
as the transcript identifier and may contain a value such as GENE|NM_000000.1.
The script calls ``bedtools multicov`` with same- or opposite-strand counting,
so the BAM must be coordinate-sorted and indexed and BED chromosome names must
match the BAM header. Paired-end libraries can be restricted to read 1.
"""

from __future__ import annotations

import argparse
import csv
import math
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


@dataclass(frozen=True)
class Transcript:
    chrom: str
    start: int
    end: int
    name: str
    score: str
    strand: str

    @property
    def length(self) -> int:
        return self.end - self.start

    @property
    def tss(self) -> int:
        return self.start if self.strand == "+" else self.end - 1

    @property
    def tes(self) -> int:
        return self.end - 1 if self.strand == "+" else self.start


@dataclass(frozen=True)
class Region:
    chrom: str
    start: int
    end: int
    region_id: str
    score: str
    strand: str
    transcript_index: int
    kind: str

    @property
    def length(self) -> int:
        return self.end - self.start

    def bed_fields(self) -> List[str]:
        return [
            self.chrom,
            str(self.start),
            str(self.end),
            self.region_id,
            self.score,
            self.strand,
        ]


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate strand-specific pausing indices from a BED6 transcript "
            "annotation and a coordinate-sorted, indexed BAM."
        )
    )
    parser.add_argument("--bed", required=True, type=Path, help="MANE BED6 file")
    parser.add_argument("--bam", required=True, type=Path, help="Sorted and indexed BAM")
    parser.add_argument("--output", required=True, type=Path, help="Output TSV file")
    parser.add_argument(
        "--pause-upstream",
        type=int,
        default=50,
        help="Bases upstream of TSS in the pause region (default: 50)",
    )
    parser.add_argument(
        "--pause-downstream",
        type=int,
        default=250,
        help="Bases downstream of TSS in the pause region (default: 250)",
    )
    parser.add_argument(
        "--body-gap",
        type=int,
        default=0,
        help=(
            "Additional uncounted bases after the pause region before the gene "
            "body (default: 0; use 250 to start the body at TSS+500)"
        ),
    )
    parser.add_argument(
        "--min-transcript-length",
        type=int,
        default=2000,
        help="Discard transcripts shorter than this many bp (default: 2000; 0 disables)",
    )
    parser.add_argument(
        "--mapq",
        type=int,
        default=20,
        help="Minimum BAM mapping quality passed to bedtools (default: 20)",
    )
    parser.add_argument(
        "--library-direction",
        choices=("same", "opposite"),
        default="same",
        help=(
            "Relationship between BAM alignment strand and transcript strand "
            "(default: same)"
        ),
    )
    parser.add_argument(
        "--read-selection",
        choices=("all", "r1"),
        default="all",
        help="Count all alignments or only read 1 from paired-end BAMs (default: all)",
    )
    parser.add_argument(
        "--pseudocount",
        type=float,
        default=1.0,
        help="Count pseudocount for stabilized PI columns (default: 1)",
    )
    parser.add_argument(
        "--keep-temporary-bed",
        type=Path,
        help="Optionally save the generated pause/body BED6 regions here",
    )
    return parser.parse_args(argv)


def require_inputs(args: argparse.Namespace) -> None:
    if not args.bed.is_file():
        raise ValueError(f"BED file does not exist: {args.bed}")
    if not args.bam.is_file():
        raise ValueError(f"BAM file does not exist: {args.bam}")
    bai_candidates = [
        Path(str(args.bam) + ".bai"),
        args.bam.with_suffix(".bai"),
        Path(str(args.bam) + ".csi"),
        args.bam.with_suffix(".csi"),
    ]
    if not any(path.is_file() for path in bai_candidates):
        raise ValueError(
            "BAM index was not found. Run: samtools index " + str(args.bam)
        )
    if shutil.which("bedtools") is None:
        raise ValueError("bedtools was not found in PATH")
    if shutil.which("samtools") is None:
        raise ValueError("samtools was not found in PATH")
    if args.pause_upstream < 0 or args.pause_downstream <= 0:
        raise ValueError("Pause offsets must be non-negative and downstream must be > 0")
    if args.body_gap < 0 or args.min_transcript_length < 0:
        raise ValueError("Body gap and minimum transcript length cannot be negative")
    if args.mapq < 0 or args.pseudocount < 0:
        raise ValueError("MAPQ and pseudocount cannot be negative")


def read_bam_contigs(bam: Path) -> Dict[str, int]:
    result = subprocess.run(
        ["samtools", "view", "-H", str(bam)],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    contigs: Dict[str, int] = {}
    for line in result.stdout.splitlines():
        if not line.startswith("@SQ\t"):
            continue
        fields = dict(
            token.split(":", 1) for token in line.split("\t")[1:] if ":" in token
        )
        if "SN" in fields and "LN" in fields:
            contigs[fields["SN"]] = int(fields["LN"])
    if not contigs:
        raise ValueError("No @SQ chromosome records were found in the BAM header")
    return contigs


def read_transcripts(path: Path) -> List[Transcript]:
    transcripts: List[Transcript] = []
    seen_names = set()
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("track"):
                continue
            fields = line.split("\t")
            if len(fields) < 6:
                raise ValueError(f"BED line {line_number} has fewer than 6 columns")
            chrom, start_text, end_text, name, score, strand = fields[:6]
            try:
                start, end = int(start_text), int(end_text)
            except ValueError as exc:
                raise ValueError(f"Invalid BED coordinates on line {line_number}") from exc
            if start < 0 or end <= start:
                raise ValueError(f"Invalid BED interval on line {line_number}: {start}-{end}")
            if strand not in {"+", "-"}:
                raise ValueError(f"Invalid strand on BED line {line_number}: {strand}")
            if name in seen_names:
                raise ValueError(
                    f"Duplicate transcript identifier on BED line {line_number}: {name}"
                )
            seen_names.add(name)
            transcripts.append(Transcript(chrom, start, end, name, score, strand))
    if not transcripts:
        raise ValueError("No transcript records were found in the BED file")
    return transcripts


def make_regions(
    tx: Transcript,
    tx_index: int,
    chrom_length: int,
    pause_upstream: int,
    pause_downstream: int,
    body_gap: int,
) -> Tuple[Region, Region]:
    """Create adjacent strand-aware BED regions.

    For a plus-strand TSS at x, the pause interval is
    [x-upstream, x+downstream). For a minus-strand TSS at x, the same
    transcription-relative interval is [x-downstream+1, x+upstream+1).
    """
    if tx.strand == "+":
        pause_start = max(0, tx.tss - pause_upstream)
        pause_end = min(chrom_length, tx.tss + pause_downstream)
        body_start = tx.tss + pause_downstream + body_gap
        body_end = tx.end
    else:
        pause_start = max(0, tx.tss - pause_downstream + 1)
        pause_end = min(chrom_length, tx.tss + pause_upstream + 1)
        body_start = tx.start
        body_end = tx.tss - pause_downstream - body_gap + 1

    if pause_end <= pause_start:
        raise ValueError(f"Empty pause region for {tx.name}")
    if body_end <= body_start:
        raise ValueError(f"Empty gene-body region for {tx.name}")

    pause = Region(
        tx.chrom,
        pause_start,
        pause_end,
        f"{tx_index}|pause",
        "0",
        tx.strand,
        tx_index,
        "pause",
    )
    body = Region(
        tx.chrom,
        body_start,
        body_end,
        f"{tx_index}|body",
        "0",
        tx.strand,
        tx_index,
        "body",
    )
    return pause, body


def write_regions(path: Path, regions: Iterable[Region]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for region in regions:
            writer.writerow(region.bed_fields())


def count_regions(
    region_bed: Path,
    bam: Path,
    mapq: int,
    library_direction: str,
) -> Dict[str, int]:
    strand_option = "-s" if library_direction == "same" else "-S"
    command = [
        "bedtools",
        "multicov",
        strand_option,
        "-q",
        str(mapq),
        "-bams",
        str(bam),
        "-bed",
        str(region_bed),
    ]
    result = subprocess.run(
        command,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    counts: Dict[str, int] = {}
    for line_number, line in enumerate(result.stdout.splitlines(), 1):
        fields = line.split("\t")
        if len(fields) < 7:
            raise ValueError(f"Unexpected bedtools output on line {line_number}")
        region_id = fields[3]
        counts[region_id] = int(fields[-1])
    return counts


def select_read1(source_bam: Path, output_bam: Path) -> None:
    """Create a coordinate-sorted BAM containing only read 1 alignments."""
    subprocess.run(
        [
            "samtools",
            "view",
            "-bh",
            "-f",
            "64",
            "-o",
            str(output_bam),
            str(source_bam),
        ],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    subprocess.run(
        ["samtools", "index", str(output_bam)],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def split_name(name: str) -> Tuple[str, str]:
    if "|" not in name:
        return name, ""
    return tuple(name.split("|", 1))  # type: ignore[return-value]


def finite_or_text(value: float) -> str:
    if math.isnan(value):
        return "NA"
    if math.isinf(value):
        return "Inf" if value > 0 else "-Inf"
    return f"{value:.10g}"


def write_results(
    output: Path,
    kept: List[Tuple[Transcript, Region, Region]],
    counts: Dict[str, int],
    pseudocount: float,
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "gene",
        "transcript",
        "chrom",
        "strand",
        "transcript_start",
        "transcript_end",
        "tss_bed0",
        "tes_bed0",
        "transcript_length",
        "pause_start",
        "pause_end",
        "pause_length",
        "pause_count",
        "pause_density",
        "body_start",
        "body_end",
        "body_length",
        "body_count",
        "body_density",
        "pi_raw",
        "pi_pseudocount",
        "log2_pi_pseudocount",
    ]
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        for tx, pause, body in kept:
            pause_count = counts[pause.region_id]
            body_count = counts[body.region_id]
            pause_density = pause_count / pause.length
            body_density = body_count / body.length
            if body_density == 0:
                pi_raw = math.inf if pause_density > 0 else math.nan
            else:
                pi_raw = pause_density / body_density

            pause_density_pc = (pause_count + pseudocount) / pause.length
            body_density_pc = (body_count + pseudocount) / body.length
            if body_density_pc == 0:
                pi_pc = math.inf if pause_density_pc > 0 else math.nan
            else:
                pi_pc = pause_density_pc / body_density_pc
            log2_pi_pc = math.log2(pi_pc) if pi_pc > 0 and math.isfinite(pi_pc) else math.nan
            gene, transcript = split_name(tx.name)
            writer.writerow(
                [
                    gene,
                    transcript,
                    tx.chrom,
                    tx.strand,
                    tx.start,
                    tx.end,
                    tx.tss,
                    tx.tes,
                    tx.length,
                    pause.start,
                    pause.end,
                    pause.length,
                    pause_count,
                    finite_or_text(pause_density),
                    body.start,
                    body.end,
                    body.length,
                    body_count,
                    finite_or_text(body_density),
                    finite_or_text(pi_raw),
                    finite_or_text(pi_pc),
                    finite_or_text(log2_pi_pc),
                ]
            )


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    try:
        require_inputs(args)
        contigs = read_bam_contigs(args.bam)
        transcripts = read_transcripts(args.bed)

        kept: List[Tuple[Transcript, Region, Region]] = []
        skipped_short = 0
        skipped_contig = 0
        skipped_body = 0
        regions: List[Region] = []
        for tx in transcripts:
            if tx.length < args.min_transcript_length:
                skipped_short += 1
                continue
            if tx.chrom not in contigs:
                skipped_contig += 1
                continue
            try:
                pause, body = make_regions(
                    tx,
                    len(kept),
                    contigs[tx.chrom],
                    args.pause_upstream,
                    args.pause_downstream,
                    args.body_gap,
                )
            except ValueError:
                skipped_body += 1
                continue
            kept.append((tx, pause, body))
            regions.extend((pause, body))

        if not kept:
            raise ValueError("No transcripts remained after filtering")

        with tempfile.TemporaryDirectory(prefix="calculate_pi_") as temp_dir:
            region_bed = Path(temp_dir) / "pi_regions.bed"
            write_regions(region_bed, regions)
            if args.keep_temporary_bed:
                args.keep_temporary_bed.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(region_bed, args.keep_temporary_bed)
            count_bam = args.bam
            if args.read_selection == "r1":
                count_bam = Path(temp_dir) / "read1.bam"
                select_read1(args.bam, count_bam)
            counts = count_regions(
                region_bed,
                count_bam,
                args.mapq,
                args.library_direction,
            )

        expected = len(regions)
        if len(counts) != expected:
            raise ValueError(
                f"Expected {expected} regional counts but bedtools returned {len(counts)}"
            )
        write_results(args.output, kept, counts, args.pseudocount)

        print(f"Input transcripts: {len(transcripts)}", file=sys.stderr)
        print(f"Retained transcripts: {len(kept)}", file=sys.stderr)
        print(f"Skipped as short: {skipped_short}", file=sys.stderr)
        print(f"Skipped because chromosome was absent from BAM: {skipped_contig}", file=sys.stderr)
        print(f"Skipped because gene body was empty: {skipped_body}", file=sys.stderr)
        print(f"Library direction: {args.library_direction}", file=sys.stderr)
        print(f"Read selection: {args.read_selection}", file=sys.stderr)
        print(f"Output: {args.output}", file=sys.stderr)
        return 0
    except (ValueError, OSError, subprocess.CalledProcessError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        if isinstance(exc, subprocess.CalledProcessError) and exc.stderr:
            print(exc.stderr.strip(), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
