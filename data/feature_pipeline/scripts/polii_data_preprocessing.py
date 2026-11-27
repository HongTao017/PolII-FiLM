#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"Construct gene-level pausing predictors in the specified feature order.\n\nInputs comprise MANE GTF, GRCh38 FASTA/FAI and optional PC1 and PI tables. Outputs include gene-by-feature and feature-by-gene matrices, a scaled architecture-feature view and annotation matching records. Features cover transcript/exon/intron lengths, UTR/CDS composition, CpG statistics, overlapping k-mer frequencies and in-frame codon frequencies. Missing external assays remain unavailable rather than becoming artificial zeros. PI sample columns must be excluded when predicting their derived PC1 target.\n"

from __future__ import annotations

import argparse
import gzip
import itertools
import math
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


DNA_ALPHABET = "ACGT"
STOP_CODONS = {"TAA", "TAG", "TGA"}
SENSE_CODONS = [
    "".join(parts)
    for parts in itertools.product(DNA_ALPHABET, repeat=3)
    if "".join(parts) not in STOP_CODONS
]


BASE_FEATURE_ORDER = [
    "transcript",
    "chrom",
    "strand",
    "tss_bed0",
    "tes_bed0",
    "gene_biotype",
    "protein_coding",
    "transcript_span_bp",
    "log10_transcript_span_bp",
    "spliced_length_bp",
    "exon_total_length_bp",
    "exon_count",
    "mean_exon_length_bp",
    "min_exon_length_bp",
    "max_exon_length_bp",
    "first_exon_length_bp",
    "last_exon_length_bp",
    "intron_count",
    "intron_total_length_bp",
    "mean_intron_length_bp",
    "min_intron_length_bp",
    "max_intron_length_bp",
    "first_intron_length_bp",
    "last_intron_length_bp",
    "tss_to_first_splice_bp",
    "cds_length_bp",
    "utr5_length_bp",
    "utr3_length_bp",
    "cds_fraction",
    "utr5_fraction",
    "utr3_fraction",
    "exon_fraction",
    "intron_fraction",
    "orf_exon_junction_count",
    "orf_exon_junction_density_per_kb",
    "promoter_gc_fraction",
    "promoter_at_fraction",
    "pause_gc_fraction",
    "pause_at_fraction",
    "early_body_gc_fraction",
    "early_body_at_fraction",
    "spliced_gc_fraction",
    "utr5_gc_fraction",
    "cds_gc_fraction",
    "utr3_gc_fraction",
    "promoter_cpg_count",
    "pause_cpg_count",
    "early_body_cpg_count",
    "promoter_cpg_oe",
    "pause_cpg_oe",
    "early_body_cpg_oe",
]


@dataclass
class TranscriptRecord:
    "Transcript annotation intervals parsed from GTF, stored as zero-based half-open coordinates."

    transcript_id: str
    gene_id: str = ""
    gene_name: str = ""
    chrom: str = ""
    strand: str = "."
    start0: Optional[int] = None
    end0: Optional[int] = None
    gene_biotype: str = ""
    transcript_biotype: str = ""
    tags: List[str] = field(default_factory=list)
    aliases: List[str] = field(default_factory=list)
    exons: List[Tuple[int, int]] = field(default_factory=list)
    cds: List[Tuple[int, int]] = field(default_factory=list)

    @property
    def span_bp(self) -> int:
        if self.start0 is None or self.end0 is None:
            return 0
        return max(0, self.end0 - self.start0)

    @property
    def protein_coding(self) -> int:
        "Unify coding-transcript labels across annotation providers using CDS evidence."
        values = {
            self.gene_biotype.strip().lower(),
            self.transcript_biotype.strip().lower(),
        }
        coding_labels = {"protein_coding", "protein-coding", "mrna"}
        return int(bool(values & coding_labels))


class IndexedFasta:
    "Read genomic intervals with a FAI index without loading the entire reference into memory."

    def __init__(self, fasta_path: Path, fai_path: Optional[Path] = None) -> None:
        self.fasta_path = fasta_path
        self.fai_path = fai_path or Path(f"{fasta_path}.fai")
        if not self.fasta_path.exists():
            raise FileNotFoundError(f"FASTA file not found: {self.fasta_path}")
        if not self.fai_path.exists():
            raise FileNotFoundError(f"FAI file not found: {self.fai_path}")

        self.index: Dict[str, Tuple[int, int, int, int]] = {}
        with self.fai_path.open("r", encoding="utf-8") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 5:
                    continue
                name = fields[0]
                self.index[name] = tuple(map(int, fields[1:5]))
        self.handle = self.fasta_path.open("rb")

    def close(self) -> None:
        self.handle.close()

    def __enter__(self) -> "IndexedFasta":
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        self.close()

    def resolve_contig(self, chrom: str) -> str:
        "Resolve common chromosome aliases to the actual FASTA contig name."
        candidates = [chrom]
        if chrom.startswith("chr"):
            candidates.append(chrom[3:])
        else:
            candidates.append(f"chr{chrom}")
        if chrom in {"chrM", "M"}:
            candidates.extend(["MT", "chrMT"])
        for candidate in candidates:
            if candidate in self.index:
                return candidate
        raise KeyError(f"No FASTA chromosome corresponding to {chrom!r}")

    def fetch(self, chrom: str, start0: int, end0: int) -> str:
        "Fetch an interval and clip positions outside chromosome boundaries."
        contig = self.resolve_contig(chrom)
        length, offset, line_bases, line_width = self.index[contig]
        start0 = max(0, int(start0))
        end0 = min(length, int(end0))
        if end0 <= start0:
            return ""

        pieces: List[bytes] = []
        position = start0
        while position < end0:
            line_index, in_line = divmod(position, line_bases)
            take = min(end0 - position, line_bases - in_line)
            byte_offset = offset + line_index * line_width + in_line
            self.handle.seek(byte_offset)
            pieces.append(self.handle.read(take))
            position += take
        return b"".join(pieces).decode("ascii").upper()


def open_text(path: Path):
    "Open plain or gzip-compressed GTF text according to its extension."
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def strip_version(identifier: str) -> str:
    "Remove terminal accession-version suffixes for table matching."
    return re.sub(r"\.\d+$", "", str(identifier or ""))


def parse_gtf_attributes(text: str) -> Dict[str, List[str]]:
    "Parse GTF attributes while retaining repeated entries such as tags."
    attributes: Dict[str, List[str]] = defaultdict(list)
    for item in text.strip().strip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        match = re.match(r"([^\s=]+)[\s=]+\"?([^\"]*)\"?$", item)
        if match:
            key, value = match.groups()
            attributes[key].append(value.strip())
    return dict(attributes)


def first_attr(attributes: Mapping[str, Sequence[str]], *keys: str) -> str:
    "Return the first available GTF attribute in candidate-field order."
    for key in keys:
        values = attributes.get(key)
        if values:
            return str(values[0])
    return ""


def merge_intervals(intervals: Iterable[Tuple[int, int]]) -> List[Tuple[int, int]]:
    "Merge overlapping intervals to avoid duplicate length contributions."
    ordered = sorted((int(start), int(end)) for start, end in intervals if end > start)
    merged: List[List[int]] = []
    for start, end in ordered:
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


def parse_mane_gtf(gtf_path: Path) -> Dict[str, TranscriptRecord]:
    "Parse MANE transcript records keyed by transcript identifier."
    if not gtf_path.exists():
        raise FileNotFoundError(f"GTF file not found: {gtf_path}")

    records: Dict[str, TranscriptRecord] = {}
    with open_text(gtf_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            chrom, _, feature, start, end, _, strand, _, attr_text = fields
            if feature not in {"transcript", "exon", "CDS"}:
                continue
            attributes = parse_gtf_attributes(attr_text)
            transcript_id = first_attr(attributes, "transcript_id")
            if not transcript_id:
                continue

            start0 = int(start) - 1
            end0 = int(end)
            record = records.setdefault(transcript_id, TranscriptRecord(transcript_id=transcript_id))
            record.gene_id = first_attr(attributes, "gene_id") or record.gene_id
            record.gene_name = first_attr(attributes, "gene_name", "gene") or record.gene_name
            record.chrom = chrom or record.chrom
            record.strand = strand or record.strand
            record.gene_biotype = first_attr(attributes, "gene_type", "gene_biotype") or record.gene_biotype
            record.transcript_biotype = first_attr(
                attributes, "transcript_type", "transcript_biotype"
            ) or record.transcript_biotype
            record.tags.extend(attributes.get("tag", []))

            for key in (
                "transcript_name",
                "refseq_transcript_id",
                "RefSeq",
                "havana_transcript",
            ):
                record.aliases.extend(attributes.get(key, []))

            if feature == "transcript":
                record.start0 = start0
                record.end0 = end0
            elif feature == "exon":
                record.exons.append((start0, end0))
            elif feature == "CDS":
                record.cds.append((start0, end0))

    for record in records.values():
        record.exons = merge_intervals(record.exons)
        record.cds = merge_intervals(record.cds)
        record.tags = sorted(set(record.tags))
        record.aliases = sorted(set(record.aliases))
        if record.start0 is None and record.exons:
            record.start0 = min(start for start, _ in record.exons)
        if record.end0 is None and record.exons:
            record.end0 = max(end for _, end in record.exons)

    return records


def transcript_selection_key(record: TranscriptRecord) -> Tuple[int, int, int, int]:
    "Rank candidate transcripts for one gene using the saved selection policy."
    mane_select = int(any("MANE_Select" in tag or "MANE Select" in tag for tag in record.tags))
    cds_length = interval_total_length(record.cds)
    spliced_length = interval_total_length(record.exons)
    return mane_select, cds_length, spliced_length, record.span_bp


def select_transcripts_by_gene(
    records: Mapping[str, TranscriptRecord],
    target_table: Optional[pd.DataFrame],
    gene_column: str,
    transcript_column: Optional[str],
) -> Tuple[List[str], Dict[str, TranscriptRecord], pd.DataFrame]:
    "Select MANE transcripts in the target PC1 table gene order."
    by_gene: Dict[str, List[TranscriptRecord]] = defaultdict(list)
    for record in records.values():
        keys = {record.gene_name, record.gene_id, strip_version(record.gene_id)}
        for key in keys:
            if key:
                by_gene[key].append(record)

    if target_table is None:
        gene_order = sorted({record.gene_name or strip_version(record.gene_id) for record in records.values()})
        target_transcripts: Dict[str, str] = {}
    else:
        gene_order = target_table[gene_column].astype(str).tolist()
        target_transcripts = {}
        if transcript_column:
            target_transcripts = dict(
                zip(
                    target_table[gene_column].astype(str),
                    target_table[transcript_column].fillna("").astype(str),
                )
            )

    selected: Dict[str, TranscriptRecord] = {}
    report_rows: List[Dict[str, object]] = []
    for gene in gene_order:
        candidates = list({record.transcript_id: record for record in by_gene.get(gene, [])}.values())
        requested = target_transcripts.get(gene, "")
        requested_no_version = strip_version(requested)

        exact_matches: List[TranscriptRecord] = []
        if requested_no_version:
            for candidate in candidates:
                candidate_ids = {
                    strip_version(candidate.transcript_id),
                    *(strip_version(alias) for alias in candidate.aliases),
                }
                if requested_no_version in candidate_ids:
                    exact_matches.append(candidate)

        pool = exact_matches or candidates
        chosen = max(pool, key=transcript_selection_key) if pool else None
        if chosen is not None:
            selected[gene] = chosen

        report_rows.append(
            {
                "gene": gene,
                "requested_transcript": requested,
                "selected_gtf_transcript": chosen.transcript_id if chosen else "",
                "candidate_count": len(candidates),
                "selection_method": (
                    "requested_transcript_match"
                    if exact_matches
                    else "MANE_Select_then_longest_CDS"
                    if chosen
                    else "not_found"
                ),
                "status": "matched" if chosen else "missing_annotation",
            }
        )

    return gene_order, selected, pd.DataFrame(report_rows)


def interval_total_length(intervals: Sequence[Tuple[int, int]]) -> int:
    "Calculate the total length of merged zero-based half-open intervals."
    return sum(end - start for start, end in merge_intervals(intervals))


def intervals_in_transcription_order(
    intervals: Sequence[Tuple[int, int]], strand: str
) -> List[Tuple[int, int]]:
    "Sort intervals in transcription order, reversing genomic order on the negative strand."
    return sorted(intervals, key=lambda item: item[0], reverse=(strand == "-"))


def reverse_complement(sequence: str) -> str:
    "Reverse-complement genomic sequence into transcription orientation."
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return sequence.translate(table)[::-1].upper()


def fetch_intervals_sequence(
    fasta: IndexedFasta,
    chrom: str,
    intervals: Sequence[Tuple[int, int]],
    strand: str,
) -> str:
    "Fetch and join intervals in transcription-oriented 5-prime to 3-prime order."
    merged = merge_intervals(intervals)
    sequence = "".join(fasta.fetch(chrom, start, end) for start, end in merged)
    return reverse_complement(sequence) if strand == "-" else sequence


def intersect_interval(
    interval: Tuple[int, int], boundary: Tuple[int, int]
) -> Optional[Tuple[int, int]]:
    "Intersect two half-open intervals, returning None when they do not overlap."
    start = max(interval[0], boundary[0])
    end = min(interval[1], boundary[1])
    return (start, end) if end > start else None


def infer_utr_intervals(record: TranscriptRecord) -> Tuple[List[Tuple[int, int]], List[Tuple[int, int]]]:
    "Infer 5-prime and 3-prime UTR intervals from exon and CDS boundaries."
    if not record.cds:
        return [], []
    cds_start = min(start for start, _ in record.cds)
    cds_end = max(end for _, end in record.cds)
    low_side: List[Tuple[int, int]] = []
    high_side: List[Tuple[int, int]] = []
    for exon in record.exons:
        low = intersect_interval(exon, (exon[0], min(exon[1], cds_start)))
        high = intersect_interval(exon, (max(exon[0], cds_end), exon[1]))
        if low:
            low_side.append(low)
        if high:
            high_side.append(high)
    if record.strand == "-":
        return merge_intervals(high_side), merge_intervals(low_side)
    return merge_intervals(low_side), merge_intervals(high_side)


def derive_introns(record: TranscriptRecord) -> List[Tuple[int, int]]:
    "Derive introns from gaps between adjacent exons."
    exons = sorted(merge_intervals(record.exons))
    return [
        (left[1], right[0])
        for left, right in zip(exons, exons[1:])
        if right[0] > left[1]
    ]


def safe_fraction(numerator: float, denominator: float) -> float:
    "Return a fraction or a missing value when its denominator is zero."
    return float(numerator / denominator) if denominator else float("nan")


def summarize_lengths(lengths: Sequence[int]) -> Tuple[float, float, float]:
    "Summarize mean, minimum and maximum interval lengths; return NA for empty inputs."
    if not lengths:
        return float("nan"), float("nan"), float("nan")
    return float(np.mean(lengths)), float(min(lengths)), float(max(lengths))


def build_annotation_and_structure_features(
    gene: str,
    record: TranscriptRecord,
    output_transcript: Optional[str] = None,
) -> Tuple[Dict[str, object], Dict[str, List[Tuple[int, int]]]]:
    "Build the annotation and gene-architecture feature set."
    exons = merge_intervals(record.exons)
    introns = derive_introns(record)
    cds = merge_intervals(record.cds)
    utr5, utr3 = infer_utr_intervals(record)

    exon_lengths = [end - start for start, end in exons]
    intron_lengths_genomic = [end - start for start, end in introns]
    exon_ordered = intervals_in_transcription_order(exons, record.strand)
    intron_ordered = intervals_in_transcription_order(introns, record.strand)

    span = record.span_bp
    exon_total = sum(exon_lengths)
    spliced_length = exon_total
    intron_total = sum(intron_lengths_genomic)
    cds_length = interval_total_length(cds)
    utr5_length = interval_total_length(utr5)
    utr3_length = interval_total_length(utr3)

    mean_exon, min_exon, max_exon = summarize_lengths(exon_lengths)
    mean_intron, min_intron, max_intron = summarize_lengths(intron_lengths_genomic)

    coding_exon_count = sum(
        1
        for exon in exons
        if any(intersect_interval(exon, cds_interval) is not None for cds_interval in cds)
    )
    orf_junction_count = max(coding_exon_count - 1, 0)

    if record.strand == "+":
        tss = int(record.start0)
        tes = int(record.end0) - 1
    else:
        tss = int(record.end0) - 1
        tes = int(record.start0)

    first_exon_length = (
        exon_ordered[0][1] - exon_ordered[0][0] if exon_ordered else float("nan")
    )
    last_exon_length = (
        exon_ordered[-1][1] - exon_ordered[-1][0] if exon_ordered else float("nan")
    )
    first_intron_length = (
        intron_ordered[0][1] - intron_ordered[0][0] if intron_ordered else float("nan")
    )
    last_intron_length = (
        intron_ordered[-1][1] - intron_ordered[-1][0] if intron_ordered else float("nan")
    )

    values: Dict[str, object] = {
        "gene": gene,
        "transcript": output_transcript or record.transcript_id,
        "chrom": record.chrom,
        "strand": record.strand,
        "tss_bed0": tss,
        "tes_bed0": tes,
        "gene_biotype": record.gene_biotype or record.transcript_biotype,
        "protein_coding": record.protein_coding,
        "transcript_span_bp": span,
        "log10_transcript_span_bp": math.log10(span + 0.1),
        "spliced_length_bp": spliced_length,
        "exon_total_length_bp": exon_total,
        "exon_count": len(exons),
        "mean_exon_length_bp": mean_exon,
        "min_exon_length_bp": min_exon,
        "max_exon_length_bp": max_exon,
        "first_exon_length_bp": first_exon_length,
        "last_exon_length_bp": last_exon_length,
        "intron_count": len(introns),
        "intron_total_length_bp": intron_total,
        "mean_intron_length_bp": mean_intron,
        "min_intron_length_bp": min_intron,
        "max_intron_length_bp": max_intron,
        "first_intron_length_bp": first_intron_length,
        "last_intron_length_bp": last_intron_length,
        "tss_to_first_splice_bp": first_exon_length if introns else float("nan"),
        "cds_length_bp": cds_length,
        "utr5_length_bp": utr5_length,
        "utr3_length_bp": utr3_length,
        "cds_fraction": safe_fraction(cds_length, spliced_length),
        "utr5_fraction": safe_fraction(utr5_length, spliced_length),
        "utr3_fraction": safe_fraction(utr3_length, spliced_length),
        "exon_fraction": safe_fraction(exon_total, span),
        "intron_fraction": safe_fraction(intron_total, span),
        "orf_exon_junction_count": orf_junction_count,
        "orf_exon_junction_density_per_kb": safe_fraction(
            orf_junction_count, cds_length / 1000.0
        ),
    }
    regions = {
        "exons": exons,
        "introns": introns,
        "cds": cds,
        "utr5": utr5,
        "utr3": utr3,
    }
    return values, regions


def relative_tss_window(
    tss_bed0: int,
    strand: str,
    relative_start: int,
    relative_end: int,
) -> Tuple[int, int]:
    "Convert transcription-oriented TSS windows to zero-based genomic intervals."
    if strand == "-":
        return tss_bed0 - relative_end + 1, tss_bed0 - relative_start + 1
    return tss_bed0 + relative_start, tss_bed0 + relative_end


def fetch_oriented_window(
    fasta: IndexedFasta,
    chrom: str,
    tss_bed0: int,
    strand: str,
    relative_start: int,
    relative_end: int,
) -> str:
    "Fetch a TSS-relative window in transcription orientation."
    start0, end0 = relative_tss_window(tss_bed0, strand, relative_start, relative_end)
    sequence = fasta.fetch(chrom, start0, end0)
    return reverse_complement(sequence) if strand == "-" else sequence


def sequence_composition(sequence: str) -> Tuple[float, float]:
    "Calculate GC and AT fractions using only valid A/C/G/T bases."
    counts = Counter(base for base in sequence.upper() if base in DNA_ALPHABET)
    total = sum(counts.values())
    if total == 0:
        return float("nan"), float("nan")
    return (counts["G"] + counts["C"]) / total, (counts["A"] + counts["T"]) / total


def cpg_features(sequence: str) -> Tuple[float, float]:
    "Calculate CpG counts and observed-to-expected ratios."
    sequence = sequence.upper()
    valid = "".join(base for base in sequence if base in DNA_ALPHABET)
    cpg_count = valid.count("CG")
    c_count = valid.count("C")
    g_count = valid.count("G")
    if not valid or c_count == 0 or g_count == 0:
        return float(cpg_count), float("nan")
    return float(cpg_count), cpg_count * len(valid) / (c_count * g_count)


def all_kmers(min_k: int = 1, max_k: int = 4) -> List[str]:
    "Enumerate k-mers by length and lexicographic order."
    return [
        "".join(parts)
        for k in range(min_k, max_k + 1)
        for parts in itertools.product(DNA_ALPHABET, repeat=k)
    ]


def kmer_frequencies(sequence: str, min_k: int = 1, max_k: int = 4) -> Dict[str, float]:
    "Calculate overlapping k-mer frequencies with a one-base sliding step."
    sequence = sequence.upper()
    output: Dict[str, float] = {}
    for k in range(min_k, max_k + 1):
        names = ["".join(parts) for parts in itertools.product(DNA_ALPHABET, repeat=k)]
        windows = [
            sequence[index : index + k]
            for index in range(max(0, len(sequence) - k + 1))
        ]
        valid = [window for window in windows if set(window) <= set(DNA_ALPHABET)]
        counts = Counter(valid)
        denominator = len(valid)
        for name in names:
            output[name] = safe_fraction(counts[name], denominator)
    return output


def codon_frequencies(cds_sequence: str) -> Dict[str, float]:
    "Calculate frequencies for 61 sense codons using width three and step three."
    sequence = cds_sequence.upper()
    codons = [sequence[index : index + 3] for index in range(0, len(sequence) - 2, 3)]
    valid = [
        codon
        for codon in codons
        if set(codon) <= set(DNA_ALPHABET) and codon not in STOP_CODONS
    ]
    counts = Counter(valid)
    return {codon: safe_fraction(counts[codon], len(valid)) for codon in SENSE_CODONS}


def build_sequence_features(
    record: TranscriptRecord,
    base_values: Mapping[str, object],
    regions: Mapping[str, List[Tuple[int, int]]],
    fasta: IndexedFasta,
    windows: Mapping[str, Tuple[int, int]],
    expand_kmers: bool,
) -> Dict[str, object]:
    "Build regional sequence, CpG, k-mer and codon features."
    tss = int(base_values["tss_bed0"])
    region_sequences = {
        name: fetch_oriented_window(
            fasta,
            record.chrom,
            tss,
            record.strand,
            relative_start,
            relative_end,
        )
        for name, (relative_start, relative_end) in windows.items()
    }
    body_start = windows["pause"][1]
    body_end = record.span_bp
    region_sequences["body"] = (
        fetch_oriented_window(
            fasta,
            record.chrom,
            tss,
            record.strand,
            body_start,
            body_end,
        )
        if body_end > body_start
        else ""
    )
    spliced_sequence = fetch_intervals_sequence(
        fasta, record.chrom, regions["exons"], record.strand
    )
    utr5_sequence = fetch_intervals_sequence(
        fasta, record.chrom, regions["utr5"], record.strand
    )
    cds_sequence = fetch_intervals_sequence(
        fasta, record.chrom, regions["cds"], record.strand
    )
    utr3_sequence = fetch_intervals_sequence(
        fasta, record.chrom, regions["utr3"], record.strand
    )

    output: Dict[str, object] = {}
    for name in ("promoter", "pause", "early_body"):
        gc_fraction, at_fraction = sequence_composition(region_sequences[name])
        cpg_count, cpg_oe = cpg_features(region_sequences[name])
        output[f"{name}_gc_fraction"] = gc_fraction
        output[f"{name}_at_fraction"] = at_fraction
        output[f"{name}_cpg_count"] = cpg_count
        output[f"{name}_cpg_oe"] = cpg_oe

    for prefix, sequence in (
        ("spliced", spliced_sequence),
        ("utr5", utr5_sequence),
        ("cds", cds_sequence),
        ("utr3", utr3_sequence),
    ):
        gc_fraction, _ = sequence_composition(sequence)
        output[f"{prefix}_gc_fraction"] = gc_fraction

    if expand_kmers:
        for region_name in ("pause", "body"):
            frequencies = kmer_frequencies(region_sequences[region_name], 1, 4)
            for kmer, value in frequencies.items():
                output[f"{region_name}_kmer_{kmer}_frequency"] = value
        for codon, value in codon_frequencies(cds_sequence).items():
            output[f"cds_codon_{codon}_frequency"] = value

    return output


def detect_column(frame: pd.DataFrame, candidates: Sequence[str], required: bool = True) -> Optional[str]:
    "Find common column names case-insensitively across target and PI tables."
    lookup = {str(column).lower(): str(column) for column in frame.columns}
    for candidate in candidates:
        if candidate.lower() in lookup:
            return lookup[candidate.lower()]
    if required:
        raise ValueError(f"Required column not found; candidates: {', '.join(candidates)}")
    return None


def load_target_table(path: Optional[Path]) -> Tuple[Optional[pd.DataFrame], Optional[str], Optional[str], Optional[str]]:
    "Load a gene-level PC1 target table and identify gene, transcript and response fields."
    if path is None:
        return None, None, None, None
    frame = pd.read_csv(path, sep="\t", low_memory=False)
    gene_column = detect_column(frame, ["gene", "gene_name", "symbol"])
    transcript_column = detect_column(
        frame, ["transcript", "transcript_id", "mane_transcript"], required=False
    )
    y_column = detect_column(frame, ["Y_PC1", "PC1_Y", "gene_PC1", "PC1"])
    frame = frame.drop_duplicates(subset=[gene_column], keep="first").reset_index(drop=True)
    return frame, gene_column, transcript_column, y_column


def parse_labelled_paths(items: Sequence[str]) -> List[Tuple[str, Path]]:
    "Parse repeatable LABEL=PATH arguments."
    output: List[Tuple[str, Path]] = []
    for item in items:
        if "=" not in item:
            raise ValueError(f"PI matrix argument must use LABEL=PATH: {item}")
        label, raw_path = item.split("=", 1)
        output.append((label.strip(), Path(raw_path).expanduser()))
    return output


def load_pi_features(
    labelled_paths: Sequence[Tuple[str, Path]], gene_order: Sequence[str]
) -> Tuple[pd.DataFrame, List[str]]:
    "Load PI matrices and align rows to the target gene order."
    merged = pd.DataFrame(index=pd.Index(gene_order, name="gene"))
    output_columns: List[str] = []
    for label, path in labelled_paths:
        if not path.exists():
            raise FileNotFoundError(f"PI matrix not found: {path}")
        frame = pd.read_csv(path, sep="\t", low_memory=False)
        gene_column = detect_column(
            frame, ["gene", "gene_name", "symbol"], required=False
        )
        combined_column = detect_column(
            frame, ["gene_transcript"], required=False
        )
        if gene_column is None and combined_column is None:
            raise ValueError(
                f"PI matrix {path} has neither gene nor gene_transcript columns"
            )


        if gene_column is None:
            gene_column = "__parsed_gene"
            frame[gene_column] = (
                frame[combined_column]
                .astype(str)
                .str.split("|", n=1, regex=False)
                .str[0]
            )
        else:
            frame[gene_column] = frame[gene_column].astype(str)

        frame = frame.drop_duplicates(subset=[gene_column], keep="first").set_index(gene_column)


        metadata_columns = {
            "gene",
            "gene_name",
            "symbol",
            "gene_transcript",
            "transcript",
            "transcript_id",
            "chrom",
            "strand",
            "tss_bed0",
            "tes_bed0",
        }
        value_columns = [
            column for column in frame.columns if str(column).lower() not in metadata_columns
        ]
        numeric = frame[value_columns].apply(pd.to_numeric, errors="coerce")
        numeric = numeric.dropna(axis=1, how="all")
        renamed = {
            column: f"PI_{label}_{str(column)}"
            for column in numeric.columns
        }
        numeric = numeric.rename(columns=renamed)
        numeric = numeric.reindex(gene_order)
        merged = merged.join(numeric, how="left")
        output_columns.extend(numeric.columns.tolist())
    return merged, output_columns


def minmax_scale(series: pd.Series) -> pd.Series:
    "Rescale observed values linearly to the interval from zero to one."
    numeric = pd.to_numeric(series, errors="coerce")
    minimum = numeric.min(skipna=True)
    maximum = numeric.max(skipna=True)
    if pd.isna(minimum) or pd.isna(maximum) or maximum == minimum:
        return pd.Series(np.nan, index=series.index, dtype=float)
    return (numeric - minimum) / (maximum - minimum)


def chromosome_numeric(series: pd.Series) -> pd.Series:
    "Encode standard chromosomes numerically and retain missing values for other contigs."
    mapping = {str(index): index for index in range(1, 23)}
    mapping.update({"X": 23, "Y": 24, "M": 25, "MT": 25})
    cleaned = series.astype(str).str.replace(r"^chr", "", regex=True)
    return cleaned.map(mapping).astype(float)


def build_scaled_architecture_model_view(frame: pd.DataFrame) -> pd.DataFrame:
    "Build a scaled architecture-feature view from the available predictors."
    output = pd.DataFrame(index=frame.index)
    output["gene"] = frame["gene"]


    output["tx.strand"] = frame["strand"].map({"+": 1, "-": 2, "*": 3})
    output["tx.gc.seq"] = frame["transcript_genomic_gc_fraction"]
    scaled_length = minmax_scale(frame["transcript_span_bp"])
    scaled_exon_count = minmax_scale(frame["exon_count"])
    output["tx.len"] = scaled_length
    output["tx.ex.count"] = scaled_exon_count


    scaled_ratio = scaled_length / scaled_exon_count.replace(0, np.nan)
    output["tx.ex.ratio"] = minmax_scale(scaled_ratio)
    output["tx.ex.width"] = minmax_scale(frame["mean_exon_length_bp"])
    output["tx.ex.seq"] = frame["exon_fraction"]


    genomic_start = pd.concat(
        [
            pd.to_numeric(frame["tss_bed0"], errors="coerce"),
            pd.to_numeric(frame["tes_bed0"], errors="coerce"),
        ],
        axis=1,
    ).min(axis=1)
    output["tx.loc"] = minmax_scale(genomic_start)
    output["tx.chr.loc"] = chromosome_numeric(frame["chrom"])
    return output


def ordered_feature_columns(
    frame: pd.DataFrame,
    pi_columns: Sequence[str],
    include_y: bool,
) -> List[str]:
    "Construct the final predictor order, expanding k-mer, codon and PI columns in place."
    columns = [column for column in BASE_FEATURE_ORDER if column in frame.columns]
    pause_kmers = [f"pause_kmer_{kmer}_frequency" for kmer in all_kmers(1, 4)]
    body_kmers = [f"body_kmer_{kmer}_frequency" for kmer in all_kmers(1, 4)]
    codons = [f"cds_codon_{codon}_frequency" for codon in SENSE_CODONS]
    columns.extend(column for column in pause_kmers if column in frame.columns)
    columns.extend(column for column in body_kmers if column in frame.columns)
    columns.extend(column for column in codons if column in frame.columns)
    columns.extend(column for column in pi_columns if column in frame.columns)
    if include_y and "Y_PC1" in frame.columns:
        columns.append("Y_PC1")
    return columns


def write_outputs(
    frame: pd.DataFrame,
    feature_columns: Sequence[str],
    output_prefix: Path,
    report: pd.DataFrame,
) -> None:
    "Write both matrix orientations, the scaled architecture view and annotation records."
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    gene_by_feature_path = Path(f"{output_prefix}.gene_by_feature.tsv")
    feature_by_gene_path = Path(f"{output_prefix}.feature_by_gene.tsv")
    scaled_architecture_path = Path(f"{output_prefix}.scaled_architecture_model_features.tsv")
    report_path = Path(f"{output_prefix}.annotation_report.tsv")

    gene_by_feature = frame[["gene", *feature_columns]].copy()
    gene_by_feature.to_csv(gene_by_feature_path, sep="\t", index=False, na_rep="NA")


    transposed = gene_by_feature.set_index("gene").T
    transposed.index.name = "feature_or_columns"
    transposed.to_csv(feature_by_gene_path, sep="\t", index=True, na_rep="NA")

    build_scaled_architecture_model_view(frame).to_csv(scaled_architecture_path, sep="\t", index=False, na_rep="NA")
    report.to_csv(report_path, sep="\t", index=False, na_rep="NA")

    print(f"Written: {gene_by_feature_path}")
    print(f"Written: {feature_by_gene_path}")
    print(f"Written: {scaled_architecture_path}")
    print(f"Written: {report_path}")


def build_argument_parser() -> argparse.ArgumentParser:
    "Define command-line options and explicit genomic-window boundaries."
    parser = argparse.ArgumentParser(
        description="Build gene features from MANE GTF and GRCh38 FASTA in the specified feature order."
    )
    parser.add_argument("--gtf", type=Path, required=True, help="Path to MANE GTF or GTF.gz")
    parser.add_argument("--fasta", type=Path, required=True, help="Path to GRCh38 FASTA")
    parser.add_argument("--fai", type=Path, default=None, help="Path to the FASTA .fai index; defaults to FASTA+.fai")
    parser.add_argument(
        "--pc1",
        type=Path,
        default=None,
        help="Gene PC1 table specifying gene order; PC1_Y is written as the final Y_PC1 column",
    )
    parser.add_argument(
        "--pi-matrix",
        action="append",
        default=[],
        metavar="LABEL=PATH",
        help="Repeatable PI matrix argument, for example GRO_seq=/path/gro.tsv",
    )
    parser.add_argument(
        "--output-prefix",
        type=Path,
        required=True,
        help="Output file prefix, for example results/POLII_XY",
    )
    parser.add_argument("--promoter-upstream", type=int, default=1000)
    parser.add_argument("--promoter-downstream", type=int, default=100)
    parser.add_argument("--pause-upstream", type=int, default=50)
    parser.add_argument("--pause-downstream", type=int, default=250)
    parser.add_argument("--early-body-start", type=int, default=250)
    parser.add_argument("--early-body-end", type=int, default=2250)
    parser.add_argument(
        "--no-expand-sequence-features",
        action="store_true",
        help="Skip expansion of the 680 pause/body 1-4-mer and 61 codon columns for a quick annotation-only run",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    "Run the feature-processing workflow and write its outputs."
    args = build_argument_parser().parse_args(argv)
    target_table, gene_column, transcript_column, y_column = load_target_table(args.pc1)

    print("Reading MANE GTF and constructing transcript regions...")
    records = parse_mane_gtf(args.gtf)
    if not records:
        raise RuntimeError("No transcripts were parsed from the GTF; check the format and transcript_id attributes")

    if target_table is None:
        gene_column = "gene"
    gene_order, selected, report = select_transcripts_by_gene(
        records,
        target_table,
        gene_column or "gene",
        transcript_column,
    )

    requested_transcripts: Dict[str, str] = {}
    y_values: Dict[str, float] = {}
    if target_table is not None and gene_column:
        if transcript_column:
            requested_transcripts = dict(
                zip(
                    target_table[gene_column].astype(str),
                    target_table[transcript_column].fillna("").astype(str),
                )
            )
        if y_column:
            y_values = dict(
                zip(
                    target_table[gene_column].astype(str),
                    pd.to_numeric(target_table[y_column], errors="coerce"),
                )
            )

    windows = {
        "promoter": (-args.promoter_upstream, args.promoter_downstream),
        "pause": (-args.pause_upstream, args.pause_downstream),
        "early_body": (args.early_body_start, args.early_body_end),
    }

    rows: List[Dict[str, object]] = []
    print(f"Computing {len(gene_order):,} target genes: structure and sequence features...")
    with IndexedFasta(args.fasta, args.fai) as fasta:
        for index, gene in enumerate(gene_order, start=1):
            record = selected.get(gene)
            if record is None:
                row: Dict[str, object] = {"gene": gene}
                if gene in requested_transcripts:
                    row["transcript"] = requested_transcripts[gene]
                if gene in y_values:
                    row["Y_PC1"] = y_values[gene]
                rows.append(row)
                continue

            output_transcript = requested_transcripts.get(gene) or record.transcript_id
            row, regions = build_annotation_and_structure_features(
                gene, record, output_transcript=output_transcript
            )
            row.update(
                build_sequence_features(
                    record,
                    row,
                    regions,
                    fasta,
                    windows,
                    expand_kmers=not args.no_expand_sequence_features,
                )
            )


            genomic_transcript_sequence = fetch_intervals_sequence(
                fasta,
                record.chrom,
                [(int(record.start0), int(record.end0))],
                record.strand,
            )
            row["transcript_genomic_gc_fraction"] = sequence_composition(
                genomic_transcript_sequence
            )[0]
            if gene in y_values:
                row["Y_PC1"] = y_values[gene]
            rows.append(row)

            if index % 1000 == 0:
                print(f"  Processed {index:,}/{len(gene_order):,} genes")

    frame = pd.DataFrame(rows)
    frame = frame.set_index("gene").reindex(gene_order).reset_index()

    labelled_pi_paths = parse_labelled_paths(args.pi_matrix)
    pi_frame, pi_columns = load_pi_features(labelled_pi_paths, gene_order)
    if pi_columns:
        frame = frame.join(pi_frame.reset_index(drop=True))

    feature_columns = ordered_feature_columns(
        frame,
        pi_columns,
        include_y=bool(y_values),
    )
    write_outputs(frame, feature_columns, args.output_prefix, report)

    matched = int((report["status"] == "matched").sum())
    print(f"Annotation matches: {matched:,}/{len(report):,} genes")
    if matched < len(report):
        print("Unmatched genes are retained with NA for unavailable features; see annotation_report.")
    print("PI sample columns are directly associated with Y_PC1 and must be excluded from predictors when training a Y_PC1 model.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"Error: {error}", file=sys.stderr)
        raise
