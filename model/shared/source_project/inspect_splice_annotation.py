"""Read-only audit: match existing representative transcripts to MANE exons."""
import argparse
import gzip
import json
import re
from collections import Counter
from pathlib import Path

import pandas as pd


ATTR = re.compile(r'(\w+) "([^"]*)";')


def read_exons(gtf_path, transcripts):
    exons = {}
    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    with opener(gtf_path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9 or parts[2] != "exon":
                continue
            attr = dict(ATTR.findall(parts[8]))
            transcript = attr.get("transcript_id")
            if transcript not in transcripts:
                continue
            record = (parts[0], parts[6], int(parts[3]) - 1, int(parts[4]),
                      int(attr["exon_number"]))
            exons.setdefault(transcript, []).append(record)
    return exons


def describe_transcript(row, exons):
    records = sorted(set(exons.get(row.transcript, [])), key=lambda item: item[4])
    result = {"exon_count": len(records)}
    if not records:
        return {**result, "splice_status": "missing_transcript_annotation"}
    if any(chrom != row.chrom or strand != row.strand
           for chrom, strand, _, _, _ in records):
        return {**result, "splice_status": "chromosome_or_strand_mismatch"}
    if [record[4] for record in records] != list(range(1, len(records) + 1)):
        return {**result, "splice_status": "invalid_exon_numbering"}
    coordinates = [(start, end) for _, _, start, end, _ in records]
    if any(start >= end for start, end in coordinates):
        return {**result, "splice_status": "invalid_exon_length"}
    oriented = sorted(coordinates, reverse=row.strand == "-")
    if coordinates != oriented:
        return {**result, "splice_status": "exon_order_mismatch"}
    tss = records[0][2] if row.strand == "+" else records[0][3] - 1
    tes = records[-1][3] - 1 if row.strand == "+" else records[-1][2]
    if tss != row.tss_bed0 or tes != row.tes_bed0:
        return {**result, "splice_status": "transcript_endpoint_mismatch",
                "annotation_tss_bed0": tss, "annotation_tes_bed0": tes}
    if len(records) == 1:
        return {**result, "splice_status": "intronless"}
    if any(left[1] >= right[0] for left, right in zip(
            sorted(coordinates)[:-1], sorted(coordinates)[1:])):
        return {**result, "splice_status": "overlapping_or_adjacent_exons"}
    first, second = records[:2]
    if row.strand == "+":
        donor, acceptor = first[3], second[2]
        intron_start, intron_end = donor, acceptor
    else:
        donor, acceptor = first[2], second[3]
        intron_start, intron_end = acceptor, donor
    if intron_end <= intron_start:
        return {**result, "splice_status": "invalid_first_intron"}
    return {**result, "splice_status": "paired_first_intron",
            "first_donor_boundary0": donor, "first_acceptor_boundary0": acceptor,
            "first_intron_start0": intron_start, "first_intron_end0": intron_end,
            "first_intron_length": intron_end - intron_start}


def audit(genes_path, gtf_path):
    genes = pd.read_csv(genes_path, sep="\t", dtype={"gene": str, "transcript": str})
    exons = read_exons(gtf_path, set(genes.transcript))
    rows = [describe_transcript(row, exons) for row in genes.itertuples(index=False)]
    details = pd.DataFrame(rows)
    statuses = Counter(details.splice_status)
    failures = ~details.splice_status.isin(["paired_first_intron", "intronless"])
    summary = {"genes": len(genes), "status_counts": dict(statuses),
               "failed_examples": genes.loc[failures, ["gene", "transcript"]]
                    .join(details.loc[failures]).head(10).to_dict("records"),
               "genes_path": str(genes_path), "gtf_path": str(gtf_path)}
    return genes, details, summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genes", required=True, type=Path)
    parser.add_argument("--gtf", required=True, type=Path)
    args = parser.parse_args()
    _, _, summary = audit(args.genes, args.gtf)
    print(json.dumps(summary, indent=2))
    if summary["failed_examples"]:
        raise SystemExit("Annotation audit failed: do not fabricate splice coordinates")


if __name__ == "__main__":
    main()
