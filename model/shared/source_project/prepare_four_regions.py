"""Build new splice-centred inputs without changing the original project."""
import argparse
import json
import shutil
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from pyfaidx import Fasta

from corgi_polii.utils import sha256, split_indices, write_json
from inspect_splice_annotation import audit


REGIONS = {
    "tss_sequence": {"anchor": "tss_bed0", "anchor_type": "nucleotide",
                     "relative_start": -300, "relative_end": 300},
    "first_donor_sequence": {"anchor": "first_donor_boundary0", "anchor_type": "boundary",
                             "relative_start": -1000, "relative_end": 1000},
    "first_acceptor_sequence": {"anchor": "first_acceptor_boundary0", "anchor_type": "boundary",
                                "relative_start": -1000, "relative_end": 1000},
    "pre_tes_sequence": {"anchor": "tes_bed0", "anchor_type": "nucleotide",
                         "relative_start": -1000, "relative_end": 0},
}
LOOKUP = np.full(256, 4, dtype=np.uint8)
for index, base in enumerate(b"ACGT"):
    LOOKUP[base] = index
    LOOKUP[base + 32] = index
COMPLEMENT = np.array([3, 2, 1, 0, 4], dtype=np.uint8)
ONEHOT = np.concatenate([np.eye(4, dtype=np.bool_), np.zeros((1, 4), dtype=np.bool_)])


def boundary_window(fasta, chrom, boundary, strand, start=-1000, end=1000):
    """Half-open window around a boundary; junction lies between indices 999/1000."""
    length = end - start
    chrom_length = len(fasta[chrom])
    if not 0 <= boundary <= chrom_length or strand not in {"+", "-"} or length <= 0:
        raise ValueError("Invalid splice boundary")
    requested_start = boundary + start if strand == "+" else boundary - end
    requested_end = requested_start + length
    clipped_start, clipped_end = max(0, requested_start), min(chrom_length, requested_end)
    sequence = str(fasta[chrom][clipped_start:clipped_end]).upper()
    raw = np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)
    codes = np.full(length, 4, dtype=np.uint8)
    codes[clipped_start - requested_start:clipped_end - requested_start] = LOOKUP[raw]
    if strand == "-":
        codes = COMPLEMENT[codes[::-1]]
    return codes, {
        "requested_start0": requested_start, "requested_end0": requested_end,
        "clipped_start0": clipped_start, "clipped_end0": clipped_end,
        "left_pad_bp": max(0, -requested_start),
        "right_pad_bp": max(0, requested_end - chrom_length),
        "n_fraction": float(np.mean(codes == 4)),
    }


def decode(codes):
    return "".join("ACGTN"[int(code)] for code in codes)


def audit_overlap(intervals, genes):
    """Enumerate overlapping windows of different genes, independent of strand."""
    overlaps = []
    within_gene_pairs = 0
    for chrom, chunk in intervals.groupby("chrom", sort=True):
        active = []
        for row in chunk.sort_values(["start0", "end0", "gene_index"]).itertuples(index=False):
            active = [other for other in active if other.end0 > row.start0]
            for other in active:
                bp = min(row.end0, other.end0) - max(row.start0, other.start0)
                if bp <= 0:
                    continue
                if row.gene_index == other.gene_index:
                    within_gene_pairs += 1
                    continue
                if row.fold == other.fold:
                    continue
                overlaps.append({
                    "chrom": chrom, "gene_index_a": int(other.gene_index),
                    "gene_a": other.gene, "region_a": other.region, "fold_a": int(other.fold),
                    "gene_index_b": int(row.gene_index), "gene_b": row.gene,
                    "region_b": row.region, "fold_b": int(row.fold), "overlap_bp": int(bp),
                })
            active.append(row)
    columns = ["chrom", "gene_index_a", "gene_a", "region_a", "fold_a",
               "gene_index_b", "gene_b", "region_b", "fold_b", "overlap_bp"]
    pairs = pd.DataFrame(overlaps, columns=columns)
    summary = []
    for test_fold in range(1, 11):
        train, valid, test = split_indices(genes.fold.to_numpy(), test_fold)
        splits = {"train": set(train), "validation": set(valid), "test": set(test)}
        for left_role, right_role in [("train", "validation"), ("train", "test"),
                                      ("validation", "test")]:
            left, right = splits[left_role], splits[right_role]
            left_hits, right_hits, matching = set(), set(), 0
            for pair in pairs.itertuples(index=False):
                a, b = pair.gene_index_a, pair.gene_index_b
                if a in right and b in left:
                    a, b = b, a
                if a in left and b in right:
                    matching += 1
                    left_hits.add(a)
                    right_hits.add(b)
            summary.append({"test_fold": test_fold, "split_a": left_role, "split_b": right_role,
                            "overlap_window_pairs": matching, "genes_a": len(left_hits),
                            "genes_b": len(right_hits), "total_genes_a": len(left),
                            "total_genes_b": len(right)})
    return pairs, pd.DataFrame(summary), within_gene_pairs


def independent_sequence_qa(fasta, genes, prepared):
    """Independently fetch coordinates and reverse-complement as strings."""
    selected = []
    for strand in ("+", "-"):
        matching = genes.index[(genes.strand == strand) &
                               (genes.splice_status == "paired_first_intron")].to_numpy()
        selected.extend(matching[:64].tolist())
        selected.extend(matching[-16:].tolist())
    selected.extend(genes.index[genes.splice_status == "intronless"].tolist())
    selected = sorted(set(selected))
    complement = str.maketrans("ACGTN", "TGCAN")
    counts = Counter()
    for name in REGIONS:
        array = np.load(prepared / f"{name}.npy", mmap_mode="r")
        for index in selected:
            row = genes.iloc[index]
            if name.startswith("first_") and row.splice_status == "intronless":
                if array[index].any():
                    raise AssertionError("Intronless splice windows must be all zero")
                counts["intronless_zero_checks"] += 1
                continue
            start = int(row[f"{name}_requested_start0"])
            end = int(row[f"{name}_requested_end0"])
            size = len(fasta[row.chrom])
            text = str(fasta[row.chrom][max(0, start):min(size, end)]).upper()
            text = "N" * max(0, -start) + text + "N" * max(0, end - size)
            text = "".join(base if base in "ACGT" else "N" for base in text)
            if row.strand == "-":
                text = text.translate(complement)[::-1]
            data = array[index]
            actual = np.full(len(data), "N", dtype="<U1")
            known = data.any(axis=1)
            actual[known] = np.array(list("ACGT"))[data[known].argmax(axis=1)]
            if "".join(actual) != text:
                raise AssertionError(f"Independent reference mismatch: {row.gene} {name}")
            if name.startswith("first_"):
                boundary = int(row[REGIONS[name]["anchor"]])
                expected = boundary + REGIONS[name]["relative_start"] if row.strand == "+" else \
                    boundary - REGIONS[name]["relative_end"]
                if start != expected or end - start != 2000:
                    raise AssertionError("Splice window boundary off by one")
            counts[f"reference_checks_{row.strand}"] += 1
    return {"sampled_genes": len(selected), **dict(counts)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-root", required=True, type=Path)
    parser.add_argument("--gtf", required=True, type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    source = args.source_root.resolve()
    original = source / "data/bpnet_prepared"
    prepared = root / "data/bpnet_prepared"
    if root == source or source in root.parents and root.name == "data":
        raise ValueError("Output must be an isolated experiment directory")
    if (root / "data").is_symlink() or prepared.is_symlink():
        raise ValueError("New prepared directory must not alias original data")
    if prepared.exists() and any(prepared.iterdir()):
        raise FileExistsError("New prepared directory is nonempty; refusing to overwrite")
    prepared.mkdir(parents=True, exist_ok=True)
    old_manifest = json.loads((original / "manifest.json").read_text())
    for filename, expected in old_manifest["prepared_sha256"].items():
        if sha256(original / filename) != expected:
            raise ValueError(f"Original artifact hash changed: {filename}")
    genes, splice, annotation_audit = audit(original / "genes.tsv", args.gtf)
    if annotation_audit["failed_examples"]:
        write_json(root / "logs/annotation_failure.json", annotation_audit)
        raise ValueError("Representative transcript annotation mismatch")
    original_genes = genes.copy()
    genes = genes.loc[:, [column for column in genes if not column.startswith("early_5prime_sequence_")]]
    genes = pd.concat([genes, splice], axis=1)
    genes = genes.rename(columns={"window_group": "original_window_group"})
    genes["splice_junction_index0"] = 1000
    for name in ("feature_manifest.tsv", "module_sizes.json", "regulatory_log10_unscaled.npy",
                 "target.npy", "tss_sequence.npy", "pre_tes_sequence.npy"):
        shutil.copy2(original / name, prepared / name)
        if sha256(prepared / name) != old_manifest["prepared_sha256"][name]:
            raise AssertionError(f"Preserved artifact mismatch: {name}")
    shutil.copytree(source / "data/bpnet_splits", root / "data/bpnet_splits")
    for path in sorted((root / "data/bpnet_splits").glob("*.tsv")):
        if sha256(path) != sha256(source / "data/bpnet_splits" / path.name):
            raise AssertionError("Original fold file changed during copying")
    fasta_path = Path(next(item["path"] for item in old_manifest["source_records"]
                           if item["source"] == "fasta"))
    fasta = Fasta(str(fasta_path), as_raw=True, sequence_always_upper=True, rebuild=False)
    arrays = {name: np.lib.format.open_memmap(prepared / f"{name}.npy", mode="w+",
                                            dtype=np.bool_, shape=(len(genes), 2000, 4))
              for name in ("first_donor_sequence", "first_acceptor_sequence")}
    records, donor_motifs, acceptor_motifs = [], Counter(), Counter()
    for index, row in enumerate(genes.itertuples(index=False)):
        fields = {}
        for name in arrays:
            if row.splice_status == "intronless":
                arrays[name][index] = False
                info = {"requested_start0": None, "requested_end0": None,
                        "clipped_start0": None, "clipped_end0": None,
                        "left_pad_bp": 0, "right_pad_bp": 0, "n_fraction": 1.0}
            else:
                codes, info = boundary_window(fasta, row.chrom,
                    int(getattr(row, REGIONS[name]["anchor"])), row.strand)
                arrays[name][index] = ONEHOT[codes]
                if name == "first_donor_sequence":
                    donor_motifs[decode(codes[1000:1002])] += 1
                else:
                    acceptor_motifs[decode(codes[998:1000])] += 1
            fields.update({f"{name}_{key}": value for key, value in info.items()})
        records.append(fields)
        if (index + 1) % 3000 == 0:
            print(f"prepared splice sequences {index + 1}/{len(genes)}", flush=True)
    for array in arrays.values():
        array.flush()
    del arrays
    genes = pd.concat([genes, pd.DataFrame(records)], axis=1)
    interval_rows = []
    for index, row in enumerate(genes.itertuples(index=False)):
        for name in REGIONS:
            start = getattr(row, f"{name}_clipped_start0")
            end = getattr(row, f"{name}_clipped_end0")
            if pd.notna(start) and pd.notna(end) and end > start:
                interval_rows.append({"gene_index": index, "gene": row.gene,
                    "transcript": row.transcript, "chrom": row.chrom, "strand": row.strand,
                    "region": name, "start0": int(start), "end0": int(end), "fold": int(row.fold)})
    intervals = pd.DataFrame(interval_rows)
    pairs, overlap_summary, same_gene_pairs = audit_overlap(intervals, genes)
    intervals.to_csv(prepared / "region_intervals.tsv", sep="\t", index=False)
    pairs.to_csv(prepared / "cross_fold_overlap_pairs.tsv", sep="\t", index=False)
    overlap_summary.to_csv(prepared / "four_region_window_overlap_audit.tsv", sep="\t", index=False)
    genes.to_csv(prepared / "genes.tsv", sep="\t", index=False)
    invariant_columns = ["gene", "transcript", "chrom", "strand", "tss_bed0", "tes_bed0", "Y_PC1", "fold"]
    pd.testing.assert_frame_equal(genes[invariant_columns], original_genes[invariant_columns])
    validation_checks = independent_sequence_qa(fasta, genes, prepared)
    for name, spec in REGIONS.items():
        array = np.load(prepared / f"{name}.npy", mmap_mode="r")
        expected_shape = (len(genes), spec["relative_end"] - spec["relative_start"], 4)
        if array.shape != expected_shape or array.dtype != np.bool_:
            raise AssertionError(f"Invalid {name} array")
        for offset in range(0, len(genes), 512):
            if (array[offset:offset + 512].sum(axis=2) > 1).any():
                raise AssertionError("Invalid one-hot encoding")
        validation_checks[f"{name}_shape"] = list(array.shape)
    y = np.load(prepared / "target.npy")
    x = np.load(prepared / "regulatory_log10_unscaled.npy", mmap_mode="r")
    if x.shape != (len(genes), 2554) or not np.allclose(genes.Y_PC1, y, rtol=0, atol=1e-10):
        raise AssertionError("Original feature/target alignment lost")
    folds_hash = sha256(root / "data/bpnet_splits/fold_assignments.tsv")
    if folds_hash != old_manifest["prepared_sha256"]["../bpnet_splits/fold_assignments.tsv"]:
        raise AssertionError("Fold assignment hash differs")
    manifest = {
        "gene_count": len(genes), "feature_count": 2554, "module_sizes": old_manifest["module_sizes"],
        "region_lengths": {name: spec["relative_end"] - spec["relative_start"] for name, spec in REGIONS.items()},
        "region_definitions": REGIONS, "sequence_bp_length": 5600,
        "sequence_input_mode": "four_independent_dilated_encoders",
        "dna_sequences_concatenated": False, "separator_length": 0,
        "dna_encoding": old_manifest["dna_encoding"],
        "dna_orientation": old_manifest["dna_orientation"],
        "splice_definition": "exon_number 1 -> exon_number 2 of the original representative transcript",
        "splice_boundary_convention": "0-based interbase boundary; 1000 bases on either side; junction at index 1000",
        "intronless_policy": "keep original gene/target/fold; both splice windows all zero; status in metadata only",
        "annotation_status_counts": annotation_audit["status_counts"],
        "donor_dinucleotides": dict(donor_motifs), "acceptor_dinucleotides": dict(acceptor_motifs),
        "regulatory_transform": old_manifest["regulatory_transform"],
        "target_transform": old_manifest["target_transform"],
        "split_strategy": "inherited original folds, NOT regrouped for extended four-region windows",
        "split_seed": old_manifest["split_seed"],
        "fold_gene_counts": [int((genes.fold == fold).sum()) for fold in range(1, 11)],
        "cross_fold_overlap_pairs": len(pairs),
        "cross_fold_overlap_genes": len(set(pairs.gene_index_a) | set(pairs.gene_index_b)),
        "same_gene_overlapping_window_pairs": same_gene_pairs,
        "evaluation_scope": "exploratory inherited-fold CV; cross-fold DNA overlap present" if len(pairs)
                            else "four-region non-overlapping inherited gene CV; not homology-disjoint",
        "source_records": old_manifest["source_records"] + [
            {"source": "exon_gtf", "path": str(args.gtf), "sha256": sha256(args.gtf)},
            {"source": "original_prepared_manifest", "path": str(original / "manifest.json"),
             "sha256": sha256(original / "manifest.json")}],
        "fasta_stat": {"size": fasta_path.stat().st_size, "mtime_ns": fasta_path.stat().st_mtime_ns},
        "prepared_sha256": {path.name: sha256(path) for path in sorted(prepared.iterdir()) if path.is_file()},
    }
    for path in sorted((root / "data/bpnet_splits").glob("*.tsv")):
        manifest["prepared_sha256"][f"../bpnet_splits/{path.name}"] = sha256(path)
    write_json(prepared / "manifest.json", manifest)
    validation = {"status": "PASS", "manifest_sha256": sha256(prepared / "manifest.json"),
                  "annotation_audit": annotation_audit, "sequence_checks": validation_checks,
                  "original_data_unchanged": all(sha256(original / name) == expected
                        for name, expected in old_manifest["prepared_sha256"].items()),
                  "cross_fold_overlap_pairs": len(pairs),
                  "unbiased_cv_eligible_with_original_folds": len(pairs) == 0}
    if not validation["original_data_unchanged"]:
        raise AssertionError("Original source data changed")
    write_json(root / "logs/bpnet_validation.json", validation)
    write_json(root / "logs/four_region_preparation_summary.json", {
        "status": "PASS", "genes": len(genes), "features": 2554,
        "annotation_status_counts": annotation_audit["status_counts"],
        "regions": manifest["region_lengths"], "sequence_input_mode": manifest["sequence_input_mode"],
        "donor_dinucleotides": dict(donor_motifs), "acceptor_dinucleotides": dict(acceptor_motifs),
        "cross_fold_overlap_pairs": len(pairs), "cross_fold_overlap_genes": manifest["cross_fold_overlap_genes"],
        "sequence_checks": validation_checks, "original_data_unchanged": True})
    fasta.close()
    print(json.dumps({"status": "PASS", "annotation": annotation_audit["status_counts"],
                      "regions": manifest["region_lengths"], "cross_fold_overlap_pairs": len(pairs),
                      "cross_fold_overlap_genes": manifest["cross_fold_overlap_genes"],
                      "donor_dinucleotides": dict(donor_motifs),
                      "acceptor_dinucleotides": dict(acceptor_motifs)}, indent=2), flush=True)


if __name__ == "__main__":
    main()
