"""Coordinate-only DNA-overlap components and independently audited CV splits.

Coordinates are zero-based, half-open. Strand and window type do not exempt an
overlap. No targets or regulatory values are used to assign folds.
"""
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import random

REGIONS = ("tss_sequence", "first_donor_sequence", "first_acceptor_sequence", "pre_tes_sequence")


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path):
    with Path(path).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader), reader.fieldnames


def write_tsv(path, rows, fields):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def local_path(root, relative):
    root = Path(root).resolve()
    path = (root / relative).resolve()
    if not path.is_relative_to(root):
        raise ValueError(f"Path escapes experiment: {relative}")
    return path


def verify_source(source):
    source = Path(source).resolve()
    prepared = source / "data/bpnet_prepared"
    manifest_path = prepared / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    for name, expected in manifest["prepared_sha256"].items():
        path = (prepared / name).resolve()
        if not path.is_relative_to(source) or sha256(path) != expected:
            raise ValueError(f"Original data hash mismatch: {name}")
    return manifest, sha256(manifest_path)


def coordinate(value):
    if value in ("", "nan", "NaN", "None"):
        return None
    number = float(value)
    if not math.isfinite(number) or number < 0 or number != int(number):
        raise ValueError(f"Invalid genomic coordinate: {value}")
    return int(number)


def reconstruct_intervals(genes):
    intervals = []
    for index, row in enumerate(genes):
        present = set()
        for region in REGIONS:
            start = coordinate(row[f"{region}_clipped_start0"])
            end = coordinate(row[f"{region}_clipped_end0"])
            if start is None and end is None:
                if region not in ("first_donor_sequence", "first_acceptor_sequence") or row["splice_status"] != "intronless":
                    raise ValueError(f"Unexpected missing window: {row['gene']} {region}")
                continue
            if start is None or end is None or end <= start:
                raise ValueError(f"Invalid window: {row['gene']} {region}")
            present.add(region)
            intervals.append(dict(gene_index=index, gene=row["gene"], transcript=row["transcript"],
                                  chrom=row["chrom"], strand=row["strand"], region=region,
                                  start0=start, end0=end, fold=int(row["fold"])))
        expected = set(REGIONS) if row["splice_status"] != "intronless" else {"tss_sequence", "pre_tes_sequence"}
        if present != expected:
            raise ValueError(f"Window availability mismatch: {row['gene']}")
    return intervals


def check_original_intervals(intervals, original):
    fields = ("gene_index", "gene", "transcript", "chrom", "strand", "region", "start0", "end0", "fold")
    keys = lambda rows: sorted(tuple(str(row[field]) for field in fields) for row in rows)
    if keys(intervals) != keys(original):
        raise ValueError("Coordinates reconstructed from genes differ from original region_intervals.tsv")


class UnionFind:
    def __init__(self, size):
        self.parent = list(range(size))
        self.size = [1] * size

    def find(self, item):
        while self.parent[item] != item:
            self.parent[item] = self.parent[self.parent[item]]
            item = self.parent[item]
        return item

    def union(self, first, second):
        first, second = self.find(first), self.find(second)
        if first == second:
            return
        if self.size[first] < self.size[second]:
            first, second = second, first
        self.parent[second] = first
        self.size[first] += self.size[second]


def overlap_components(gene_count, intervals):
    union = UnionFind(gene_count)
    by_chrom = defaultdict(list)
    for interval in intervals:
        by_chrom[interval["chrom"]].append(interval)
    edges, same_gene_pairs, old_cross_pairs, old_cross_genes = [], 0, 0, set()
    for chrom in sorted(by_chrom):
        active = []
        for current in sorted(by_chrom[chrom], key=lambda row: (row["start0"], row["end0"], row["gene_index"], row["region"])):
            active = [row for row in active if row["end0"] > current["start0"]]
            for previous in active:
                overlap = min(previous["end0"], current["end0"]) - current["start0"]
                if overlap <= 0:
                    continue
                first, second = previous["gene_index"], current["gene_index"]
                if first == second:
                    same_gene_pairs += 1
                    continue
                union.union(first, second)
                if previous["fold"] != current["fold"]:
                    old_cross_pairs += 1
                    old_cross_genes.update((first, second))
                edges.append(dict(chrom=chrom, gene_index_a=first, gene_a=previous["gene"],
                                  region_a=previous["region"], gene_index_b=second, gene_b=current["gene"],
                                  region_b=current["region"], overlap_bp=overlap))
            active.append(current)
    groups = defaultdict(list)
    for index in range(gene_count):
        groups[union.find(index)].append(index)
    components = sorted(groups.values(), key=lambda rows: rows[0])
    summary = dict(same_gene_overlapping_window_pairs=same_gene_pairs,
                   original_cross_fold_overlap_pairs=old_cross_pairs,
                   original_cross_fold_overlap_genes=len(old_cross_genes),
                   inter_gene_overlapping_window_pairs=len(edges))
    return components, edges, summary


def assign_folds(components, chromosomes, seed, folds=10):
    if len(components) < folds:
        raise ValueError("Too few independent components for requested fold count")
    rng = random.Random(seed)
    order = list(range(len(components)))
    rng.shuffle(order)
    order.sort(key=lambda group: -len(components[group]))
    counts = [0] * folds
    chrom_counts = [Counter() for _ in range(folds)]
    gene_folds = [0] * len(chromosomes)
    gene_groups = [-1] * len(chromosomes)
    for group in order:
        members = components[group]
        group_chroms = Counter(chromosomes[index] for index in members)
        scores = [(counts[fold], sum(chrom_counts[fold][chrom] * n for chrom, n in group_chroms.items()))
                  for fold in range(folds)]
        best = min(scores)
        fold = rng.choice([i for i, score in enumerate(scores) if score == best])
        counts[fold] += len(members)
        chrom_counts[fold].update(group_chroms)
        for index in members:
            if gene_folds[index]:
                raise ValueError("Gene assigned more than once")
            gene_folds[index], gene_groups[index] = fold + 1, group
    if not all(gene_folds):
        raise ValueError("Gene omitted from fold assignment")
    return gene_folds, gene_groups


def merge_coverage(intervals, selected_folds):
    """Independent geometry audit: union coverage, not the graph's overlap edges."""
    selected_folds = set(selected_folds)
    by_chrom = defaultdict(list)
    for row in intervals:
        if int(row["fold"]) in selected_folds:
            by_chrom[row["chrom"]].append((int(row["start0"]), int(row["end0"])))
    result = {}
    for chrom, windows in by_chrom.items():
        merged = []
        for start, end in sorted(windows):
            if merged and start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            else:
                merged.append((start, end))
        result[chrom] = merged
    return result


def intersect_coverage(first, second):
    segments, bp = 0, 0
    for chrom in first.keys() & second.keys():
        a, b = first[chrom], second[chrom]
        i = j = 0
        while i < len(a) and j < len(b):
            overlap = min(a[i][1], b[j][1]) - max(a[i][0], b[j][0])
            if overlap > 0:
                segments += 1
                bp += overlap
            if a[i][1] <= b[j][1]:
                i += 1
            else:
                j += 1
    return segments, bp


def geometric_audits(intervals):
    coverage = {fold: merge_coverage(intervals, [fold]) for fold in range(1, 11)}
    fold_audit = []
    for first in range(1, 11):
        for second in range(first + 1, 11):
            segments, bp = intersect_coverage(coverage[first], coverage[second])
            fold_audit.append(dict(fold_a=first, fold_b=second, overlapping_segments=segments, overlapping_bp=bp))
    role_audit = []
    for test in range(1, 11):
        validation = test % 10 + 1
        roles = {"test": [test], "validation": [validation],
                 "train": [fold for fold in range(1, 11) if fold not in (test, validation)]}
        merged = {role: merge_coverage(intervals, folds) for role, folds in roles.items()}
        for first, second in (("train", "validation"), ("train", "test"), ("validation", "test")):
            segments, bp = intersect_coverage(merged[first], merged[second])
            role_audit.append(dict(round=test, role_a=first, role_b=second,
                                   overlapping_segments=segments, overlapping_bp=bp))
    if any(row["overlapping_bp"] for row in fold_audit + role_audit):
        raise ValueError("DNA overlap remains between folds or split roles")
    return fold_audit, role_audit


def verify_experiment(root, source, hash_source=True):
    root, source = Path(root).resolve(), Path(source).resolve()
    manifest_path = root / "split_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    if Path(manifest["source_root"]).resolve() != source:
        raise ValueError("Split manifest source differs from training source")
    if hash_source:
        source_manifest, source_hash = verify_source(source)
    else:
        source_manifest = json.loads((source / "data/bpnet_prepared/manifest.json").read_text())
        source_hash = sha256(source / "data/bpnet_prepared/manifest.json")
    if source_hash != manifest["source_manifest_sha256"]:
        raise ValueError("Original prepared manifest changed")
    if source_manifest["prepared_sha256"] != manifest["source_prepared_sha256"]:
        raise ValueError("Source artifact hash list changed")
    for relative, expected in manifest["generated_sha256"].items():
        if sha256(local_path(root, relative)) != expected:
            raise ValueError(f"Generated split hash mismatch: {relative}")
    genes, original_fields = read_tsv(source / "data/bpnet_prepared/genes.tsv")
    new, _ = read_tsv(root / "data/genes_repartitioned.tsv")
    assignments, _ = read_tsv(root / "data/bpnet_splits/fold_assignments.tsv")
    if len(new) != len(genes) or len(assignments) != len(genes) or len(genes) != manifest["gene_count"]:
        raise ValueError("Gene count changed")
    if len({row["gene"] for row in genes}) != len(genes):
        raise ValueError("Original gene IDs are not unique")
    group_folds = defaultdict(set)
    for index, (old, row, assignment) in enumerate(zip(genes, new, assignments)):
        if any(old[field] != row[field] for field in original_fields if field != "fold"):
            raise ValueError(f"Original metadata or row order changed: row {index}")
        if row["original_fold"] != old["fold"] or int(row["gene_index"]) != index:
            raise ValueError("Original fold/index preservation mismatch")
        if int(row["fold"]) not in range(1, 11):
            raise ValueError("Invalid fold")
        expected = dict(gene_index=str(index), gene=row["gene"], fold=row["fold"],
                        window_group=row["dna_overlap_group"], original_fold=old["fold"])
        if assignment != expected:
            raise ValueError("Assignments differ from training metadata")
        group_folds[row["dna_overlap_group"]].add(row["fold"])
    if any(len(folds) != 1 for folds in group_folds.values()):
        raise ValueError("Component split across folds")
    original_intervals = reconstruct_intervals(genes)
    original_table, _ = read_tsv(source / "data/bpnet_prepared/region_intervals.tsv")
    check_original_intervals(original_intervals, original_table)
    intervals, _ = read_tsv(root / "data/region_intervals_repartitioned.tsv")
    if len(intervals) != len(original_intervals):
        raise ValueError("Interval count changed")
    for old, row in zip(original_intervals, intervals):
        if any(str(old[key]) != row[key] for key in old if key != "fold"):
            raise ValueError("Coordinates or interval ordering changed")
        index = int(row["gene_index"])
        if row["fold"] != new[index]["fold"] or row["original_fold"] != str(old["fold"]) or row["dna_overlap_group"] != new[index]["dna_overlap_group"]:
            raise ValueError("Interval fold/group mismatch")
    fold_audit, role_audit = geometric_audits(intervals)
    for name, calculated in (("fold_pair_overlap_audit.tsv", fold_audit), ("window_overlap_audit.tsv", role_audit)):
        saved, _ = read_tsv(root / "data/bpnet_splits" / name)
        if saved != [{key: str(value) for key, value in row.items()} for row in calculated]:
            raise ValueError(f"Saved audit differs from independent geometric audit: {name}")
    for test in range(1, 11):
        rows, _ = read_tsv(root / f"data/bpnet_splits/round_{test:02d}.tsv")
        if len(rows) != len(genes):
            raise ValueError("Round missing/duplicate genes")
        for row, gene in zip(rows, new):
            role = "test" if int(gene["fold"]) == test else "validation" if int(gene["fold"]) == test % 10 + 1 else "train"
            if row != dict(gene=gene["gene"], window_group=gene["dna_overlap_group"], role=role):
                raise ValueError("Round roles differ from new metadata")
    counts = Counter(int(row["fold"]) for row in new)
    if set(counts) != set(range(1, 11)) or [counts[i] for i in range(1, 11)] != manifest["fold_gene_counts"]:
        raise ValueError("Fold count mismatch")
    return dict(status="PASS", gene_count=len(genes), interval_count=len(intervals),
                component_count=len(group_folds), fold_gene_counts=[counts[i] for i in range(1, 11)],
                cross_fold_overlap_pairs=0, cross_fold_overlap_bp=0,
                fold_pairs_checked=len(fold_audit), round_role_pairs_checked=len(role_audit),
                original_data_hashes_verified=hash_source, split_manifest_sha256=sha256(manifest_path),
                source_manifest_sha256=source_hash, metadata_sha256=sha256(root / "data/genes_repartitioned.tsv"))
