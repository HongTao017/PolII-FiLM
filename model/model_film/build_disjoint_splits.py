"""Create a new genomic-coordinate-disjoint 10-fold overlay; never edit inputs."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path

from disjoint_splits import (assign_folds, check_original_intervals, geometric_audits,
                             overlap_components, read_tsv, reconstruct_intervals,
                             sha256, verify_experiment, verify_source, write_tsv)


def build(root, source, seed):
    root, source = Path(root).resolve(), Path(source).resolve()
    if root == source or root.is_relative_to(source) or source.is_relative_to(root):
        raise ValueError("Output experiment must be independent of original source")
    if (root / "split_manifest.json").exists() or (root / "data").exists():
        raise FileExistsError("Refusing to overwrite an existing split/data directory")
    original_manifest, original_hash = verify_source(source)
    genes, fields = read_tsv(source / "data/sequence_prepared/genes.tsv")
    if len(genes) != original_manifest["gene_count"] or len({row["gene"] for row in genes}) != len(genes):
        raise ValueError("Original gene counts or uniqueness invalid")
    if set(fields) & {"original_fold", "gene_index", "dna_overlap_group"}:
        raise ValueError("Reserved overlay column already present")
    intervals = reconstruct_intervals(genes)
    saved_intervals, _ = read_tsv(source / "data/sequence_prepared/region_intervals.tsv")
    check_original_intervals(intervals, saved_intervals)
    components, edges, old_audit = overlap_components(len(genes), intervals)
    if old_audit["original_cross_fold_overlap_pairs"] != original_manifest["cross_fold_overlap_pairs"] or old_audit["original_cross_fold_overlap_genes"] != original_manifest["cross_fold_overlap_genes"]:
        raise ValueError("Recomputed original overlap differs from source manifest")
    folds, groups = assign_folds(components, [row["chrom"] for row in genes], seed)
    new_genes = [{**row, "fold": folds[i], "original_fold": row["fold"], "gene_index": i,
                  "dna_overlap_group": groups[i]} for i, row in enumerate(genes)]
    new_intervals = [{**row, "original_fold": row["fold"], "fold": folds[row["gene_index"]],
                      "dna_overlap_group": groups[row["gene_index"]]} for row in intervals]
    fold_audit, role_audit = geometric_audits(new_intervals)
    for edge in edges:
        if folds[edge["gene_index_a"]] != folds[edge["gene_index_b"]]:
            raise ValueError("Overlap edge crosses folds")
    root.mkdir(parents=True, exist_ok=True)
    write_tsv(root / "data/genes_repartitioned.tsv", new_genes, fields + ["original_fold", "gene_index", "dna_overlap_group"])
    write_tsv(root / "data/region_intervals_repartitioned.tsv", new_intervals,
              list(intervals[0]) + ["original_fold", "dna_overlap_group"])
    split_dir = root / "data/gene_splits"
    assignment = [dict(gene_index=i, gene=row["gene"], fold=folds[i], window_group=groups[i], original_fold=row["fold"])
                  for i, row in enumerate(genes)]
    write_tsv(split_dir / "fold_assignments.tsv", assignment, list(assignment[0]))
    summaries = []
    for test in range(1, 11):
        validation = test % 10 + 1
        round_rows = [dict(gene=row["gene"], window_group=groups[i],
                           role="test" if folds[i] == test else "validation" if folds[i] == validation else "train")
                      for i, row in enumerate(genes)]
        write_tsv(split_dir / f"round_{test:02d}.tsv", round_rows, ["gene", "window_group", "role"])
        counts = Counter(row["role"] for row in round_rows)
        summaries.append(dict(round=test, test_fold=test, validation_fold=validation,
                              train_genes=counts["train"], validation_genes=counts["validation"], test_genes=counts["test"]))
    write_tsv(split_dir / "split_summary.tsv", summaries, list(summaries[0]))
    group_rows = [dict(window_group=group, gene_count=len(members), fold=folds[members[0]],
                       chrom=genes[members[0]]["chrom"], genes=";".join(genes[i]["gene"] for i in members))
                  for group, members in enumerate(components)]
    write_tsv(split_dir / "overlap_group_sizes.tsv", group_rows, list(group_rows[0]))
    write_tsv(split_dir / "overlap_component_edges.tsv", edges,
              ["chrom", "gene_index_a", "gene_a", "region_a", "gene_index_b", "gene_b", "region_b", "overlap_bp"])
    write_tsv(split_dir / "fold_pair_overlap_audit.tsv", fold_audit, list(fold_audit[0]))
    write_tsv(split_dir / "window_overlap_audit.tsv", role_audit, list(role_audit[0]))
    write_tsv(split_dir / "cross_fold_overlap_pairs.tsv", [],
              ["chrom", "gene_index_a", "gene_a", "region_a", "gene_index_b", "gene_b", "region_b", "overlap_bp"])
    counts = Counter(folds)
    manifest = dict(created_utc=datetime.now(timezone.utc).isoformat(), source_root=str(source),
                    source_manifest_sha256=original_hash, source_prepared_sha256=original_manifest["prepared_sha256"],
                    gene_count=len(genes), feature_count=original_manifest["feature_count"],
                    region_lengths=original_manifest["region_lengths"], interval_count=len(intervals), split_seed=seed,
                    split_strategy="All positive genomic overlaps across any of four windows, both strands; transitive connected components; largest-first size balance and chromosome tie balance",
                    coordinate_convention="0-based half-open; touching endpoints are not overlaps",
                    targets_used_for_fold_assignment=False, original_gene_order_preserved=True,
                    original_arrays_reused_read_only=True, original_metadata_preserved_except_fold_in_new_overlay=True,
                    original_files_modified=False, gene_count_deleted=0,
                    component_count=len(components), multi_gene_components=sum(len(group) > 1 for group in components),
                    genes_in_multi_gene_components=sum(len(group) for group in components if len(group) > 1),
                    largest_component_gene_count=max(map(len, components)),
                    fold_gene_counts=[counts[i] for i in range(1, 11)],
                    genes_changing_fold=sum(folds[i] != int(row["fold"]) for i, row in enumerate(genes)),
                    **old_audit, cross_fold_overlap_pairs=0, cross_fold_overlap_genes=0,
                    cross_fold_overlap_bp=0, fold_pairs_checked=45, round_role_pairs_checked=30,
                    limitation="Input genomic loci are disjoint across folds; homologous DNA at distinct loci and all other forms of leakage are not automatically controlled",
                    generated_sha256={str(path.relative_to(root)): sha256(path) for path in sorted((root / "data").rglob("*.tsv"))})
    with (root / "split_manifest.json").open("x") as handle:
        json.dump(manifest, handle, indent=2)
    # Rehash the source after writing outputs and redo the independent geometric audit.
    report = verify_experiment(root, source, hash_source=True)
    with (root / "split_validation.json").open("x") as handle:
        json.dump(report, handle, indent=2)
    return {**report, **{key: manifest[key] for key in ("multi_gene_components", "genes_in_multi_gene_components",
                                                      "largest_component_gene_count", "genes_changing_fold",
                                                      "original_cross_fold_overlap_pairs", "original_cross_fold_overlap_genes")}}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--source", type=Path)
    parser.add_argument("--seed", type=int, default=20260916)
    args = parser.parse_args()
    source = args.source or Path(json.loads((args.root / "training_config.json").read_text())["source_root"])
    print(json.dumps(build(args.root, source, args.seed), indent=2))
