import csv
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
from pyfaidx import Fasta

from .utils import sha256, split_indices, write_json

DNA_LOOKUP = np.full(256, 4, dtype=np.uint8)
for _index, _base in enumerate(b"ACGT"):
    DNA_LOOKUP[_base] = _index
    DNA_LOOKUP[_base + 32] = _index
COMPLEMENT = np.array([3, 2, 1, 0, 4], dtype=np.uint8)
ANNOTATION_COLUMNS = ["gene", "transcript", "chrom", "strand", "tss_bed0", "tes_bed0"]


def read_table(path, usecols=None):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        header = next(csv.reader(handle, delimiter="\t"))
    if len(header) != len(set(header)):
        raise ValueError(f"Duplicate header: {path}")
    table = pd.read_csv(path, sep="\t", dtype={"gene": str}, usecols=usecols)
    if "gene" not in table or table.gene.isna().any() or table.gene.duplicated().any():
        raise ValueError(f"Invalid gene IDs: {path}")
    return table.set_index("gene")


def align(table, genes, source):
    missing = genes.difference(table.index)
    if len(missing):
        raise ValueError(f"{source}: {len(missing)} missing genes; first: {missing[:5].tolist()}")
    return table.loc[genes]


def dna_window(fasta, chrom, tss, strand, length, relative_start=-300, relative_end=2000):
    if relative_end - relative_start != length or not relative_start <= 0 < relative_end:
        raise ValueError("Invalid transcription-relative window")
    if strand not in {"+", "-"} or chrom not in fasta:
        raise ValueError(f"Invalid reference coordinate: {chrom} {strand}")
    chrom_length = len(fasta[chrom])
    if not 0 <= tss < chrom_length:
        raise ValueError(f"TSS outside chromosome: {chrom}:{tss}")
    start = tss + relative_start if strand == "+" else tss - relative_end + 1
    end = start + length
    clipped_start, clipped_end = max(0, start), min(chrom_length, end)
    text = fasta[chrom][clipped_start:clipped_end].seq
    codes = np.full(length, 4, dtype=np.uint8)
    raw = np.frombuffer(text.encode("ascii"), dtype=np.uint8)
    codes[clipped_start - start:clipped_end - start] = DNA_LOOKUP[raw]
    if strand == "-":
        codes = COMPLEMENT[codes[::-1]]
    record = {
        "window_start0_requested": start, "window_end0_requested": end,
        "window_start0_clipped": clipped_start, "window_end0_clipped": clipped_end,
        "reference_left_pad_bp": max(0, -start),
        "reference_right_pad_bp": max(0, end - chrom_length),
        "tss_sequence_index0": -relative_start,
        "n_fraction": float(np.mean(codes == 4)),
    }
    return codes, record


def assign_overlap_folds(metadata, seed):
    """Keep connected overlapping intervals together; balance without looking at Y."""
    groups = np.full(len(metadata), -1, dtype=np.int64)
    group = -1
    for _, chunk in metadata.assign(row_index=np.arange(len(metadata))).groupby("chrom", sort=True):
        current_end = -1
        for row in chunk.sort_values(["window_start0_clipped", "window_end0_clipped"]).itertuples():
            if row.window_start0_clipped >= current_end:
                group += 1
            groups[row.row_index] = group
            current_end = max(current_end, row.window_end0_clipped)
    unique, counts = np.unique(groups, return_counts=True)
    if len(unique) < 10:
        raise ValueError("Fewer than 10 independent window groups; choose another evaluation design")
    rng = np.random.default_rng(seed)
    tie_order = rng.permutation(len(unique))
    order = tie_order[np.argsort(-counts[tie_order], kind="stable")]
    sizes = np.zeros(10, dtype=np.int64)
    assignments = np.zeros(len(unique), dtype=np.int64)
    for index in order:
        smallest = np.flatnonzero(sizes == sizes.min())
        fold = int(rng.choice(smallest))
        assignments[index] = fold + 1
        sizes[fold] += counts[index]
    return assignments[groups], groups


def audit_window_overlap(genes):
    """Counts test windows touching another split, not independent sequence blocks."""
    output = []
    for fold in range(1, 11):
        train_idx, valid_idx, test_idx = split_indices(genes.fold.values, fold)
        for role, reference_idx in (("train", train_idx), ("validation", valid_idx)):
            ref = genes.iloc[reference_idx]
            query = genes.iloc[test_idx]
            hit_count, best_bp = 0, []
            for chrom, chunk in query.groupby("chrom"):
                ref_chr = ref.loc[ref.chrom == chrom]
                starts = ref_chr.window_start0_clipped.to_numpy()
                ends = ref_chr.window_end0_clipped.to_numpy()
                for row in chunk.itertuples():
                    overlap = np.maximum(0, np.minimum(ends, row.window_end0_clipped)
                                         - np.maximum(starts, row.window_start0_clipped))
                    best = int(overlap.max()) if len(overlap) else 0
                    hit_count += int(best > 0)
                    best_bp.append(best)
            output.append({"test_fold": fold, "reference_split": role,
                           "test_gene_count": len(query), "overlapping_test_genes": hit_count,
                           "overlapping_fraction": hit_count / len(query),
                           "mean_max_overlap_bp": float(np.mean(best_bp))})
    return pd.DataFrame(output)


def prepare_project(project, sources_path):
    project, sources_path = Path(project).resolve(), Path(sources_path).resolve()
    cfg = json.loads(sources_path.read_text(encoding="utf-8"))
    prepared, inputs = project / "data/prepared", project / "data/input"
    prepared.mkdir(parents=True, exist_ok=True)
    inputs.mkdir(parents=True, exist_ok=True)
    if (prepared / "manifest.json").exists():
        raise FileExistsError("Prepared data already complete; use validate.py instead of overwriting")
    annotation_with_target = read_table(cfg["annotation"], ANNOTATION_COLUMNS + ["Y_PC1"])
    genes = annotation_with_target.index
    if len(genes) != cfg["expected_genes"]:
        raise ValueError(f"Expected {cfg['expected_genes']} genes, got {len(genes)}")
    target = pd.to_numeric(annotation_with_target.Y_PC1, errors="raise").to_numpy(dtype=np.float64)
    if not np.isfinite(target).all():
        raise ValueError("Nonfinite Y_PC1")
    parts = []
    modules = {}
    feature_records = []
    provenance = []
    if cfg["split_strategy"] != "overlap_grouped":
        raise ValueError("This project requires overlap-grouped folds")
    for name in ("annotation",):
        path = Path(cfg[name])
        provenance.append({"source": name, "path": str(path), "sha256": sha256(path)})
        annotation_with_target.reset_index().to_csv(inputs / "annotation.tsv", sep="\t", index=False)
    for module, spec in cfg["assays"].items():
        path = Path(spec["path"])
        table = align(read_table(path), genes, path)
        if "Y_PC1" not in table or not np.allclose(table.Y_PC1, target, rtol=0, atol=1e-10):
            raise ValueError(f"Target mismatch in {module}")
        columns = [c for c in table.columns if c != "Y_PC1"]
        if len(columns) != spec["features"]:
            raise ValueError(f"Feature count mismatch in {module}: {len(columns)}")
        if any("Y_PC1" in c or not ("_gene_body_" in c or "_pause_" in c) for c in columns):
            raise ValueError(f"Unexpected feature name in {module}")
        values = table[columns].to_numpy(dtype=np.float64)
        finite = np.isfinite(values)
        if np.isinf(values).any() or np.any(values[finite] < 0) or np.any(values[finite] % 1 != 0):
            raise ValueError(f"Invalid counts in {module}")
        parts.append(table[columns])
        modules[module] = len(columns)
        feature_records.extend({"module": module, "feature": f, "transform": "log10(count+1)"}
                               for f in columns)
        table.reset_index().to_csv(inputs / f"{module}_{len(columns)}X_with_Y_PC1.tsv",
                                  sep="\t", index=False)
        provenance.append({"source": module, "path": str(path), "sha256": sha256(path)})
    raw = pd.concat(parts, axis=1)
    expected_features = sum(spec["features"] for spec in cfg["assays"].values())
    if raw.columns.duplicated().any() or len(raw.columns) != expected_features:
        raise ValueError(f"Expected {expected_features} unique regulatory features")
    x = raw.to_numpy(dtype=np.float64, copy=True)
    if np.isinf(x).any():
        raise ValueError("Infinite raw features")
    x = np.log10(x + 1)
    np.save(prepared / "features_unscaled.npy", x.astype(np.float32))
    np.save(prepared / "target.npy", target)
    combined = pd.concat([raw, pd.Series(target, index=genes, name="Y_PC1")], axis=1)
    combined.reset_index().to_csv(inputs / f"regulatory_counts_{expected_features}X_with_Y_PC1.tsv",
                                  sep="\t", index=False)
    pd.DataFrame(feature_records).assign(feature_index=np.arange(len(feature_records))).to_csv(
        prepared / "feature_manifest.tsv", sep="\t", index=False)
    write_json(prepared / "module_sizes.json", modules)
    annotation = annotation_with_target.drop(columns="Y_PC1")
    if annotation.isna().any().any():
        raise ValueError("Missing annotation")
    for name in ("tss_bed0", "tes_bed0"):
        values = pd.to_numeric(annotation[name], errors="raise")
        if not np.isfinite(values).all() or np.any(values % 1 != 0):
            raise ValueError(f"Invalid coordinates: {name}")
        annotation[name] = values.astype(np.int64)
    metadata = annotation.assign(Y_PC1=target).reset_index()
    (project / "data/splits").mkdir(parents=True, exist_ok=True)
    fasta_path = Path(cfg["fasta"])
    if not Path(str(fasta_path) + ".fai").is_file():
        raise FileNotFoundError("Reference FASTA index is required; source reference is read-only")
    (project / "data/reference").mkdir(parents=True, exist_ok=True)
    for source in (fasta_path, Path(str(fasta_path) + ".fai")):
        destination = project / "data/reference" / source.name
        if not destination.exists():
            destination.symlink_to(source)
    print(f"Reference checksum: {fasta_path}", flush=True)
    provenance.append({"source": "fasta", "path": str(fasta_path), "sha256": sha256(fasta_path)})
    dna_path = prepared / "dna_onehot.npy"
    temporary = prepared / "dna_onehot.partial.npy"
    dna = np.lib.format.open_memmap(temporary, mode="w+", dtype=np.bool_,
                                    shape=(len(metadata), cfg["sequence_length"], 4))
    onehot = np.concatenate([np.eye(4, dtype=np.bool_), np.zeros((1, 4), dtype=np.bool_)])
    window_records = []
    with Fasta(str(fasta_path), build_index=False, rebuild=False) as fasta:
        for index, row in enumerate(metadata.itertuples(index=False)):
            codes, record = dna_window(fasta, row.chrom, int(row.tss_bed0), row.strand,
                                       cfg["sequence_length"], cfg["window_start_relative"],
                                       cfg["window_end_relative"])
            dna[index] = onehot[codes]
            window_records.append(record)
            if index % 500 == 0 or index + 1 == len(metadata):
                print(f"DNA extracted {index + 1}/{len(metadata)}", flush=True)
    dna.flush()
    del dna
    temporary.replace(dna_path)
    metadata = pd.concat([metadata, pd.DataFrame(window_records)], axis=1)
    metadata["fold"], metadata["window_group"] = assign_overlap_folds(metadata, cfg["split_seed"])
    metadata[["gene", "fold", "window_group"]].to_csv(project / "data/splits/fold_assignments.tsv",
                                                       sep="\t", index=False)
    splits = []
    for test_fold in range(1, 11):
        train, valid, test = split_indices(metadata.fold.values, test_fold)
        split_file = metadata[["gene", "window_group"]].copy()
        roles = np.full(len(metadata), "train", dtype=object)
        roles[valid], roles[test] = "validation", "test"
        split_file["role"] = roles
        split_file.to_csv(project / f"data/splits/round_{test_fold:02d}.tsv", sep="\t", index=False)
        splits.append({"round": test_fold, "test_fold": test_fold,
                       "validation_fold": test_fold % 10 + 1,
                       "train_genes": len(train), "validation_genes": len(valid), "test_genes": len(test)})
    pd.DataFrame(splits).to_csv(project / "data/splits/split_summary.tsv", sep="\t", index=False)
    metadata.to_csv(prepared / "genes.tsv", sep="\t", index=False)
    metadata[["chrom", "window_start0_clipped", "window_end0_clipped", "gene"]].assign(
        score=0, strand=metadata.strand).to_csv(prepared / "gene_windows.bed", sep="\t",
                                               index=False, header=False)
    print("Auditing DNA window overlap across folds", flush=True)
    audit = audit_window_overlap(metadata)
    audit.to_csv(project / "data/splits/window_overlap_audit.tsv", sep="\t", index=False)
    if audit.overlapping_test_genes.sum() != 0:
        raise AssertionError("Overlapping DNA windows occur across folds")
    print("Hashing prepared artifacts", flush=True)
    artifact_names = (
        "feature_manifest.tsv", "genes.tsv", "module_sizes.json", "dna_onehot.npy",
        "features_unscaled.npy", "gene_windows.bed", "target.npy",
    )
    hashes = {name: sha256(prepared / name) for name in artifact_names}
    hashes["../splits/fold_assignments.tsv"] = sha256(project / "data/splits/fold_assignments.tsv")
    for path in sorted((project / "data/splits").glob("*.tsv")):
        hashes[f"../splits/{path.name}"] = sha256(path)
    manifest = {
        "gene_count": len(metadata), "feature_count": raw.shape[1], "module_sizes": modules,
        "sequence_length": cfg["sequence_length"], "dna_encoding": "bool one-hot A,C,G,T; N=0000",
        "window_start_relative": cfg["window_start_relative"], "window_end_relative": cfg["window_end_relative"],
        "split_strategy": cfg["split_strategy"], "split_seed": cfg["split_seed"],
        "independent_window_groups": int(metadata.window_group.nunique()),
        "dna_orientation": "transcriptional; negative strand reverse complemented",
        "anchor": f"representative-transcript TSS at sequence index {-cfg['window_start_relative']}",
        "missing_feature_cells": int(np.isnan(x).sum()), "scaled_before_split": False,
        "source_records": provenance, "prepared_sha256": hashes,
        "overlap_audit": audit.to_dict("records"),
        "evaluation_scope": "overlap-grouped gene CV on fixed precomputed Y_PC1; not homology-disjoint or de-novo PCA",
        "feature_scope": "regulatory assays only; no BPQI/basic/sequence/CpG/codon/k-mer gene features",
    }
    write_json(prepared / "manifest.json", manifest)
    shutil.copy2(sources_path, inputs / "data_sources.json")
    print(json.dumps({k: manifest[k] for k in ("gene_count", "feature_count", "module_sizes",
                                              "missing_feature_cells")}, indent=2), flush=True)
