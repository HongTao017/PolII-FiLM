#!/usr/bin/env python3
"""Regroup cleaned POLII features into B/Q/I/C/P/E/T matrices.

Existing values are copied without recalculation, normalization, or imputation.
G5 is split into pause and gene-body k-mer groups using its column prefixes.
"""

from __future__ import annotations

import argparse
from collections import OrderedDict
from pathlib import Path

import pandas as pd


SOURCE_DIR = Path(
    "/mnt/HDD8TB/hongtao/Data_Seq/DATA/POLII_feature_values/"
    "04_model_matrices/feature_groups_cleaned"
)
OUTPUT_DIR = Path(
    "/mnt/HDD8TB/hongtao/Data_Seq/DATA/POLII_feature_values/"
    "04_model_matrices/feature_groups_BQICPET"
)
ID_COLUMN = "gene"
TARGET_COLUMN = "Y_PC1"

SOURCES = OrderedDict(
    [
        ("B", "G1_basic_gene_features_cleaned_22X.tsv"),
        ("Q", "G2_sequence_composition_cleaned_13X.tsv"),
        ("I", "G3_CpG_island_cleaned_8X.tsv"),
        ("C", "G4_codon_frequency_cleaned_61X.tsv"),
        ("K", "G5_pause_body_kmer_cleaned_680X.tsv"),
        ("T", "G6_CTCF_intron_3prime.tsv"),
    ]
)
EXPECTED_COUNTS = OrderedDict(
    [("B", 22), ("Q", 13), ("I", 8), ("C", 61), ("P", 340), ("E", 340), ("T", 2)]
)
GROUP_NAMES = OrderedDict(
    [
        ("B", ("basic_gene_structure", "基础基因结构")),
        ("Q", ("sequence_composition", "序列组成")),
        ("I", ("CpG_island", "CpG岛")),
        ("C", ("codon_frequency", "密码子频率")),
        ("P", ("pause_region_1to4mer", "pause区域1-4mer")),
        ("E", ("gene_body_region_1to4mer", "gene-body区域1-4mer")),
        ("T", ("CTCF_predicted_binding", "CTCF预测结合")),
    ]
)
OUTPUT_FILES = OrderedDict(
    [
        ("B", "B_basic_gene_structure_22X.tsv"),
        ("Q", "Q_sequence_composition_13X.tsv"),
        ("I", "I_CpG_island_8X.tsv"),
        ("C", "C_codon_frequency_61X.tsv"),
        ("P", "P_pause_region_1to4mer_340X.tsv"),
        ("E", "E_gene_body_region_1to4mer_340X.tsv"),
        ("T", "T_CTCF_predicted_binding_2X.tsv"),
    ]
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Regroup cleaned POLII features as B/Q/I/C/P/E/T."
    )
    parser.add_argument("--source-dir", type=Path, default=SOURCE_DIR)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    return parser.parse_args()


def read_table(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Missing source table: {path}")
    frame = pd.read_csv(path, sep="\t", low_memory=False)
    if ID_COLUMN not in frame.columns:
        raise ValueError(f"{path.name} does not contain '{ID_COLUMN}'")
    if frame[ID_COLUMN].isna().any() or frame[ID_COLUMN].duplicated().any():
        raise ValueError(f"{path.name} has empty or duplicated gene identifiers")
    return frame


def feature_columns(frame: pd.DataFrame) -> list[str]:
    return [
        column
        for column in frame.columns
        if column not in {ID_COLUMN, TARGET_COLUMN}
    ]


def validate_alignment(
    reference: pd.DataFrame, current: pd.DataFrame, source_name: str
) -> None:
    if not reference[ID_COLUMN].equals(current[ID_COLUMN]):
        raise ValueError(f"Gene order differs in {source_name}")
    if TARGET_COLUMN not in current.columns:
        return
    left = pd.to_numeric(reference[TARGET_COLUMN], errors="coerce")
    right = pd.to_numeric(current[TARGET_COLUMN], errors="coerce")
    if not (left.eq(right) | (left.isna() & right.isna())).all():
        raise ValueError(f"Y_PC1 differs in {source_name}")


def main() -> None:
    args = parse_args()
    tables = OrderedDict(
        (group, read_table(args.source_dir / filename))
        for group, filename in SOURCES.items()
    )
    reference = tables["B"]
    if TARGET_COLUMN not in reference.columns:
        raise ValueError(f"Reference table must contain {TARGET_COLUMN}")
    for group, frame in tables.items():
        validate_alignment(reference, frame, SOURCES[group])

    kmer_columns = feature_columns(tables["K"])
    pause_columns = [c for c in kmer_columns if c.lower().startswith("pause_kmer_")]
    body_columns = [c for c in kmer_columns if c.lower().startswith("body_kmer_")]
    classified = set(pause_columns + body_columns)
    unclassified = [c for c in kmer_columns if c not in classified]
    if unclassified:
        raise ValueError(f"Unclassified G5 columns: {unclassified[:10]}")

    groups: OrderedDict[str, list[str]] = OrderedDict(
        [
            ("B", feature_columns(tables["B"])),
            ("Q", feature_columns(tables["Q"])),
            ("I", feature_columns(tables["I"])),
            ("C", feature_columns(tables["C"])),
            ("P", pause_columns),
            ("E", body_columns),
            ("T", feature_columns(tables["T"])),
        ]
    )
    observed = {group: len(columns) for group, columns in groups.items()}
    if observed != dict(EXPECTED_COUNTS):
        raise ValueError(
            f"Unexpected group sizes. Expected {dict(EXPECTED_COUNTS)}, "
            f"observed {observed}"
        )

    all_features = [column for columns in groups.values() for column in columns]
    index = pd.Index(all_features)
    duplicate_features = index[index.duplicated()].unique().tolist()
    if duplicate_features:
        raise ValueError(f"Features assigned to multiple groups: {duplicate_features[:10]}")
    if len(all_features) != 786:
        raise ValueError(f"Expected 786 X features, observed {len(all_features)}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    output_groups: OrderedDict[str, pd.DataFrame] = OrderedDict()
    for group, columns in groups.items():
        source_group = "K" if group in {"P", "E"} else group
        frame = pd.concat(
            [
                tables[source_group][[ID_COLUMN]],
                tables[source_group][columns],
                reference[[TARGET_COLUMN]],
            ],
            axis=1,
        )
        frame.to_csv(
            args.output_dir / OUTPUT_FILES[group],
            sep="\t",
            index=False,
            na_rep="NA",
        )
        output_groups[group] = frame

    combined = pd.concat(
        [reference[[ID_COLUMN]]]
        + [output_groups[group][groups[group]] for group in groups]
        + [reference[[TARGET_COLUMN]]],
        axis=1,
    )
    combined.to_csv(
        args.output_dir / "POLII_BQICPET_786X.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )
    combined.head(50).to_csv(
        args.output_dir / "POLII_BQICPET_786X.preview_50_genes.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )

    group_rows = []
    feature_rows = []
    validation_rows = []
    for group, columns in groups.items():
        english_name, chinese_name = GROUP_NAMES[group]
        group_rows.append(
            {
                "letter": group,
                "group_name": english_name,
                "中文含义": chinese_name,
                "feature_count": len(columns),
                "output_file": OUTPUT_FILES[group],
            }
        )
        for order, column in enumerate(columns, start=1):
            feature_rows.append(
                {
                    "letter": group,
                    "group_name": english_name,
                    "中文含义": chinese_name,
                    "feature_order_within_group": order,
                    "feature_name": column,
                }
            )
        validation_rows.append(
            {
                "letter": group,
                "gene_count": len(output_groups[group]),
                "feature_count": len(columns),
                "missing_X_values": int(
                    output_groups[group][columns].isna().sum().sum()
                ),
                "duplicated_gene_count": int(
                    output_groups[group][ID_COLUMN].duplicated().sum()
                ),
            }
        )

    pd.DataFrame(group_rows).to_csv(
        args.output_dir / "feature_group_dictionary.tsv", sep="\t", index=False
    )
    pd.DataFrame(feature_rows).to_csv(
        args.output_dir / "feature_dictionary_long.tsv", sep="\t", index=False
    )
    validation_rows.append(
        {
            "letter": "ALL",
            "gene_count": len(combined),
            "feature_count": len(all_features),
            "missing_X_values": int(combined[all_features].isna().sum().sum()),
            "duplicated_gene_count": int(combined[ID_COLUMN].duplicated().sum()),
        }
    )
    pd.DataFrame(validation_rows).to_csv(
        args.output_dir / "validation_summary.tsv", sep="\t", index=False
    )

    readme = """# POLII B/Q/I/C/P/E/T 特征分组

该目录只重新划分已有清洗后特征，不重新计算、不归一化、不填补任何数值。

| 字母 | 含义 | X特征数 |
|---|---|---:|
| B | 基础基因结构 | 22 |
| Q | 序列组成：GC及局部CpG | 13 |
| I | CpG岛注释 | 8 |
| C | CDS密码子频率 | 61 |
| P | pause区域1-4mer频率 | 340 |
| E | gene-body区域1-4mer频率 | 340 |
| T | CTCF预测结合 | 2 |

`gene`是行标识，`Y_PC1`是预测目标，两者均不计入X特征。
`transcript/chrom/strand/tss_bed0/tes_bed0/gene_biotype`属于注释信息，
不会重新加入清洗后的数值模型矩阵。

`feature_dictionary_long.tsv`逐列列出了每一个特征实际属于哪个字母。
`validation_summary.tsv`记录每组的基因数、特征数和缺失值数量。
"""
    (args.output_dir / "README.md").write_text(readme, encoding="utf-8")

    print(f"Output directory: {args.output_dir}")
    print(f"Genes: {len(combined):,}")
    print(f"X features: {len(all_features)}")
    print("Groups: " + ", ".join(f"{g}={len(c)}" for g, c in groups.items()))


if __name__ == "__main__":
    main()
