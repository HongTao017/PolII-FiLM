#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""生成TBP在core promoter、promoter、pause和gene body中的motif特征。

四个区域均按5'到3'转录方向定义：
    core_promoter: TSS [-100, +50)
    promoter:      TSS [-1000, +100)
    pause:         TSS [-50, +250)
    body:          TSS [+250, TES)

每个区域输出三个特征：FIMO显著命中数、每kb命中数和最高FIMO score。
未达到FIMO阈值的基因，其命中数和最高分记为0；长度为0时每kb命中数记NA。
"""

from __future__ import annotations

import argparse
import shlex
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from pandas.errors import EmptyDataError


BASE = Path("/mnt/HDD8TB/hongtao/Data_Seq/DATA/POLII_feature_values")
REFERENCE = Path.home() / "sequence features/groseq_project/reference"

REGIONS = {
    "core_promoter": (-100, 50),
    "promoter": (-1000, 100),
    "pause": (-50, 250),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="生成TBP四区域FIMO motif特征表。"
    )
    parser.add_argument(
        "--matrix",
        type=Path,
        default=BASE / "04_model_matrices/POLII_XY.gene_by_feature.tsv",
        help="含gene/chrom/strand/tss_bed0/tes_bed0的基因表。",
    )
    parser.add_argument(
        "--fasta",
        type=Path,
        default=REFERENCE / "GRCh38.p14.genome.fa",
        help="GRCh38基因组FASTA。",
    )
    parser.add_argument(
        "--fai",
        type=Path,
        default=REFERENCE / "GRCh38.p14.genome.fa.fai",
        help="基因组FASTA索引。",
    )
    parser.add_argument(
        "--motif",
        type=Path,
        default=BASE
        / "02_external_regulatory_features/TF_motif/TBP/MA0108.3_TBP.meme",
        help="TBP的MEME格式motif文件。",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=BASE
        / "02_external_regulatory_features/TF_motif/TBP/four_regions",
        help="BED、FASTA、FIMO结果和最终特征表的输出目录。",
    )
    parser.add_argument(
        "--meme-runner",
        default="",
        help=(
            "fasta-get-markov和fimo之前的命令，例如"
            "'/path/to/micromamba run -n meme_tools'；若工具已在PATH中则留空。"
        ),
    )
    parser.add_argument("--threshold", type=float, default=1e-4)
    parser.add_argument("--max-stored-scores", type=int, default=1_000_000)
    parser.add_argument(
        "--reuse",
        action="store_true",
        help="若对应fimo.tsv已经存在则直接汇总，不重新运行该区域。",
    )
    return parser.parse_args()


def require_file(path: Path, label: str) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"找不到{label}：{path}")


def run(command: list[str]) -> None:
    print("运行：", shlex.join(command), flush=True)
    subprocess.run(command, check=True)


def load_chrom_sizes(fai_path: Path) -> dict[str, int]:
    sizes: dict[str, int] = {}
    with fai_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                sizes[fields[0]] = int(fields[1])
    return sizes


def relative_window(
    tss_bed0: int, strand: str, relative_start: int, relative_end: int
) -> tuple[int, int]:
    """将相对TSS窗口换成BED的0-based半开坐标。"""
    if strand == "-":
        return tss_bed0 - relative_end + 1, tss_bed0 - relative_start + 1
    return tss_bed0 + relative_start, tss_bed0 + relative_end


def build_region_beds(
    matrix: pd.DataFrame, chrom_sizes: dict[str, int], output_dir: Path
) -> tuple[dict[str, Path], dict[str, pd.Series]]:
    """生成四个区域的BED6，并记录每个基因实际用于归一化的区域长度。"""
    required = ["gene", "chrom", "strand", "tss_bed0", "tes_bed0"]
    missing = [column for column in required if column not in matrix.columns]
    if missing:
        raise ValueError(f"输入矩阵缺少列：{', '.join(missing)}")
    if matrix["gene"].duplicated().any():
        examples = matrix.loc[matrix["gene"].duplicated(), "gene"].head().tolist()
        raise ValueError(f"输入矩阵存在重复基因，例如：{examples}")

    bed_paths: dict[str, Path] = {}
    lengths: dict[str, pd.Series] = {}
    region_rows: dict[str, list[tuple[object, ...]]] = {
        name: [] for name in (*REGIONS, "body")
    }
    region_lengths: dict[str, dict[str, float]] = {
        name: {} for name in region_rows
    }

    for row in matrix.loc[:, required].itertuples(index=False):
        gene = str(row.gene)
        chrom = str(row.chrom)
        strand = str(row.strand)
        if chrom not in chrom_sizes or strand not in {"+", "-"}:
            for name in region_rows:
                region_lengths[name][gene] = np.nan
            continue

        tss = int(row.tss_bed0)
        tes = int(row.tes_bed0)
        chrom_length = chrom_sizes[chrom]
        span = abs(tes - tss) + 1

        definitions = dict(REGIONS)
        definitions["body"] = (250, span)
        for name, (relative_start, relative_end) in definitions.items():
            if relative_end <= relative_start:
                region_lengths[name][gene] = 0.0
                continue
            start, end = relative_window(
                tss, strand, relative_start, relative_end
            )
            start = max(0, start)
            end = min(chrom_length, end)
            length = max(0, end - start)
            region_lengths[name][gene] = float(length)
            if length == 0:
                continue
            sequence_name = f"{gene}|{name}({strand})"
            region_rows[name].append(
                (chrom, start, end, sequence_name, 0, strand)
            )

    genes = matrix["gene"].astype(str).reset_index(drop=True)
    for name, rows in region_rows.items():
        path = output_dir / f"gene_{name}.bed"
        with path.open("w", encoding="utf-8", newline="\n") as handle:
            for fields in sorted(rows, key=lambda item: (str(item[0]), int(item[1]))):
                handle.write("\t".join(map(str, fields)) + "\n")
        bed_paths[name] = path
        lengths[name] = genes.map(region_lengths[name]).astype(float)

    return bed_paths, lengths


def resolve_tool(prefix: list[str], tool: str) -> list[str]:
    if prefix:
        return [*prefix, tool]
    executable = shutil.which(tool)
    if executable is None:
        raise FileNotFoundError(
            f"找不到{tool}。请通过--meme-runner指定micromamba运行前缀。"
        )
    return [executable]


def run_region_scan(
    region: str,
    bed_path: Path,
    fasta_path: Path,
    motif_path: Path,
    output_dir: Path,
    meme_prefix: list[str],
    threshold: float,
    max_stored_scores: int,
    reuse: bool,
) -> Path:
    region_fasta = output_dir / f"gene_{region}.fa"
    background = output_dir / f"{region}_background.txt"
    fimo_dir = output_dir / f"fimo_{region}"
    fimo_tsv = fimo_dir / "fimo.tsv"

    if reuse and fimo_tsv.is_file() and fimo_tsv.stat().st_size > 0:
        print(f"复用已有结果：{fimo_tsv}")
        return fimo_tsv

    if fimo_dir.exists():
        raise FileExistsError(
            f"输出目录已经存在：{fimo_dir}。删除旧目录后重跑，或使用--reuse。"
        )

    run(
        [
            "bedtools",
            "getfasta",
            "-fi",
            str(fasta_path),
            "-bed",
            str(bed_path),
            "-fo",
            str(region_fasta),
            "-name",
            "-s",
        ]
    )
    with background.open("w", encoding="utf-8", newline="\n") as handle:
        command = [
            *resolve_tool(meme_prefix, "fasta-get-markov"),
            "-m",
            "1",
            str(region_fasta),
        ]
        print("运行：", shlex.join(command), ">", background, flush=True)
        subprocess.run(command, check=True, stdout=handle)

    run(
        [
            *resolve_tool(meme_prefix, "fimo"),
            "--bgfile",
            str(background),
            "--thresh",
            str(threshold),
            "--max-stored-scores",
            str(max_stored_scores),
            "--oc",
            str(fimo_dir),
            str(motif_path),
            str(region_fasta),
        ]
    )
    require_file(fimo_tsv, f"{region} FIMO结果")
    return fimo_tsv


def summarize_fimo(
    fimo_path: Path,
    genes: pd.Series,
    lengths: pd.Series,
    region: str,
) -> pd.DataFrame:
    """把FIMO逐位点结果汇总成每个基因的三个数值特征。"""
    try:
        frame = pd.read_csv(
            fimo_path, sep="\t", comment="#", low_memory=False
        )
    except EmptyDataError:
        frame = pd.DataFrame(columns=["sequence_name", "score"])

    required = {"sequence_name", "score"}
    if not required.issubset(frame.columns):
        # 某些FIMO版本在零命中时只生成空文件或注释，不写表头。
        if frame.empty:
            frame = pd.DataFrame(columns=["sequence_name", "score"])
        else:
            raise ValueError(f"{fimo_path}缺少FIMO列：{sorted(required)}")

    frame = frame.loc[frame["sequence_name"].notna()].copy()
    frame["gene"] = frame["sequence_name"].astype(str).str.split("|", n=1).str[0]
    frame["score"] = pd.to_numeric(frame["score"], errors="coerce")
    frame = frame.loc[frame["score"].notna()]
    if frame.empty:
        print(
            f"警告：{region}区域没有达到FIMO阈值的TBP motif，"
            "该区域的count和max_score将记为0。"
        )

    grouped = frame.groupby("gene", sort=False)["score"].agg(["size", "max"])
    result = pd.DataFrame({"gene": genes.astype(str).reset_index(drop=True)})
    result = result.join(grouped, on="gene")
    result["size"] = result["size"].fillna(0).astype(int)
    result["max"] = result["max"].fillna(0.0).astype(float)

    length_values = pd.to_numeric(lengths, errors="coerce").reset_index(drop=True)
    per_kb = np.where(
        length_values > 0,
        result["size"].to_numpy(dtype=float) * 1000.0 / length_values,
        np.nan,
    )
    return pd.DataFrame(
        {
            "gene": result["gene"],
            f"TBP_{region}_motif_count": result["size"],
            f"TBP_{region}_motif_count_per_kb": per_kb,
            f"TBP_{region}_motif_max_score": result["max"],
        }
    )


def main() -> None:
    args = parse_args()
    require_file(args.matrix, "基因矩阵")
    require_file(args.fasta, "GRCh38 FASTA")
    require_file(args.fai, "FASTA索引")
    require_file(args.motif, "TBP MEME motif")
    if shutil.which("bedtools") is None:
        raise FileNotFoundError("找不到bedtools")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    matrix = pd.read_csv(args.matrix, sep="\t", low_memory=False)
    matrix["gene"] = matrix["gene"].astype(str)
    chrom_sizes = load_chrom_sizes(args.fai)
    bed_paths, lengths = build_region_beds(matrix, chrom_sizes, args.output_dir)
    meme_prefix = shlex.split(args.meme_runner)

    region_tables: list[pd.DataFrame] = []
    for region in (*REGIONS, "body"):
        fimo_path = run_region_scan(
            region=region,
            bed_path=bed_paths[region],
            fasta_path=args.fasta,
            motif_path=args.motif,
            output_dir=args.output_dir,
            meme_prefix=meme_prefix,
            threshold=args.threshold,
            max_stored_scores=args.max_stored_scores,
            reuse=args.reuse,
        )
        region_tables.append(
            summarize_fimo(
                fimo_path,
                matrix["gene"],
                lengths[region],
                region,
            )
        )

    merged = region_tables[0]
    for table in region_tables[1:]:
        merged = merged.merge(table, on="gene", how="left", validate="one_to_one")

    output_path = args.output_dir / "TBP_four_region_motif_features.tsv"
    merged.to_csv(output_path, sep="\t", index=False, na_rep="NA")
    print(f"完成：{len(merged):,}个基因，12个TBP特征")
    print(f"输出：{output_path}")


if __name__ == "__main__":
    main()
