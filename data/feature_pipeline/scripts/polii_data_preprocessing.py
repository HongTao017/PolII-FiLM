#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""按照用户的X标签顺序构建Pol II暂停指数预测所需的基因特征表。

脚本借鉴两套公开方法，但不会混淆它们的定义：

1. Akcan等人的方法：从转录本基因组区间构建长度、外显子数量、平均
   外显子长度、外显子比例和完整转录本GC等特征，并将部分连续变量
   缩放到0到1。
2. Agarwal和Kelley的方法：计算5'UTR、CDS/ORF、3'UTR的长度和GC、
   内含子总长度、ORF外显子连接密度、k-mer频率和密码子频率。

输入文件：
    MANE GTF、GRCh38 FASTA及其FAI、可选的基因PC1表和PI矩阵。

主要输出：
    <prefix>.gene_by_feature.tsv
        每行一个基因，每列一个特征，适合机器学习程序读取。
    <prefix>.feature_by_gene.tsv
        每行一个特征，每列一个基因，顺序对应用户现有数据表。
    <prefix>.akcan_model_features.tsv
        将当前可以复现的特征转换成Akcan模型使用的表达方式。
    <prefix>.annotation_report.tsv
        记录每个目标基因是否成功匹配MANE转录本。

没有CAGE、CpG岛、RNA-seq、ATAC、ChIP或eCLIP文件时，脚本不会把这些
特征伪造为0；它们不会进入当前可计算X表。

/home/user_011/miniconda3/envs/RAN1/bin/python \
  "/mnt/HDD8TB/hongtao/Data_Seq/DATA/POLII_feature_values/scripts/polii_data_preprocessing.py" \
  --gtf "/home/user_011/sequence features/groseq_project/reference/MANE.GRCh38.v1.5.standard_chromosomes.gtf.gz" \
  --fasta "/home/user_011/sequence features/groseq_project/reference/GRCh38.p14.genome.fa" \
  --fai "/home/user_011/sequence features/groseq_project/reference/GRCh38.p14.genome.fa.fai" \
  --pc1 "/mnt/HDD8TB/hongtao/Data_Seq/DATA/POLII_feature_values/03_prediction_targets/PI_XY_feature_template_genes_as_rows.tsv" \
  --output-prefix "/mnt/HDD8TB/hongtao/Data_Seq/DATA/POLII_feature_values/04_model_matrices/POLII_XY"
"""

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


# 这里严格保留用户数据表中当前能够计算的基础字段顺序。
# k-mer、密码子和单样本PI会在各自占位位置展开成多列。
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
    """保存一个转录本从GTF中解析出的全部注释区域。

    坐标统一保存为BED风格的0-based半开区间[start, end)，这样区间长度
    永远是end-start，能够避免GTF的1-based闭区间产生加1或减1错误。
    """

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
        """将Ensembl和RefSeq对编码转录本的常见标签统一为1。

        Ensembl通常写作protein_coding，MANE中的RefSeq注释可能写作mRNA。
        mRNA在这里表示具有CDS的蛋白编码信使RNA，不能判成非编码转录本。
        """
        values = {
            self.gene_biotype.strip().lower(),
            self.transcript_biotype.strip().lower(),
        }
        coding_labels = {"protein_coding", "protein-coding", "mrna"}
        return int(bool(values & coding_labels))


class IndexedFasta:
    """使用FAI索引随机读取大型FASTA，不把3 GB基因组全部载入内存。

    输入是FASTA和samtools生成的.fai文件。fetch()接受0-based半开坐标，
    返回对应的基因组正链序列。链方向的反向互补在后续函数中处理。
    """

    def __init__(self, fasta_path: Path, fai_path: Optional[Path] = None) -> None:
        self.fasta_path = fasta_path
        self.fai_path = fai_path or Path(f"{fasta_path}.fai")
        if not self.fasta_path.exists():
            raise FileNotFoundError(f"找不到FASTA文件：{self.fasta_path}")
        if not self.fai_path.exists():
            raise FileNotFoundError(f"找不到FAI文件：{self.fai_path}")

        self.index: Dict[str, Tuple[int, int, int, int]] = {}
        with self.fai_path.open("r", encoding="utf-8") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 5:
                    continue
                name = fields[0]
                self.index[name] = tuple(map(int, fields[1:5]))  # type: ignore[assignment]
        self.handle = self.fasta_path.open("rb")

    def close(self) -> None:
        self.handle.close()

    def __enter__(self) -> "IndexedFasta":
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        self.close()

    def resolve_contig(self, chrom: str) -> str:
        """处理chr1与1等常见染色体命名差异，返回FASTA中的真实名称。"""
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
        raise KeyError(f"FASTA中找不到与{chrom!r}对应的染色体")

    def fetch(self, chrom: str, start0: int, end0: int) -> str:
        """读取指定区间；越过染色体边界的部分会自动裁剪。"""
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
    """根据扩展名透明读取普通GTF或gzip压缩GTF。"""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def strip_version(identifier: str) -> str:
    """删除Ensembl或RefSeq编号末尾的版本号，便于不同表格之间匹配。"""
    return re.sub(r"\.\d+$", "", str(identifier or ""))


def parse_gtf_attributes(text: str) -> Dict[str, List[str]]:
    """解析GTF第9列，保留tag等可能重复出现的属性。"""
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
    """按候选字段顺序返回第一个存在的GTF属性值。"""
    for key in keys:
        values = attributes.get(key)
        if values:
            return str(values[0])
    return ""


def merge_intervals(intervals: Iterable[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """合并重叠区间，避免同一转录本的重复注释造成长度重复计算。"""
    ordered = sorted((int(start), int(end)) for start, end in intervals if end > start)
    merged: List[List[int]] = []
    for start, end in ordered:
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


def parse_mane_gtf(gtf_path: Path) -> Dict[str, TranscriptRecord]:
    """读取MANE GTF并构建以transcript_id为键的转录本对象。

    处理结果包括转录本范围、基因名称、链方向、外显子和CDS区间。
    与Akcan的data_preprocessing.R相同，这一步只生成标准化区域；长度、
    比例和序列特征由后续函数计算。
    """
    if not gtf_path.exists():
        raise FileNotFoundError(f"找不到GTF文件：{gtf_path}")

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
    """给同一基因的多个候选转录本排序。

    首先优先MANE Select；若仍有多个候选，依次选择CDS最长、拼接长度
    最长和基因组跨度最长的转录本。这与Agarwal优先最长ORF的思想一致，
    同时保证MANE Select优先于其他转录本。
    """
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
    """按目标PC1表的基因顺序选择MANE转录本。

    如果目标表提供转录本编号，优先匹配该编号或无版本号编号；匹配不到
    时使用MANE Select和最长CDS规则。返回基因顺序、所选转录本以及匹配
    报告，确保最终X表和Y表顺序完全一致。
    """
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
    """计算一组0-based半开区间合并后的总长度。"""
    return sum(end - start for start, end in merge_intervals(intervals))


def intervals_in_transcription_order(
    intervals: Sequence[Tuple[int, int]], strand: str
) -> List[Tuple[int, int]]:
    """按照真实转录方向排列区间，负链从高坐标向低坐标排列。"""
    return sorted(intervals, key=lambda item: item[0], reverse=(strand == "-"))


def reverse_complement(sequence: str) -> str:
    """将基因组正链序列转换为负链转录方向序列。"""
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return sequence.translate(table)[::-1].upper()


def fetch_intervals_sequence(
    fasta: IndexedFasta,
    chrom: str,
    intervals: Sequence[Tuple[int, int]],
    strand: str,
) -> str:
    """提取并拼接多个区间，结果始终按5'到3'转录方向排列。"""
    merged = merge_intervals(intervals)
    sequence = "".join(fasta.fetch(chrom, start, end) for start, end in merged)
    return reverse_complement(sequence) if strand == "-" else sequence


def intersect_interval(
    interval: Tuple[int, int], boundary: Tuple[int, int]
) -> Optional[Tuple[int, int]]:
    """返回两个半开区间的交集，无交集时返回None。"""
    start = max(interval[0], boundary[0])
    end = min(interval[1], boundary[1])
    return (start, end) if end > start else None


def infer_utr_intervals(record: TranscriptRecord) -> Tuple[List[Tuple[int, int]], List[Tuple[int, int]]]:
    """通过外显子和CDS推导5'UTR及3'UTR区间。

    用户的MANE GTF没有独立UTR行，因此先确定CDS最小和最大基因组坐标，
    再截取CDS两侧的外显子部分。负链中高坐标侧是5'UTR，低坐标侧是
    3'UTR。无CDS时返回两个空列表，不把整条转录本误判为UTR。
    """
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
    """从相邻外显子的间隙构建内含子区间。"""
    exons = sorted(merge_intervals(record.exons))
    return [
        (left[1], right[0])
        for left, right in zip(exons, exons[1:])
        if right[0] > left[1]
    ]


def safe_fraction(numerator: float, denominator: float) -> float:
    """安全计算比例；分母为0时返回缺失值而不是伪造0。"""
    return float(numerator / denominator) if denominator else float("nan")


def summarize_lengths(lengths: Sequence[int]) -> Tuple[float, float, float]:
    """返回平均、最小和最大长度；没有区间时返回三个NA。"""
    if not lengths:
        return float("nan"), float("nan"), float("nan")
    return float(np.mean(lengths)), float(min(lengths)), float(max(lengths))


def build_annotation_and_structure_features(
    gene: str,
    record: TranscriptRecord,
    output_transcript: Optional[str] = None,
) -> Tuple[Dict[str, object], Dict[str, List[Tuple[int, int]]]]:
    """依照用户表格第1到36项计算注释和结构特征。

    返回值一是能够直接写入X表的数值字典；返回值二是后续序列计算所需
    的外显子、内含子、CDS和UTR区间。单外显子基因的内含子统计值设为
    NA，只有intron_count明确为0。
    """
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
    """把相对TSS的转录方向窗口转换为基因组0-based半开区间。

    例如[-50, 250)表示TSS上游50 bp到下游250 bp。负链会自动翻转，
    从而保证正负链提取的序列都按转录方向解释。
    """
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
    """提取相对TSS窗口并转换到5'到3'转录方向。"""
    start0, end0 = relative_tss_window(tss_bed0, strand, relative_start, relative_end)
    sequence = fasta.fetch(chrom, start0, end0)
    return reverse_complement(sequence) if strand == "-" else sequence


def sequence_composition(sequence: str) -> Tuple[float, float]:
    """计算GC和AT比例，只使用A/C/G/T有效碱基作为分母。"""
    counts = Counter(base for base in sequence.upper() if base in DNA_ALPHABET)
    total = sum(counts.values())
    if total == 0:
        return float("nan"), float("nan")
    return (counts["G"] + counts["C"]) / total, (counts["A"] + counts["T"]) / total


def cpg_features(sequence: str) -> Tuple[float, float]:
    """计算CpG数量和观察/期望值。

    O/E公式为CpG_count×有效序列长度/(C_count×G_count)。当C或G为0时，
    O/E返回NA，避免除零。
    """
    sequence = sequence.upper()
    valid = "".join(base for base in sequence if base in DNA_ALPHABET)
    cpg_count = valid.count("CG")
    c_count = valid.count("C")
    g_count = valid.count("G")
    if not valid or c_count == 0 or g_count == 0:
        return float(cpg_count), float("nan")
    return float(cpg_count), cpg_count * len(valid) / (c_count * g_count)


def all_kmers(min_k: int = 1, max_k: int = 4) -> List[str]:
    """按照长度优先、字母顺序生成全部k-mer名称。"""
    return [
        "".join(parts)
        for k in range(min_k, max_k + 1)
        for parts in itertools.product(DNA_ALPHABET, repeat=k)
    ]


def kmer_frequencies(sequence: str, min_k: int = 1, max_k: int = 4) -> Dict[str, float]:
    """采用Agarwal的step=1重叠滑窗计算每种k-mer频率。

    含N的窗口不进入分母。每个k长度分别归一化，所以同一k下的全部频率
    之和为1；有效窗口为0时该长度的所有频率均为NA。
    """
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
    """采用Agarwal的width=3、step=3方式计算61个非终止密码子频率。

    CDS必须先按转录方向拼接。终止密码子和含N密码子不进入分母，结果中
    61个密码子频率之和在存在有效密码子时为1。
    """
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
    """计算区域序列、CpG、pause/body k-mer和密码子特征。

    promoter、pause和early_body窗口由命令行参数统一控制。所有序列均按
    转录方向处理，负链自动反向互补。函数返回可直接合并到该基因X行的
    特征字典。1-4-mer只在pause与完整gene body中计算；gene body从
    pause窗口终点开始，一直延伸到MANE转录本TES。
    """
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
    """不区分大小写地寻找常见列名，适配不同版本的PC1和PI表。"""
    lookup = {str(column).lower(): str(column) for column in frame.columns}
    for candidate in candidates:
        if candidate.lower() in lookup:
            return lookup[candidate.lower()]
    if required:
        raise ValueError(f"没有找到必需列，候选列为：{', '.join(candidates)}")
    return None


def load_target_table(path: Optional[Path]) -> Tuple[Optional[pd.DataFrame], Optional[str], Optional[str], Optional[str]]:
    """读取基因PC1表并识别gene、transcript和Y列。

    返回的表会去除重复基因，只保留首次出现，从而保证每个基因对应一个
    最终Y_PC1值。
    """
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
    """解析重复的LABEL=PATH参数，例如GRO_seq=/data/gro.tsv。"""
    output: List[Tuple[str, Path]] = []
    for item in items:
        if "=" not in item:
            raise ValueError(f"PI矩阵参数必须是LABEL=PATH格式：{item}")
        label, raw_path = item.split("=", 1)
        output.append((label.strip(), Path(raw_path).expanduser()))
    return output


def load_pi_features(
    labelled_paths: Sequence[Tuple[str, Path]], gene_order: Sequence[str]
) -> Tuple[pd.DataFrame, List[str]]:
    """读取一个或多个基因×样本PI矩阵并按目标基因顺序对齐。

    每个样本列重命名为PI_<测序标签>_<原样本名>，防止不同矩阵中出现
    同名样本。PI只用于保留现有数据或复核，不应作为预测Y_PC1的X输入。
    """
    merged = pd.DataFrame(index=pd.Index(gene_order, name="gene"))
    output_columns: List[str] = []
    for label, path in labelled_paths:
        if not path.exists():
            raise FileNotFoundError(f"找不到PI矩阵：{path}")
        frame = pd.read_csv(path, sep="\t", low_memory=False)
        gene_column = detect_column(
            frame, ["gene", "gene_name", "symbol"], required=False
        )
        combined_column = detect_column(
            frame, ["gene_transcript"], required=False
        )
        if gene_column is None and combined_column is None:
            raise ValueError(
                f"PI矩阵{path}中既没有gene列，也没有gene_transcript列"
            )

        # GRO-seq和PRO-seq矩阵有独立gene列；当前mNET-seq矩阵使用
        # gene_transcript，例如SAMD11|NM_001385641.1，需要先拆出基因名。
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

        # 前6个注释字段不是PI样本，必须从数值矩阵中排除。转换后整列均为
        # NA的字段也删除，防止其他文本元数据误进入模型表。
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
    """复现Akcan使用的rescale逻辑，将非缺失值线性缩放到0到1。"""
    numeric = pd.to_numeric(series, errors="coerce")
    minimum = numeric.min(skipna=True)
    maximum = numeric.max(skipna=True)
    if pd.isna(minimum) or pd.isna(maximum) or maximum == minimum:
        return pd.Series(np.nan, index=series.index, dtype=float)
    return (numeric - minimum) / (maximum - minimum)


def chromosome_numeric(series: pd.Series) -> pd.Series:
    """按Akcan思路将染色体转为数字；非标准染色体保留为NA。

    该列仅用于复现论文。正式模型更推荐对染色体使用one-hot编码。
    """
    mapping = {str(index): index for index in range(1, 23)}
    mapping.update({"X": 23, "Y": 24, "M": 25, "MT": 25})
    cleaned = series.astype(str).str.replace(r"^chr", "", regex=True)
    return cleaned.map(mapping).astype(float)


def build_akcan_model_view(frame: pd.DataFrame) -> pd.DataFrame:
    """从当前可计算列生成Akcan风格模型特征表。

    tx.len、tx.ex.count和tx.ex.width按作者代码缩放到0到1。tx.gc.seq使用
    包含内含子的完整基因组转录本序列GC。当前脚本另外保留gene作为连接
    键，但它不应作为数值特征输入模型。
    """
    output = pd.DataFrame(index=frame.index)
    output["gene"] = frame["gene"]

    # R中as.numeric(strand(GRanges))返回strand因子的内部编码，通常是
    # '+'=1、'-'=2、'*'=3；这里按作者代码复现，而不是改成二元0/1。
    output["tx.strand"] = frame["strand"].map({"+": 1, "-": 2, "*": 3})
    output["tx.gc.seq"] = frame["transcript_genomic_gc_fraction"]
    scaled_length = minmax_scale(frame["transcript_span_bp"])
    scaled_exon_count = minmax_scale(frame["exon_count"])
    output["tx.len"] = scaled_length
    output["tx.ex.count"] = scaled_exon_count

    # 作者代码先分别缩放转录本长度和外显子数，再计算二者之比并再次缩放。
    # 最小外显子数缩放后可能为0，因此这里将相应除零结果保存为NA。
    scaled_ratio = scaled_length / scaled_exon_count.replace(0, np.nan)
    output["tx.ex.ratio"] = minmax_scale(scaled_ratio)
    output["tx.ex.width"] = minmax_scale(frame["mean_exon_length_bp"])
    output["tx.ex.seq"] = frame["exon_fraction"]

    # 作者使用start(GRanges)，即不考虑链方向的基因组低坐标起点。
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
    """生成最终X表顺序，并在原占位位置展开k-mer、密码子和PI样本列。"""
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
    """写出基因×特征、特征×基因、Akcan复现表和注释匹配报告。"""
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    gene_by_feature_path = Path(f"{output_prefix}.gene_by_feature.tsv")
    feature_by_gene_path = Path(f"{output_prefix}.feature_by_gene.tsv")
    akcan_path = Path(f"{output_prefix}.akcan_model_features.tsv")
    report_path = Path(f"{output_prefix}.annotation_report.tsv")

    gene_by_feature = frame[["gene", *feature_columns]].copy()
    gene_by_feature.to_csv(gene_by_feature_path, sep="\t", index=False, na_rep="NA")

    # 转置后第一列是feature_or_columns，其余列按目标PC1表中的基因顺序排列。
    transposed = gene_by_feature.set_index("gene").T
    transposed.index.name = "feature_or_columns"
    transposed.to_csv(feature_by_gene_path, sep="\t", index=True, na_rep="NA")

    build_akcan_model_view(frame).to_csv(akcan_path, sep="\t", index=False, na_rep="NA")
    report.to_csv(report_path, sep="\t", index=False, na_rep="NA")

    print(f"已写出：{gene_by_feature_path}")
    print(f"已写出：{feature_by_gene_path}")
    print(f"已写出：{akcan_path}")
    print(f"已写出：{report_path}")


def build_argument_parser() -> argparse.ArgumentParser:
    """定义命令行参数，并把所有窗口边界明确写入帮助信息。"""
    parser = argparse.ArgumentParser(
        description="按照指定X标签顺序，从MANE GTF和GRCh38 FASTA构建基因特征表。"
    )
    parser.add_argument("--gtf", type=Path, required=True, help="MANE GTF或GTF.gz路径")
    parser.add_argument("--fasta", type=Path, required=True, help="GRCh38 FASTA路径")
    parser.add_argument("--fai", type=Path, default=None, help="FASTA的.fai路径；默认使用FASTA+.fai")
    parser.add_argument(
        "--pc1",
        type=Path,
        default=None,
        help="基因PC1表；用于确定基因顺序并将PC1_Y写为最后一列Y_PC1",
    )
    parser.add_argument(
        "--pi-matrix",
        action="append",
        default=[],
        metavar="LABEL=PATH",
        help="可重复提供PI矩阵，例如GRO_seq=/path/gro.tsv",
    )
    parser.add_argument(
        "--output-prefix",
        type=Path,
        required=True,
        help="输出文件前缀，例如results/POLII_XY",
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
        help="不展开pause/body共680个1-4-mer和61个密码子列，用于快速检查注释流程",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    """串联全部处理步骤并生成最终结果。

    处理顺序严格对应用户表格：身份信息、坐标、转录本结构、外显子、
    内含子、CDS/UTR、比例、剪接、序列组成、CpG、k-mer、密码子、PI，
    最后写入Y_PC1。
    """
    args = build_argument_parser().parse_args(argv)
    target_table, gene_column, transcript_column, y_column = load_target_table(args.pc1)

    print("正在读取MANE GTF并构建转录本区域……")
    records = parse_mane_gtf(args.gtf)
    if not records:
        raise RuntimeError("GTF中没有解析出任何转录本，请检查文件格式和transcript_id属性")

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
    print(f"正在计算{len(gene_order):,}个目标基因的结构和序列特征……")
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

            # Akcan的tx.gc.seq是完整基因组转录本范围GC，明确包含内含子。
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
                print(f"  已处理{index:,}/{len(gene_order):,}个基因")

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
    print(f"注释匹配：{matched:,}/{len(report):,}个基因")
    if matched < len(report):
        print("未匹配基因已保留在输出中，其不可计算特征为NA；详情见annotation_report。")
    print("提醒：PI样本列与Y_PC1直接相关，训练预测Y_PC1时必须从X中排除。")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"错误：{error}", file=sys.stderr)
        raise
