"""Remove explicitly superseded PDFs after verified, recoverable ZIP backup."""
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
from collections import defaultdict
import hashlib
import json
import zipfile

ROOT = Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019').resolve()
OLD_PDFS = [
    'motifStack_11_JASPAR_continuous_ring_5x4.pdf',
    'motifStack_11_JASPAR_continuous_ring_darker_grey_5x4.pdf',
    'motifStack_11_JASPAR_even_labels_5x4.pdf',
    'motifStack_11_JASPAR_hit_status_purple_ring_5x4.pdf',
    'motifStack_11_JASPAR_less_dense_5x4.pdf',
    'motifStack_11_JASPAR_straight_spaced_5x4.pdf',
    'motifStack_JASPAR_reference_matrix_panels.pdf',
    'motifStack_all_314_PFM_matrix_atlas.pdf',
    'motifStack_radial_314_default.pdf',
    'motifStack_reference_style_circular_5x4.pdf',
]


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def checked_file(path, parent):
    assert not path.is_symlink() and path.resolve().parent == parent.resolve()
    assert path.is_file()
    return path.resolve()


def write_text(path, content):
    assert ROOT in path.resolve().parents and not path.is_symlink()
    temp = path.with_name(path.name + '.writing')
    assert not temp.exists()
    temp.write_text(content, encoding='utf-8')
    temp.replace(path)


def write_guide(report, receipt):
    rows = report['naming']
    groups = defaultdict(list)
    for row in rows:
        groups[row['feature_label']].append(row)
    current = [Path(report[k]['output']) for k in ('full_tree', 'matched_tree')]
    lines = [
        '# 当前图件与序列命名说明', '',
        '当前结果为44个已命中pattern的序列特征命名。44个原pattern均有唯一展示名，共29种不同六碱基短特征；这些短特征不等于29个经验证的生物学motif家族。原有79条合格TOMTOM匹配、11种JASPAR参考矩阵继续保留。', '',
        '## 已清理的旧版', '',
        f'已从figures目录删除{receipt["deleted_old_pdf_count"]}个历史PDF，仅保留下列2个最新版。旧版PDF和旧说明已先完整压缩备份并逐文件验证SHA256，可从以下文件恢复：', '',
        '`' + receipt['recoverable_archive'] + '`', '',
        '原始H5、原motifStack树、距离矩阵、命名表、复现脚本及作为脚本输入的历史QA记录保留。删除范围限于列出的旧版PDF，没有删除原始实验或分析数据。', '',
        '## 现在看哪一张图', '',
        '|文件|内容|', '|---|---|',
        f'|`{current[0]}`|314个pattern的总览。全部314个仍在树上，仅标出11个原代表pattern的序列短名；不是总共只有11个pattern或11个序列类群。|',
        f'|`{current[1]}`|44个已命中pattern的详细命名图，共2页，每页22个名称。两页均显示同一棵44成员树，名称交替标注，44个名称各出现一次，不是两组聚类。|', '',
        '详图从原314成员树中去掉270个未命中叶节点并收缩只有一个保留子节点的路径得到，保留44个命中pattern之间的原分叉关系，没有重新聚类。详图中的44个成员均已命中，所以详图外环全为深紫色。', '',
        '## 图中颜色与符号', '',
        '|元素|含义|', '|---|---|',
        '|深紫色外环（#918BB4）|有原有合格JASPAR/TOMTOM匹配的pattern，总览中44个。|',
        '|浅紫色外环（#DDDDEA）|没有达到原匹配阈值的pattern，总览中270个；不表示没有调控功能。|',
        '|内部中灰色扇区/背景（#B0B0B0）|标出已命中的序列pattern范围，与深紫色命中状态对应。|',
        '|粉色分支|正Y_PC1归因方向的pattern或纯正归因子树。|',
        '|蓝色分支|负Y_PC1归因方向的pattern或纯负归因子树。|',
        '|中性灰色内部枝线|该内部子树同时包含正、负归因pattern；灰色枝线与灰色背景扇区的含义不同。|',
        '|黑色三角和短直线|把一个序列名称指向对应的原pattern枝端，不表示转录方向或TF激活/抑制。|',
        '|树的分叉|来自原motifStack序列相似性层次聚类。绘制的角度、枝长和扇区宽度用于排版，不是可读取的数值距离或富集量。|', '',
        '正/负归因表示对模型输出Y_PC1的贡献方向，不能直接读作实验测得的表达上调/下调。', '',
        '## 名称怎么读', '',
        '名称来自原H5完整精度的50bp sequence概率矩阵：计算每个位置相对均匀背景的信息量，取信息量总和最高的连续6bp，再生成IUPAC简并共识。它是核心特征短标签，不是完整50bp序列，不保证就是TOMTOM局部对齐所用的片段。', '',
        '连续不少于3个相同字母简写为“次数+字母”，次数只作用于紧随其后的一个字母。例如6T＝TTTTTT、5TA＝TTTTTA、CN3GN＝CNGGGN。后缀-01、-02区分共享同一短特征的不同原pattern，不表示重复次数、显著性或不同TF家族。', '',
        '|简并字母|允许的碱基|', '|---|---|',
        '|N|A/C/G/T；这里共识规则未能用更具体的代码概括，不是原序列缺失值。|',
        '|R|A或G|', '|Y|C或T|', '|M|A或C|', '|S|C或G|', '|W|A或T|', '|K|G或T|', '',
        '允许多个碱基不意味着它们概率完全相等；具体概率仍在原矩阵中。字母定义见[Biopython IUPAC映射](https://github.com/biopython/biopython/blob/master/Bio/Data/IUPACData.py)，共识判定采用[Biopython 1.84的Cavener规则](https://github.com/biopython/biopython/blob/biopython-184/Bio/motifs/matrix.py)。', '',
        '## 全部29种短特征的展开含义（覆盖44个pattern）', '',
        '|图上名称|展开后的6bp特征|pattern数量|该特征下保留的已有TF匹配|', '|---|---|---:|---|',
    ]
    for stem, members in groups.items():
        assert len({x['core_IUPAC_6bp'] for x in members}) == 1
        display = members[0]['sequence_name'] if len(members) == 1 else members[0]['sequence_name'] + '～' + members[-1]['sequence_name']
        tfs = sorted({tf for x in members for tf in x['all_TFs'].split(';')})
        lines.append('|'+ '|'.join([display, members[0]['core_IUPAC_6bp'], str(len(members)), '、'.join(tfs)])+'|')
    lines += [
        '', '以上TF名来自原数据库匹配，序列短名则来自你的矩阵；二者不是同一套名称，也不是一对一的TF身份确认。', '',
        '例如CN3GN对应F03-preTES-N04的CNGGGN短特征，原有数据库注释为CTCF，q=0.0409677。名字概括该pattern的序列偏好，CTCF是另一项数据库相似性注释。', '',
        '## 对应表的主要字段', '',
        '|字段|含义|', '|---|---|',
        '|sequence_name|44个原pattern各自唯一的展示名，包含必要的区分编号。|',
        '|feature_label / core_IUPAC_6bp|不含区分编号的特征短名，以及展开后的6bp简并共识。|',
        '|pattern_id|原pattern ID。如F03-Donor-P01：F03为第3折，Donor为供体区域，P为正归因，01为该区域方向下原pattern编号。N表示负归因。|',
        '|region|TSS、Donor、Acceptor、preTES四种输入区域，分别为转录起始位点区域、剪接供体区域、剪接受体区域、转录终止位点前区域。|',
        '|core_start_1based / core_end_1based|所取6bp在原50bp pattern矩阵中的起止位置，从1开始；不是基因组坐标。|',
        '|core_information_bits / mean_information_bits_per_position|六碱基窗口总信息量/每位置平均信息量，相对A/C/G/T各0.25背景计算；不是q值。|',
        '|core_expected_GC_fraction|根据6个位置的概率矩阵计算的平均G+C比例，不是直接数简并字母。|',
        '|composition_description|显示性成分描述。期望GC≥0.7为GC-rich，期望AT≥0.7为AT-rich，其余Mixed；纯同碱基共识另标A-run/T-run等。|',
        '|source_seqlets|形成原pattern的seqlet数量，不是独立生物重复数。|',
        '|best_TF / best_matrix_id / best_q|最低q记录对应的TF、参考矩阵ID及q值。q并列按矩阵编号排序，best_TF不表示唯一可能TF。|',
        '|all_TFs / all_matrix_ids / all_q_values|每个pattern全部已有的合格匹配，不只保留best_TF。|',
        '|full_50bp_IUPAC|整个50bp矩阵的简并共识，不等同于6bp短标签。|',
        '|reverse_complement_IUPAC / canonical_RC_signature|6bp简并共识的反向互补及用于核对方向等价性的规范化代码；没有改变原显示方向或重新聚类。|', '',
        'q值是该匹配被认定为显著时对应的最低假发现率阈值。当前保留原标准q≤0.05且重叠≥5bp。q越小，数据库相似性匹配的统计证据通常越强，但不是“这个TF真实结合概率为1-q”。见[Tomtom官方输出说明](https://meme-suite.org/meme/doc/tomtom-output-format.html)。', '',
        '## 当前文件索引与复现', '',
        '- `analysis/matched44_sequence_names.tsv`：44行完整命名与特征表。',
        '- `analysis/matched44_named_TOMTOM79.tsv`：79条原匹配附加新名称。',
        '- `analysis/matched44_sequence_naming_notes.md`：本说明及44个名称的逐一对应表。',
        '- `analysis/matched44_sequence_naming_QA.json`：原始H5、命名和图件核验记录。',
        '- `analysis/figure_cleanup_receipt.json`：本次删除范围、备份和文件哈希。',
        '- `name_and_plot_matched44.py`：命名和绘图复现入口；旧图输入依赖的QA记录仍保留。', '',
        '命名读取的24份原始H5均与注册SHA256一致，未重新训练模型、重新运行TF-MoDISco或新增TOMTOM匹配。原树SHA256：`'+report['topology_sha256']+'`。所有PDF页面保持5×4英寸，黑色Arial 12pt标签，无文字重叠或裁切，无嵌入位图。',
    ]
    guide = '\n'.join(lines) + '\n'
    write_text(ROOT / 'README_motifStack_重绘说明.md', guide)
    write_text(ROOT / 'analysis/current_figure_guide.md', guide)
    individual = ['', '## 44个名称逐一对应', '',
                  '|名称|原pattern|6bp共识|窗口（1-based）|最低q对应TF|q|', '|---|---|---|---|---|---|']
    for row in rows:
        individual.append('|'+ '|'.join(str(x) for x in [row['sequence_name'], row['pattern_id'],
                          row['core_IUPAC_6bp'], f"{row['core_start_1based']}-{row['core_end_1based']}",
                          row['best_TF'], format(row['best_q'], '.6g')])+'|')
    write_text(ROOT / 'analysis/matched44_sequence_naming_notes.md', guide + '\n'.join(individual) + '\n')


def main():
    figures = (ROOT / 'figures').resolve()
    report = json.loads((ROOT / 'analysis/matched44_sequence_naming_QA.json').read_text())
    keep = {}
    for key in ('full_tree', 'matched_tree'):
        path = checked_file(Path(report[key]['output']), figures)
        assert sha(path) == report[key]['pdf_sha256']
        keep[path.name] = report[key]['pdf_sha256']
    assert len(keep) == 2 and not set(keep) & set(OLD_PDFS)
    old = [checked_file(figures / name, figures) for name in OLD_PDFS if (figures / name).exists()]
    if not old:
        receipt = json.loads((ROOT / 'analysis/figure_cleanup_receipt.json').read_text())
        assert Path(receipt['recoverable_archive']).is_file()
        write_guide(report, receipt)
        print(json.dumps(dict(status='already_cleaned', receipt=receipt)))
        return
    assert len(old) == 10, 'Inventory changed: do not silently extend or partly delete the reviewed list.'
    now = datetime.now(ZoneInfo('Asia/Shanghai'))
    archive_dir = ROOT / 'archive'
    assert not archive_dir.is_symlink() and archive_dir.resolve().parent == ROOT
    archive_dir.mkdir(exist_ok=True)
    archive = archive_dir / ('previous_figures_' + now.strftime('%Y%m%d_%H%M%S') + '.zip')
    assert ROOT in archive.resolve().parents and not archive.exists()
    backups = old + [ROOT / 'README_motifStack_重绘说明.md', ROOT / 'analysis/matched44_sequence_naming_notes.md']
    manifest = {str(path.relative_to(ROOT)): dict(sha256=sha(path), bytes=path.stat().st_size) for path in backups}
    with zipfile.ZipFile(archive, 'x', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as handle:
        for path in backups:
            handle.write(path, str(path.relative_to(ROOT)))
        handle.writestr('MANIFEST_SHA256.json', json.dumps(manifest, indent=2))
    with zipfile.ZipFile(archive) as handle:
        assert handle.testzip() is None
        for relative, record in manifest.items():
            assert hashlib.sha256(handle.read(relative)).hexdigest() == record['sha256']
    # Only delete the ten explicitly reviewed old PDFs after verified backup.
    # No recursive operation, directory deletion, or unreviewed path is used.
    deleted = []
    for path in old:
        checked_file(path, figures)
        assert sha(path) == manifest[str(path.relative_to(ROOT))]['sha256']
        path.unlink()
        deleted.append(str(path))
    assert all(not Path(path).exists() for path in deleted)
    assert {path.name for path in figures.glob('*.pdf')} == set(keep)
    for name, expected in keep.items():
        assert sha(figures / name) == expected
    assert sha(ROOT / 'analysis/motifStack_hclust.rds') == report['topology_sha256']
    assert sha(ROOT / 'analysis/matched44_sequence_names.tsv') == report['table_sha256']
    receipt = dict(status='PASS', time_Asia_Shanghai=now.isoformat(),
                   deleted_old_pdf_count=len(deleted), deleted_old_pdf_paths=deleted,
                   recoverable_archive=str(archive), archive_sha256=sha(archive),
                   archive_verified=True, current_pdfs=keep, current_pdf_count=2,
                   current_pdfs_unchanged=True, naming_table_unchanged=True,
                   original_tree_unchanged=True, original_data_deleted=False,
                   backup_manifest=manifest)
    write_text(ROOT / 'analysis/figure_cleanup_receipt.json', json.dumps(receipt, indent=2))
    write_guide(report, receipt)
    print(json.dumps({k:v for k,v in receipt.items() if k not in ('backup_manifest','deleted_old_pdf_paths')}))
    print('GUIDE',str(ROOT / 'analysis/current_figure_guide.md'))


if __name__ == '__main__':
    main()
