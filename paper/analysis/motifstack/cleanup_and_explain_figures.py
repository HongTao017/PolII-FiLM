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
        "# Figure and sequence-label guide", '',
        "The 44 matched sequence patterns have unique display names representing 29 distinct six-base sequence features. These labels do not establish 29 validated biological motif families. The 79 qualifying TOMTOM matches and 11 JASPAR matrices are retained.", '',
        "## Superseded figure versions", '',
        f"Removed {receipt['deleted_old_pdf_count']} historical PDFs after preserving compressed backups and verifying their recorded SHA256 values. The following two final variants remain:", '',
        '`' + receipt['recoverable_archive'] + '`', '',
        "Original H5 inputs, tree topology, distances, naming tables, scripts and required QA records are retained. Cleanup concerns the listed superseded PDFs and preserves experiment and analysis inputs.", '',
        "## Current figures", '',
        "| File | Content |", '|---|---|',
        f'|`{current[0]}`|Overview of all 314 patterns, with sequence-feature labels for 11 representatives. These labels do not imply only 11 patterns or sequence groups.|',
        f'|`{current[1]}`|Two-page detail view of the same 44-member matched subtree, with 22 alternating names per page and each name shown once.|', '',
        "The detailed view prunes 270 unmatched leaves from the original 314-member tree and contracts single-child paths. It preserves branching among the 44 matched patterns without reclustering. All leaves in this view are matched, so its outer ring is dark purple.", '',
        "## Colors and symbols", '',
        "| Element | Meaning |", '|---|---|',
        "| Dark-purple outer ring (#918BB4) | Patterns with qualifying JASPAR/TOMTOM matches: 44 in the overview. |",
        "| Light-purple outer ring (#DDDDEA) | Patterns below the match threshold: 270 in the overview. This does not establish absence of regulatory function. |",
        "| Gray sectors/background (#B0B0B0) | Regions containing matched sequence patterns. |",
        "| Pink branches | Positive-PC1 attribution patterns or subtrees containing only positive patterns. |",
        "| Blue branches | Negative-PC1 attribution patterns or subtrees containing only negative patterns. |",
        "| Neutral-gray internal branches | Subtrees containing both attribution signs; distinct from gray background sectors. |",
        "| Black triangles and short lines | Links from sequence labels to pattern tips; they do not indicate transcription direction or TF activation/repression. |",
        "| Branching | Original sequence-similarity hierarchical clustering. Display angles, lengths and sector widths support layout and are not quantitative enrichment measures. |", '',
        "Attribution signs describe contributions to predicted Y_PC1, not experimentally measured expression changes.", '',
        "## Interpreting sequence labels", '',
        "Labels summarize the highest-information contiguous six-base window of the original full-precision 50-bp probability matrix using IUPAC consensus. They are short sequence features, not complete patterns or guaranteed TOMTOM alignment intervals.", '',
        "Runs of at least three identical bases use count-plus-base notation: 6T=TTTTTT, 5TA=TTTTTA and CN3GN=CNGGGN. Suffixes distinguish patterns sharing a short label; they do not indicate significance or separate TF families.", '',
        "| IUPAC code | Allowed bases |", '|---|---|',
        "| N | A/C/G/T when the consensus rule cannot use a more specific code; not missing input sequence. |",
        "| R | A or G |", "| Y | C or T |", "| M | A or C |", "| S | C or G |", "| W | A or T |", "| K | G or T |", '',
        "Allowed bases need not have equal probabilities. Full probabilities remain in the source matrices. IUPAC definitions and the consensus rule follow the documented Biopython implementation.", '',
        "## The 29 short sequence features across 44 patterns", '',
        "| Label | Expanded six-base feature | Pattern count | Retained TF matches |", '|---|---|---:|---|',
    ]
    for stem, members in groups.items():
        assert len({x['core_IUPAC_6bp'] for x in members}) == 1
        display = members[0]['sequence_name'] if len(members) == 1 else members[0]['sequence_name'] + '～' + members[-1]['sequence_name']
        tfs = sorted({tf for x in members for tf in x['all_TFs'].split(';')})
        lines.append('|'+ '|'.join([display, members[0]['core_IUPAC_6bp'], str(len(members)), '、'.join(tfs)])+'|')
    lines += [
        '', "TF labels are database-similarity annotations; short sequence labels summarize the matrices. Neither establishes a unique TF identity.", '',
        "For example, CN3GN summarizes the CNGGGN core of F03-preTES-N04, whose retained database annotation is CTCF with q=0.0409677.", '',
        "## Main table fields", '',
        "| Field | Description |", '|---|---|',
        "| sequence_name | Unique display name for each of the 44 patterns. |",
        "| feature_label / core_IUPAC_6bp | Shared short feature label and expanded six-base consensus. |",
        "| pattern_id | Original fold, region, attribution sign and within-region index. F03-Donor-P01 denotes fold 3, donor region, positive attribution, pattern 01. |",
        "| region | TSS, Donor, Acceptor or preTES input region. |",
        "| core_start_1based / core_end_1based | One-based core coordinates within the 50-bp pattern matrix, not genomic coordinates. |",
        "| core_information_bits / mean_information_bits_per_position | Core information relative to a uniform base background, not a q value. |",
        "| core_expected_GC_fraction | Mean expected G+C fraction from the six position-probability rows. |",
        "| composition_description | Display annotation: GC-rich or AT-rich at an expected fraction of at least 0.7; otherwise Mixed. Homopolymer cores receive run labels. |",
        "| source_seqlets | Seqlets contributing to the pattern, not independent biological replicates. |",
        "| best_TF / best_matrix_id / best_q | Lowest-q retained annotation, with matrix-ID tie breaking; not a unique TF assignment. |",
        "| all_TFs / all_matrix_ids / all_q_values | All qualifying database matches. |",
        "| full_50bp_IUPAC | Consensus across the full 50-bp matrix. |",
        "| reverse_complement_IUPAC / canonical_RC_signature | Reverse complement and orientation-equivalence signature; original display orientation and clustering are retained. |", '',
        "Retained matches use q<=0.05 and overlap>=5 bp. A q value describes false-discovery-rate significance, not a TF-binding probability.", '',
        "## Files and reproduction", '',
        "- `analysis/matched44_sequence_names.tsv`: full naming table for 44 patterns.",
        "- `analysis/matched44_named_TOMTOM79.tsv`: 79 retained matches with display names.",
        "- `analysis/matched44_sequence_naming_notes.md`: sequence-label documentation.",
        "- `analysis/matched44_sequence_naming_QA.json`: source, naming and figure records.",
        "- `analysis/figure_cleanup_receipt.json`: cleanup scope and file hashes.",
        "- `name_and_plot_matched44.py`: sequence naming and plotting entry point.", '',
        "The 24 source H5 files retain their registered hashes. Naming does not retrain models, rerun motif discovery or add database matches. Original tree SHA256: `"+report['topology_sha256']+"`. PDF pages retain the 5 by 4 inch canvas and editable text.",
    ]
    guide = '\n'.join(lines) + '\n'
    write_text(ROOT / "README_motifStack_redraw.md", guide)
    write_text(ROOT / 'analysis/current_figure_guide.md', guide)
    individual = ['', "## Pattern-name correspondence", '',
                  "| Label | Pattern | Six-base consensus | One-based core window | Lowest-q annotation | q |", '|---|---|---|---|---|---|']
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
    backups = old + [ROOT / "README_motifStack_redraw.md", ROOT / 'analysis/matched44_sequence_naming_notes.md']
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
