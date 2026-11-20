"""Record reproducible local Fig5I package metadata; does not run any models."""
import hashlib
import json
from pathlib import Path
import shutil
import sys

ROOT = Path(sys.argv[1]).resolve()


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


if '--manifest-only' not in sys.argv:
    snap = ROOT / 'provenance/native_analysis_snapshot'
    shutil.copy2(snap / 'paired_effects.tsv', ROOT / 'data/full_results/paired_effects.tsv')
    protocol = json.loads((snap / 'protocol.json').read_text(encoding='utf-8'))
    run = ROOT / 'provenance/upstream_code/run_native.py'
    assert digest(run) == protocol['script_sha256']
    for fold in range(1, 11):
        assert digest(snap / f'plans/fold_{fold:02}.jsonl') == protocol['plans'][str(fold)]
    source_manifest = json.loads((ROOT / 'provenance/source_plot_manifest.json').read_text(encoding='utf-8'))
    meta = dict(panel='I', filename='fig5I.pdf', entry_script='code/render.py',
        description='Four selected native-motif disruption raincloud groups with matched Background; no stars',
        kind='native_motif_knockdown_matched_raincloud', canvas_inches=[8, 2.5],
        y_limits=[-.002,.002], tick_step=.001, display_power=-3, significance_stars_displayed=False,
        labels=['KLF16','SP1','SP2','SP4','Background'], points=[335,223,268,238,389],
        total_points=1453, no_subsampling=True, delta='WT_minus_disrupted', transform='none',
        background='Controls of displayed groups: deduplicate match_id then average by source positive DNA overlap group',
        source_local_pdf='G:/sequence_analysis_figures/Png/png/fig5/fig5h/native_motif_raincloud_20261112/raincloud_selected_y002_no_stars_20261112/selected_four_y_minus2_to2_no_stars.pdf',
        legacy_server_directory='/mnt/HDD8TB/hongtao/Data_Seq/task/A_FOUR_REGION_ISM_20261108/native_motif_knockdown_matched_20261112',
        source_pdf_sha256=digest(ROOT / 'fig5I.pdf'), groups=source_manifest['groups'],
        source_data_sha256={p.relative_to(ROOT).as_posix():digest(p) for p in sorted((ROOT / 'data').rglob('*')) if p.is_file()},
        reproduce='python plot.py --output audit/redraw.pdf', no_network_for_redraw=True,
        code_adaptations=['Move frozen input paths under data/ and drawing implementation under code/',
            'Force DejaVu Sans matching the embedded source fonts on Windows',
            'Remove unused statistics reads and old server-only audit code; preserve all drawing parameters and seeds'],
        requires_external_model_for_inference=True, new_inference=False, new_statistics=False)
    (ROOT / 'package.json').write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding='utf-8')
    (ROOT / 'provenance/source_record.json').write_text(json.dumps(dict(
        local_source=meta['source_local_pdf'], server_source=meta['legacy_server_directory'],
        source_pdf_byte_identical=True, source_pdf_sha256=meta['source_pdf_sha256'],
        frozen_inference_plans_and_per_seed_predictions=True, upstream_script_hash_matches_protocol=True,
        source_genes=2000, full_background_universe=15351, source_patterns=44,
        matched_pairs=1368, unique_positive_genes=994, main_labels=meta['labels'],
        does_not_include_model_checkpoints_or_original_training_inputs=True), ensure_ascii=False, indent=2), encoding='utf-8')
    (ROOT / 'audit/requirements.txt').write_text('numpy\npandas\nPillow\npypdf\n# Poppler pdftoppm must be on PATH for pixel comparison.\n', encoding='utf-8')
    (ROOT / 'provenance/upstream_code/README.md').write_text(
        "These archived scripts retain their experiment paths and scheduling parameters as provenance; frozen-data plotting does not execute them."
        "The previous source records include protocol and ten-fold plan hashes."
        "Re-inference requires the recorded input data, checkpoints and model loader. Use the plotting entry point for local frozen-data figures."
        "plot.py reads frozen numeric values and does not require a GPU.\n", encoding='utf-8')

inventory = []
for path in sorted(ROOT.rglob('*')):
    if path.is_file() and path.name != 'MANIFEST_SHA256.txt':
        inventory.append(f'{digest(path)}  {path.relative_to(ROOT).as_posix()}')
(ROOT / 'MANIFEST_SHA256.txt').write_text('\n'.join(inventory)+'\n', encoding='utf-8')
print(json.dumps(dict(directory=str(ROOT), hashed_files=len(inventory)), ensure_ascii=False))
