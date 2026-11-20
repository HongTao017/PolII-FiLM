"""Reuse the original top-30 horizontal-bar figure with audited REG5 data."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
REFERENCE = Path('G:/sequence_analysis_data/redraw_gene_features_20261021')
sys.path.insert(0, str(REFERENCE))
import plot_rf_importance_reference as reference

FINGERPRINT = 'c13b03947d1c64bd21ad4525860c8a1528cb98ff6d9d711936451f360cbaf0a4'
DEFAULT_OUTPUT = Path('G:/sequence_analysis_figures/Png/png/biochemical_feature_figures/fig3b/fig3b_RF_BChIPeCLIP_top30_importance_reference.pdf')


def configure_reference():
    # Only data identity/schema changes in the read-only source audit.
    source = reference.REMOTE_READ
    source = source.replace('/mnt/HDD8TB/hongtao/Data_Seq/model/BQICPET/new_train',
                            '/mnt/HDD8TB/hongtao/Data_Seq/model/B_REG5_RF_20261027')
    source = source.replace("==list('BQIP')", "==['B','ChIP','eCLIP']")
    source = source.replace('BQIP', 'BChIPeCLIP').replace('383', '2512')
    source = source.replace("assert complete['status']=='complete'",
        "assert complete['status']=='complete'\nassert fingerprint=='" + FINGERPRINT + "'\n"
        "assert complete['ordered_selected_groups']==complete['greedy_selected_groups']==['B','ChIP','eCLIP']\n"
        "assert json.loads((root/'final_source_preservation.json').read_text())['status']=='PASS'")
    reference.REMOTE_READ = source
    reference.FINGERPRINT = FINGERPRINT


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument('--audit-only', action='store_true')
    args = parser.parse_args()
    source_file = REFERENCE / 'plot_rf_importance_reference.py'
    original_hash = hashlib.sha256(source_file.read_bytes()).hexdigest()
    configure_reference()
    packet = reference.retrieve()
    assert packet['provenance']['feature_count'] == 2512
    assert len(packet['provenance']['fit_npz_sha256']) == 50
    if args.audit_only:
        print(json.dumps({'provenance': {k: v for k,v in packet['provenance'].items() if k not in ['fit_npz_sha256','source_sha256']},
                          'top30': packet['top30']}, ensure_ascii=True))
        return
    for row in packet['top30']:
        feature = row['feature']
        group, base = feature.split('::', 1)
        if group == 'B':
            label = reference.LABELS.get(base, base.replace('_', ' '))
        elif group == 'ChIP':
            cell_target, region = base.split('_ChIP_seq_', 1)
            cell, target_name = cell_target.split('_', 1)
            region_label = {'gene_body_count': 'gene-body count', 'pause_count': 'pause count'}.get(region, region.replace('_', ' '))
            label = f'{cell} {target_name} ({region_label})'
        else:
            label = feature.replace('::', ': ')
        reference.LABELS[feature] = label
    # Literal original bar-figure function, with data-dependent title only.
    fig, qa = reference.build_figure(packet)
    ax = fig.axes[0]
    ax.set_title('RF model (B + ChIP-seq + eCLIP)', fontsize=14, fontweight='bold', pad=10)
    qa = reference.check_layout(fig, ax)
    from matplotlib.text import Text
    texts = [t for t in fig.findobj(Text) if t.get_visible() and t.get_text().strip()]
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for i, text in enumerate(texts):
        first = text.get_window_extent(renderer)
        for other in texts[i+1:]:
            second = other.get_window_extent(renderer)
            assert min(first.x1, second.x1)-max(first.x0,second.x0) <= 1 or min(first.y1,second.y1)-max(first.y0,second.y0) <= 1, ('text overlap',text.get_text(),other.get_text())
    target = args.output
    if target.exists():
        raise FileExistsError(target)
    assert target.suffix.lower() == '.pdf' and target.parent.is_dir()
    report = {'qa':qa, 'provenance':packet['provenance'], 'top30':packet['top30'],
              'reference_source':str(source_file), 'reference_source_sha256':original_hash,
              'adaptations':'REG5 input source, feature labels, model title only; original bar drawing function and palette unchanged. Regulatory labels are cell target (region count); all regulatory bars in Top30 are ChIP-seq.',
              'display_labels': {r['feature']: reference.LABELS[r['feature']] for r in packet['top30']},
              'interpretation':'Unsigned impurity-based importance; not causal or signed effects; Top30 not renormalized'}
    raw_save = fig.savefig
    def save_with_metadata(path, **kwargs):
        return raw_save(path, metadata={'Title':'RF B + ChIP-seq + eCLIP: top30 feature importance',
                                        'Subject':json.dumps(report,ensure_ascii=True)}, **kwargs)
    fig.savefig = save_with_metadata
    reference.save_pdf(fig, target)
    assert hashlib.sha256(source_file.read_bytes()).hexdigest() == original_hash
    print(json.dumps({'status':'PASS','file':str(target),'bytes':target.stat().st_size,
                      'qa':qa,'top30':[{'feature':r['feature'],'mean_percent':100*r['mean']} for r in packet['top30']]},ensure_ascii=True))


if __name__ == '__main__':
    main()
