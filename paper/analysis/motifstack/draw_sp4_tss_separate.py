"""Standalone SP4 TSS overview and native pattern CWM; PDF-only server output."""
from pathlib import Path
import sys, json, hashlib
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager, transforms
from matplotlib.textpath import TextPath
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.ticker import MaxNLocator, FuncFormatter
ROOT = Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019')
A = ROOT/'analysis/SP4_F04_r0_neg_pattern0'
F = ROOT/'figures'
sys.path.insert(0, str(ROOT))
from final_figure_small import apply_style, save_pdf
for p in (ROOT/'fonts').glob('*.ttf'):
    font_manager.fontManager.addfont(str(p))
apply_style()
info = json.loads((A/'selection.json').read_text())
old_qa = json.loads((A/'full_input_CWM_QA.json').read_text())
def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(2**20), b''):
            h.update(block)
    return h.hexdigest()
old = F/'SP4_MA0685_2_F04_FULL_INPUT_PATTERN_CWM.pdf'
assert old.resolve().parent == F.resolve()
if old.exists():
    assert digest(old) == old_qa['sha256'], 'Old SP4 PDF changed; retain rather than delete'
preserved = {p.name: digest(p) for p in F.glob('*.pdf') if p != old and p.name not in ('SP4_MA0685_2_F04_TSS600_IMPORTANCE.pdf', 'SP4_MA0685_2_F04_PATTERN_CWM.pdf')}
with np.load(A/'selected_real_attribution.npz') as data:
    actual = data['actual'].copy()
    cwm = data['pattern_cwm'].copy()
with np.load(A/'full_input_and_pattern_cwm.npz') as data:
    assert np.allclose(actual, data['actual'][:600], atol=1e-9, rtol=1e-7)
    assert np.array_equal(cwm, data['pattern_cwm'])
h5path = Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motif/folds/fold_04/modisco/tss_sequence.h5')
assert digest(h5path) == info['source_sha256'][str(h5path)]
with h5py.File(h5path, 'r') as h5:
    group = h5['neg_patterns/pattern_0']
    assert np.array_equal(cwm, group['contrib_scores'][:])
    assert int(group['seqlets/n_seqlets'][0]) == 166
assert actual.shape == (600, 4) and cwm.shape == (50, 4)
y = actual.sum(axis=1)
s, e = info['selected']['start'], info['selected']['end']
a, b = info['hit']['motif_start'], info['hit']['motif_end']
assert (s, e, a, b, info['tss_input_index']) == (280, 330, 296, 305, 300)

def format_axes(ax, ylabel):
    ax.spines[['top', 'right']].set_visible(False)
    ax.spines[['left', 'bottom']].set_linewidth(.55)
    ax.axhline(0, color='#B5B5B5', lw=.65, ls=(0, (5, 4)), zorder=1)
    ax.tick_params(labelsize=12, width=.55, length=3, pad=5)
    ax.set_ylabel(ylabel, fontsize=14, labelpad=10)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: '0' if v == 0 else f'{v:+.4f}'.rstrip('0').rstrip('.')))

def set_inside_yticks(ax):
    lo, hi = ax.get_ylim()
    step = .0005 if max(abs(lo), abs(hi)) <= .0013 else .001
    ax.set_ylim(-2.25*step, 2.25*step)
    ax.set_yticks(np.array([-2, -1, 0, 1, 2])*step)
def check_figure(fig):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    boxes = []
    for txt in fig.findobj(matplotlib.text.Text):
        if not txt.get_visible() or not txt.get_text():
            continue
        bb = txt.get_window_extent(renderer)
        assert bb.x0 >= -.5 and bb.y0 >= -.5 and bb.x1 <= fig.bbox.width+.5 and bb.y1 <= fig.bbox.height+.5, (txt.get_text(), bb.bounds)
        boxes.append({'text': txt.get_text(), 'bounds_px': list(bb.bounds)})
    # Axis labels and neighbouring tick labels must have distinct boxes.
    for ax in fig.axes:
        for ticks in (ax.get_xticklabels(), ax.get_yticklabels()):
            visible = [t for t in ticks if t.get_visible() and t.get_text()]
            for left, right in zip(visible, visible[1:]):
                assert not left.get_window_extent(renderer).overlaps(right.get_window_extent(renderer)), (left.get_text(), right.get_text())
    return boxes

qa = {}
fig = plt.figure(figsize=(12, 3.5), facecolor='white')
fig.text(.11, .95, 'SP4-like motif · MA0685.2', fontsize=12, va='top')
fig.text(.11, .865, 'KLHL23 · Fold 04 · complete TSS window (600 bp)', fontsize=12, va='top')
ax = fig.add_axes([.11, .36, .855, .38])
format_axes(ax, 'Importance\nscore')
ax.spines['bottom'].set_visible(False)
ax.tick_params(axis='x', length=0)
# All 600 values are retained; horizontal expansion replaces smoothing/subsampling.
line, = ax.plot(np.arange(600), y, color='#727272', lw=.6, zorder=4)
lim = float(np.abs(y).max())*1.16
ax.set_ylim(-lim, lim)
ax.set_xlim(-.5, 599.5)
ax.set_xticks([0, 150, 300, 450, 599], ['0', '150', '300', '450', '599'])
set_inside_yticks(ax)
ax.set_xlabel('Position in TSS input (bp)', fontsize=14, labelpad=26)
ax.axvspan(s-.5, e-.5, color='#D4D4D4', alpha=.45, zorder=0)
# Native pattern source seqlet and exact SP4 hit are distinct nested intervals.
ax.add_patch(Rectangle((s-.5, -.32), e-s, .075, transform=ax.get_xaxis_transform(), color='#BFBFBF', lw=0, clip_on=False))
ax.add_patch(Rectangle((a-.5, -.32), b-a, .075, transform=ax.get_xaxis_transform(), color='#787878', lw=0, clip_on=False))
ax.plot([300, 300], [0, 1.07], transform=ax.get_xaxis_transform(), color='#ADADAD', lw=.65, clip_on=False, zorder=2)
ax.text(300, 1.15, 'TSS', transform=ax.get_xaxis_transform(), fontsize=12, ha='left', va='bottom')
ax.annotate('', xy=(322, 1.08), xytext=(300, 1.08), xycoords=ax.get_xaxis_transform(), arrowprops={'arrowstyle': '-|>', 'color': '#4F4F4F', 'lw': .65, 'mutation_scale': 8}, annotation_clip=False)
fig.text(.11, .035, 'Grey: pattern source seqlet (50 bp); dark grey: SP4-like match (9 bp)', fontsize=12, va='bottom')
assert np.array_equal(line.get_ydata(), y) and len(line.get_xdata()) == 600
boxes = check_figure(fig)
path = F/'SP4_MA0685_2_F04_TSS600_IMPORTANCE.pdf'
save_pdf(fig, path)
qa['TSS'] = {'path': str(path), 'canvas_inches': [12, 3.5], 'point_count': 600, 'raw_scores_preserved': True, 'smoothing': False, 'subsampling': False, 'signed_range': [float(y.min()), float(y.max())], 'seqlet_input_half_open': [s, e], 'motif_input_half_open': [a, b], 'TSS_index0': 300, 'text_boxes': boxes}
plt.close(fig)

fig = plt.figure(figsize=(12, 3.8), facecolor='white')
fig.text(.11, .95, 'SP4-like motif · MA0685.2', fontsize=12, va='top')
fig.text(.11, .865, 'Pattern CWM · 166 seqlets · Fold 04 / r0_neg_patterns_pattern_0', fontsize=12, va='top')
ax = fig.add_axes([.11, .23, .855, .54])
format_axes(ax, 'Contribution\nscore')
font = font_manager.FontProperties(family='Arial', weight='bold')
colors = {'A': '#719879', 'C': '#658DB7', 'G': '#C5A453', 'T': '#C7767D'}
patches = []
geometry = np.zeros_like(cwm)
for pos, row in enumerate(cwm):
    for positive in (False, True):
        stack = 0.
        indices = [k for k in range(4) if row[k] != 0 and (row[k] > 0) == positive]
        for k in sorted(indices, key=lambda j: abs(row[j])):
            value = float(row[k]); height = abs(value)
            glyph = TextPath((0, 0), 'ACGT'[k], size=1, prop=font)
            box = glyph.get_extents()
            base_y = stack if positive else -stack-height
            trans = transforms.Affine2D().translate(-box.x0, -box.y0).scale(.88/box.width, height/box.height).translate(pos+1-.44, base_y)
            patch = PathPatch(glyph, transform=trans+ax.transData, facecolor=colors['ACGT'[k]], edgecolor='none', zorder=4)
            ax.add_patch(patch); patches.append(patch); stack += height; geometry[pos, k] = value
lo = min(float(np.minimum(cwm, 0).sum(1).min()), 0.)
hi = max(float(np.maximum(cwm, 0).sum(1).max()), 0.)
span = max(hi-lo, 1e-9)
ax.set_ylim(lo-.13*span, hi+.18*span)
ax.set_xlim(.5, 50.5)
ax.set_xticks([1, 10, 20, 30, 40, 50])
set_inside_yticks(ax)
ax.set_xlabel('Pattern position (bp)', fontsize=14, labelpad=9)
assert np.array_equal(geometry, cwm), 'Vector glyph heights must preserve all signed native values'
boxes = check_figure(fig)
path = F/'SP4_MA0685_2_F04_PATTERN_CWM.pdf'
save_pdf(fig, path)
qa['CWM'] = {'path': str(path), 'canvas_inches': [12, 3.8], 'matrix_shape': [50, 4], 'n_seqlets': 166, 'vector_base_letters': len(patches), 'native_H5_values_preserved': True, 'positive_negative_stacks_preserved': True, 'text_boxes': boxes}
plt.close(fig)

sys.path.insert(0, str(ROOT/'qa_deps'))
import pymupdf
for result in qa.values():
    path = Path(result['path'])
    pdf = pymupdf.open(path)
    assert len(pdf) == 1
    page = pdf[0]
    assert np.allclose([page.rect.width, page.rect.height], np.array(result['canvas_inches'])*72, atol=.02)
    assert not page.get_images()
    sizes = sorted({round(s['size'], 2) for b0 in page.get_text('dict')['blocks'] if 'lines' in b0 for l in b0['lines'] for s in l['spans']})
    assert sizes == [12., 14.], sizes
    assert all('Arial' in f[3] for f in page.get_fonts())
    pix = page.get_pixmap(matrix=pymupdf.Matrix(2, 2))
    pixels = np.frombuffer(pix.samples, dtype=np.uint8)
    assert pixels.min() < 90 and (pixels < 200).mean() > .001
    result.update(sha256=digest(path), bytes=path.stat().st_size, pages=1, font_sizes=sizes, fonts=page.get_fonts(), no_raster_images=True, text_within_canvas=True, neighbouring_tick_labels_do_not_overlap=True, rendered_on_server_in_memory=True)
    pdf.close()
# Only after BOTH final PDFs pass QA, remove the exact superseded SP4 combined file.
receipt = {'new_figures': [qa[k]['path'] for k in qa], 'removed': [], 'preserved_other_pdfs': preserved}
if old.exists():
    receipt['removed'].append({'path': str(old), 'sha256': digest(old), 'bytes': old.stat().st_size})
    old.unlink()
assert not old.exists()
assert preserved == {p.name: digest(p) for p in F.glob('*.pdf') if p.name in preserved}
qa['source'] = {'pattern': info['pattern'], 'gene': info['selected']['gene'], 'H5': str(h5path), 'H5_SHA256': digest(h5path), 'CWM_is_pattern_aggregate_not_gene_zoom': True, 'grey_track_is_single_KLHL23_TSS_actual_attribution': True, 'no_result_download': True, 'user_authorized_size_override': True}
(A/'TSS_AND_CWM_SEPARATE_QA.json').write_text(json.dumps(qa, indent=2))
(A/'TSS_AND_CWM_CLEANUP_RECEIPT.json').write_text(json.dumps(receipt, indent=2))
(A/'TSS_AND_CWM_README.md').write_text('''# Current separate SP4 figures\n\n1. ../../figures/SP4_MA0685_2_F04_TSS600_IMPORTANCE.pdf — KLHL23 original complete TSS window, 600 observed-base signed attribution values (five-seed mean). Light grey marks source seqlet input [280,330), TSS-relative [-20,+29]; dark grey underbar marks SP4-like sequence match [296,305), TSS-relative [-4,+4]. TSS is input index 300. No smoothing, downsampling, absolute value, or sign flip. This is one real sequence instance, not a pattern average.\n\n2. ../../figures/SP4_MA0685_2_F04_PATTERN_CWM.pdf — original native 50x4 TF-MoDISco contribution weight matrix for F04__r0_neg_patterns_pattern_0, 166 seqlets. All signed A/C/G/T values retained; vector letter heights stack above/below zero. Not a numerical zoom of the first PDF.\n\nGrey highlights identify the selected region; other input branches are not shown. SP4-like denotes sequence similarity to MA0685.2, not validated SP4 occupancy; SP1/SP2 similarities also exist. Target unit: signed predicted Y_PC1 attribution. PDF only, editable Arial12/14 text; user-authorized 12x3.5 and 12x3.8 inch canvases. Both generated and checked directly on server, no raster files or result downloads.\n\nThe superseded full-input combined SP4 figure was deleted after both replacements passed QA. Circular motif trees and source files were preserved. See TSS_AND_CWM_SEPARATE_QA.json and TSS_AND_CWM_CLEANUP_RECEIPT.json.\n''')
print(json.dumps({k: {key: value for key, value in result.items() if key not in ('text_boxes', 'fonts')} for k, result in qa.items()}, indent=2))
print('CLEANUP', json.dumps(receipt))

