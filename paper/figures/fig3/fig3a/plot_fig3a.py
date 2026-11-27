
from pathlib import Path
import json, math, hashlib
import numpy as np
import final_figure_small as style

style.apply_style()
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Rectangle, Ellipse, Polygon, PathPatch
from matplotlib.path import Path as MPath

ROOT = Path(__file__).resolve().parent
META = json.loads((ROOT / 'verified_features.json').read_text(encoding='utf-8'))
OUT = ROOT / 'fig3a_large_text.pdf'
W, H = 1000, 520
for font in ['arial.ttf', 'arialbd.ttf']:
    font_manager.fontManager.addfont(str((ROOT / 'fonts') / font))
fig = plt.figure(figsize=(W / 72, H / 72), dpi=144, facecolor='white')
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W); ax.set_ylim(0, H); ax.set_aspect('equal'); ax.axis('off')
INK = '#292E35'
TEXT_COLOR = '#000000'
CHIP, CLIP, ATAC, DAM, PET = '#C78593', '#9277AF', '#559AB5', '#BB985C', '#579E95'
texts, leader_segments = [], []


def text(x, y, s, size=16, bold=False, anchor='center', color=INK):

    obj = ax.text(x, y, s, fontsize=size, fontfamily='Arial',
                  fontweight='bold' if bold else 'normal', color=TEXT_COLOR,
                  ha=anchor, va='baseline', zorder=10)
    texts.append(obj)
    return obj


def line(x1, y1, x2, y2, color=INK, width=1, dash=None):
    obj, = ax.plot([x1, x2], [y1, y2], color=color, lw=width,
                   solid_capstyle='butt', zorder=3)
    if dash: obj.set_dashes(dash)
    return obj


def rect(x, y, w, h, fill, edge=None, zorder=2):
    obj = Rectangle((x, y), w, h, facecolor=fill, edgecolor=edge or 'none',
                    lw=.8, zorder=zorder)
    ax.add_patch(obj)
    return obj


def ellipse(x, y, rx, ry, fill, edge=None, zorder=5):
    obj = Ellipse((x, y), 2 * rx, 2 * ry, facecolor=fill,
                  edgecolor=edge or 'none', lw=.8, zorder=zorder)
    ax.add_patch(obj)
    return obj


def arrow(x1, y1, x2, y2, color=INK, width=1.1):
    line(x1, y1, x2, y2, color, width)
    a = math.atan2(y2 - y1, x2 - x1); length = 6
    ax.add_patch(Polygon([(x2, y2),
        (x2 - length * math.cos(a - .45), y2 - length * math.sin(a - .45)),
        (x2 - length * math.cos(a + .45), y2 - length * math.sin(a + .45))],
        closed=True, facecolor=color, edgecolor='none', zorder=4))


def leader(points, color=INK):
    for a, b in zip(points[:-1], points[1:]):
        line(*a, *b, color, .9)
        leader_segments.append((a, b))
    ellipse(*points[-1], 2, 2, color, zorder=7)


START, TSS, BOUNDARY, TES = 190, 230, 360, 930
groups = META['group_counts']


text(145, 402, f"ATAC-seq ({groups['ATAC']:,})", size=18, bold=True, color=ATAC)
text(145, 382, 'HepG2 / K562')
text(145, 365, 'Accessible DNA', color='#5D6770')
text(365, 402, f"ChIP-seq ({groups['ChIP']:,})", size=18, bold=True, color=CHIP)
text(365, 382, 'MAFK / MAFF / ATF2')
text(365, 365, 'H2AFZ / H3K9me3')
text(625, 402, f"eCLIP ({groups['eCLIP']:,})", size=18, bold=True, color=CLIP)
text(625, 382, 'PRPF8 / PRPF4 / SF3B4')
text(625, 365, 'RNA-binding proteins', color='#5D6770')
text(860, 402, f"DamID ({groups['DamID']:,})", size=18, bold=True, color=DAM)
text(860, 382, 'CTCF / LMNB1 / NUP153')
text(860, 365, 'Protein-associated DNA intervals', color='#5D6770')

Y = 250
rect(START, Y - 14, BOUNDARY - START, 28, '#FAF0F1', zorder=0)
rect(BOUNDARY, Y - 14, TES - BOUNDARY, 28, '#F0F6F8', zorder=0)
line(76, Y, 967, Y, INK, 1.4)
for x in [428, 629, 808]: arrow(x, Y, x + 15, Y, '#5C6570', 1)
rect(230, Y - 5, 28, 10, '#BAD4DF', '#5E899B', 4)
for x, width in [(258, 40), (475, 54), (666, 66), (852, 39)]:
    rect(x, Y - 9, width, 18, '#76A6B8', '#5E899B', 4)
rect(891, Y - 5, 39, 10, '#B7CF9F', '#72905B', 4)
line(TSS, Y - 15, TSS, 287); arrow(TSS, 287, TSS + 28, 287)
text(TSS, 302, 'TSS', bold=True)
line(TES, Y - 15, TES, 282); text(TES, 298, 'TES', bold=True)
for x, label in [(131, 'Promoter'), (244, "5' UTR"), (395, 'Intron'),
                 (502, 'CDS'), (911, "3' UTR")]:
    text(x, 224, label)
line(BOUNDARY, Y - 16, BOUNDARY, Y + 16, '#A4B7BF', .8, [2, 2])


for x in [159, 174, 189, 204]: arrow(x, Y + 27, x, Y + 3, ATAC, .9)
leader([(145, 352), (145, 320), (174, 266)], ATAC)


for x in [319, 568]:
    rect(x - 13, Y + 2, 26, 5, CHIP, zorder=4)
    ellipse(x, Y + 17, 17, 10, '#E5C1C9', CHIP)
leader([(365, 352), (365, 328), (319, 272)], CHIP)


vertices = [(550, Y + 11), (540, 287), (560, 304), (580, 313)]
codes = [MPath.MOVETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4]
for x in [580, 610, 640, 670]:
    vertices.extend([(x + 8, 323), (x + 17, 302), (x + 30, 313)])
    codes.extend([MPath.CURVE4] * 3)
ax.add_patch(PathPatch(MPath(vertices, codes), facecolor='none', edgecolor=CLIP,
                       lw=1.7, zorder=4))
ellipse(631, 315, 17, 10, '#DED2EC', CLIP)
ellipse(689, 315, 14, 9, '#DED2EC', CLIP)
leader([(625, 352), (630, 320)], CLIP)
text(738, 311, 'Nascent RNA', anchor='left', color=CLIP)


rect(772, Y + 22, 128, 6, '#DFCEAF', DAM, 4)
for x in [795, 824, 853, 880]:
    line(x, Y + 28, x, Y + 36, DAM, .8)
    ellipse(x, Y + 38, 2.3, 2.3, DAM)
leader([(860, 352), (860, 303), (852, Y + 25)], DAM)


rect(310, Y - 6, 13, 12, PET, zorder=5)
rect(826, Y - 6, 13, 12, PET, zorder=5)
ax.add_patch(PathPatch(MPath([(316.5, Y - 6), (390, 133), (742, 133),
                              (832.5, Y - 6)],
                             [MPath.MOVETO] + [MPath.CURVE4] * 3),
                       facecolor='none', edgecolor=PET, lw=1.5, zorder=2))
text(575, 127, f"ChIA-PET ({groups['ChIAPET']:,})", size=18, bold=True, color=PET)
text(575, 107, 'CTCF / RAD21 / POLR2A / SMC1A')
text(575, 89, 'Paired interaction anchors', color='#5D6770')


fig.canvas.draw(); renderer = fig.canvas.get_renderer()
pixel_boxes = [(obj.get_text(), obj.get_window_extent(renderer)) for obj in texts]
for obj, (_, bb) in zip(texts, pixel_boxes):
    assert bb.x0 >= 0 and bb.y0 >= 0 and bb.x1 <= fig.bbox.width and bb.y1 <= fig.bbox.height, obj.get_text()
    assert obj.get_fontsize() == (18 if obj.get_text().startswith(('ATAC-seq (','ChIP-seq (','eCLIP (','DamID (','ChIA-PET (')) else 16)
    assert obj.get_fontfamily() == ['Arial']
    assert obj.get_color() == TEXT_COLOR
for i, (s, bb) in enumerate(pixel_boxes):
    for t, other in pixel_boxes[i + 1:]:
        assert min(bb.x1, other.x1) - max(bb.x0, other.x0) <= .5 or min(bb.y1, other.y1) - max(bb.y0, other.y0) <= .5, (s, t)
for a, b in leader_segments:
    points = ax.transData.transform(np.linspace(a, b, 150))
    for s, bb in pixel_boxes:
        collision = ((points[:, 0] > bb.x0) & (points[:, 0] < bb.x1)
                     & (points[:, 1] > bb.y0) & (points[:, 1] < bb.y1)).any()
        assert not collision, ('leader crosses text', s, a, b)
assert np.allclose(fig.get_size_inches() * 72, [1000, 520])
assert style.matplotlib.rcParams['pdf.fonttype'] == 42
assert sum(groups[k] for k in ['ChIP', 'eCLIP', 'ATAC', 'DamID', 'ChIAPET']) == 2554
try:
    style.save_pdf(fig, OUT)
except PermissionError:

    OUT = ROOT / 'fig3a_large_text.pdf'
    style.save_pdf(fig, OUT)
plt.close(fig)
audit = {'figure': OUT.name, 'canvas_pt': [1000, 520],
         'canvas_inches': [W / 72, H / 72],
         'canvas_size_exception': 'User explicitly requested no size change.',
         'skill': 'final-figure-small', 'font': 'Arial', 'ordinary_text_pt': 16,
         'assay_titles_pt': 18, 'font_size_exception': 'User explicitly requested larger lettering.',
         'text_color': TEXT_COLOR, 'all_text_black': True,
         'top_region_brackets_removed': True, 'footer_removed': True,
         'display_scope': 'Only the core gene and five assay callouts shown in user screenshot.',
         'pdf_fonttype': 42, 'single_shared_gene': True,
         'raster_artists': 0, 'text_within_canvas': True, 'no_text_overlap': True,
         'no_leader_text_intersection': True, 'regulatory_feature_count': 2554,
         'group_counts': {k: groups[k] for k in ['ChIP', 'eCLIP', 'ATAC', 'DamID', 'ChIAPET']},
         'verified_source_manifest_sha256': META['sha256'],
         'input_metadata_sha256': hashlib.sha256((ROOT/'verified_features.json').read_bytes()).hexdigest(),
         'pdf_sha256': hashlib.sha256(OUT.read_bytes()).hexdigest(),
         'scientific_scope': 'Illustrative positions; no measured locus or exclusive binding sites implied.'}
(ROOT / 'fig3a_audit.json').write_text(json.dumps(audit, indent=2), encoding='utf-8')
print('PASS: screenshot core only; 1000 x 520 pt, black Arial 16/18 pt, editable Type42 PDF; all text/leader checks passed.')
