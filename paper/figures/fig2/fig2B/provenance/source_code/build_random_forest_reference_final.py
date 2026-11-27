
from pathlib import Path

root=Path(__file__).resolve().parent
source=(root/'build_random_forest_reference_layout.py').read_text(encoding='utf-8')
for before,after in {
    '#929292':'#555555',
    '#B6B6B6':'#777777',
    '#393939':'#282828',
    '#8D3395':'#6D1878',
    '#40A6C5':'#087FA1',
    '#D9E1E8':'#ADBFD0',
    '#9BAEC8':'#607E9F',
    'random_forest_reference_layout.pdf':'random_forest_reference_layout_final.pdf',
    "    texts.append(ax.text(x,y,s,ha='center'":
        "    fs=14 if s in ('Features (786)','DNA / annotations','Predicted Y_PC1') else 12\n    texts.append(ax.text(x,y,s,ha='center'",
    "text(244,414,'validation only'":"text(244,430,'validation only'",
    "text(730,241,'B + Q + I + C + P + E + T'":"text(730,241,'B Q I C P E T'",
    "text(730,274,'DNA k-mers: 1-4 bp'":"text(730,284,'DNA k-mers: 1-4 bp'",
    "text(760,421,'Random forest: 300 trees'":"text(730,421,'Random forest: 300 trees'",
}.items():
    assert before in source,before
    source=source.replace(before,after)

start=source.index('# Regression outputs are averaged')
end=source.index('fig.canvas.draw();renderer=',start)
source=source[:start]+'''# Regression mean, with a single uncluttered aggregation junction.
for cx in [175,500,825]:
    line(cx,671,cx,690,GREY,1.1)
line(175,690,825,690,GREY,1.1)
line(500,690,500,695,GREY,1.1)
text(500,714,'Mean',12,col=GREY)
arrow(500,737,500,747)
text(500,778,'Predicted Y_PC1',14,True,col=PURPLE)

'''+source[end:]


source=source.replace('save_pdf(fig,OUT)', '''overlaps=[]
for i,t in enumerate(texts):
    a=t.get_window_extent(renderer)
    for u in texts[i+1:]:
        if a.overlaps(u.get_window_extent(renderer)):
            overlaps.append((t.get_text(),u.get_text()))
assert not overlaps,overlaps
save_pdf(fig,OUT)''')
exec(compile(source,str(__file__),'exec'),{'__file__':str(__file__),'__name__':'__main__'})
