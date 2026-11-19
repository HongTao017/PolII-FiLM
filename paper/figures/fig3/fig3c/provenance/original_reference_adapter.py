"""Execute final referenced-chat plot functions on real REG5 results.

Original files are frozen byte-for-byte copies from the referenced task. The
scatter function is called unmodified. Violin functions receive only necessary
data-schema adaptations: nine combinations instead of twelve, six feature-key
entries instead of seven, test r range through .84, and original REG5 Bonferroni
family of14. Original final colours, artist styles, numeric brackets and font
sizes are retained. Explicit current user override: both canvases9.4x5.3 inches.
Only two final PDFs are written to the specified output directory.
"""
import argparse
import csv
import hashlib
import inspect
import io
import json
from pathlib import Path
import sys

sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
ORIGINAL=HERE/'original_code'
sys.path.insert(0,str(ORIGINAL))
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from matplotlib.text import Text
import numpy as np
import pandas as pd
from scipy.stats import pearsonr,spearmanr
from final_figure_small import new_figure,save_pdf
import plot_test_compact_legend as reference_core
import plot_test_muted_filled_reference as reference_fill
import plot_test_custom_palette as reference_custom
import plot_test_pvalues_only as reference_numeric
import redraw_figures as reference_scatter

ROOT=HERE.parent
CANVAS=(9.4,5.3)
SOURCE_THREAD='01a13f40-5270-7590-9efd-4b4dad994c44'
PINNED={
 'redraw_figures.py':'2af4d3f9f1f2e4f4767080f4d2b92bedd4936e1ecea51f7723bf4c57508db675',
 'plot_test_compact_legend.py':'fa327a87aa1d34cfdd057af72fdced40562fe78cbd1a97a0437dfbc8fb9e6e2c',
 'plot_test_muted_filled_reference.py':'54d5f63dc3b5056505554745e3ddc56686f3cfb7ba0c3016afc7591e9e1a2fc9',
 'plot_test_custom_palette.py':'7512e77edd083ade95ca5a755ef41c9ae90933b1e8800e894598efe2a6dfa958',
 'plot_test_pvalues_only.py':'26ad736f5fa4b0e45253325f64abf886a21376b4fae2425a60c4125e1377ffad'}
DISPLAY=['B','BL','BP','BD','BA','BC','BCL','BCP','BCLP']
PALETTE=['#FEFCFA','#B7BEC9','#4D5058','#F5DCD1','#F0BDBD','#928C95','#B2CEF7']
COLORBAR=['#1e367b','#703b90','#c93659','#f37c22','#f7df54']
LETTERS={'B':'B','ChIP':'C','eCLIP':'L','ChIAPET':'P','DamID':'D','ATAC':'A'}
KEY=[('B','Basic gene\nfeatures (22)',.25,.175),('C','ChIP-seq\n(2,284)',.55,.175),
     ('L','eCLIP\n(206)',.85,.175),('P','ChIA-PET\n(50)',.25,.075),
     ('D','DamID\n(10)',.55,.075),('A','ATAC-seq\n(4)',.85,.075)]
TIERS={('B','BC'):4,('B','BL'):0,('B','BP'):1,('B','BD'):2,('B','BA'):3,
       ('BC','BCP'):1,('BCL','BCLP'):4}
ADAPTATIONS=[]

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def adapted(function,namespace,replacements,name=None):
    source=inspect.getsource(function)
    for before,after in replacements:
        assert source.count(before)==1,('source mismatch',function.__name__,before)
        source=source.replace(before,after)
        ADAPTATIONS.append({'function':function.__name__,'before':before,'after':after})
    scope=dict(namespace);scope['print']=lambda *a,**k:None
    exec(compile(source,'<reference-data-schema-adapter>','exec'),scope)
    return scope[name or function.__name__]

def check_text(fig):
    fig.canvas.draw();renderer=fig.canvas.get_renderer()
    texts=[t for t in fig.findobj(Text) if t.get_visible() and t.get_text().strip()]
    for text in texts:
        b=text.get_window_extent(renderer)
        assert b.x0>=-1 and b.y0>=-1 and b.x1<=fig.bbox.width+1 and b.y1<=fig.bbox.height+1,('clipping',text.get_text())
    for i,text in enumerate(texts):
        b=text.get_window_extent(renderer)
        for other in texts[i+1:]:
            c=other.get_window_extent(renderer)
            assert min(b.x1,c.x1)-max(b.x0,c.x0)<=1 or min(b.y1,c.y1)-max(b.y0,c.y0)<=1,('text overlap',text.get_text(),other.get_text())
    assert np.allclose(fig.get_size_inches(),CANVAS)
    assert matplotlib.rcParams['pdf.fonttype']==42

def build_violin(packet,all_comparisons):
    pairs=[(r['base_code'],r['candidate_code'],TIERS[(r['base_code'],r['candidate_code'])])
           for r in all_comparisons if float(r['adjusted_p'])<.05]
    assert len(all_comparisons)==14 and len(pairs)==7
    rows=[]
    for model in packet['models']:
        code=''.join(LETTERS[g] for g in packet['groups'][model])
        for fold,value in enumerate(packet['scores'][model]['test'],1):
            rows.append({'model_id':code,'scope':'test','fold':fold,'pearson_r':value})
    data=pd.DataFrame(rows)
    common={'MODELS':DISPLAY,'CANVAS':CANVAS,'COMPARISONS':pairs,
            'FEATURE_KEY':KEY,'LETTER_COLORS':dict.fromkeys('BCLPDA','black')}
    ns=dict(vars(reference_core));ns.update(common)
    core=adapted(reference_core.make_figure,ns,[
        ('positions=np.arange(12)','positions=np.arange(len(MODELS))'),
        ("assert len(boxes['boxes'])==12","assert len(boxes['boxes'])==len(MODELS)"),
        ('assert y.min()>.54 and y.max()<.78','assert y.min()>.54 and y.max()<.84'),
        ('adjusted=min(1.,raw_p*len(COMPARISONS))','adjusted=min(1.,raw_p*14)'),
        ('ax.set_xlim(-.4,11.4);ax.set_ylim(.54,.78)','ax.set_xlim(-.4,len(MODELS)-.6);ax.set_ylim(.54,.84)'),
        ("ax.set_xticks(positions);ax.set_xticklabels(['']*12)","ax.set_xticks(positions);ax.set_xticklabels(['']*len(MODELS))"),
        ('ax.set_yticks([.55,.60,.65,.70,.75])','ax.set_yticks([.55,.60,.65,.70,.75,.80])')])
    ns=dict(vars(reference_fill));ns.update(common)
    # Intermediate legacy assertions use fixed old counts; final QA below checks
    # the actual new counts, fonts, exact colours, numeric tests and all bounds.
    ns.update(original_make_figure=core,validate_layout=lambda *args:None)
    fill=adapted(reference_fill.make_figure,ns,[
        ('assert len(ax.collections) == 12 and len(ax.patches) == 12','assert len(ax.collections) == len(MODELS) and len(ax.patches) == len(MODELS)'),
        ('assert len(fig.texts) == 7 and len(fig.artists) == 1','assert len(fig.texts) == 6 and len(fig.artists) == 1'),
        ('[([0, 1, 2, 3], .127), ([4, 5, 6], .078)]','[([0, 1, 2], .127), ([3, 4, 5], .078)]')])
    ns=dict(vars(reference_custom));ns.update(common);ns.update(checked_make_figure=fill)
    custom=adapted(reference_custom.recolored_make_figure,ns,[
        ('assert len(ax.collections) == len(MODELS) == 12','assert len(ax.collections) == len(MODELS)')])
    source=inspect.getsource(reference_numeric.main)
    start=source.index('    fig,stats=make_figure(retrieve(False));ax=fig.axes[0]')
    end=source.index('    if not args.layout_only:')
    source='def build_numeric(data):\n'+source[start:end]+'    return fig,stats\n'
    replacements=[('make_figure(retrieve(False))','make_figure(data)'),
        ('assert len(old_brackets)==11','assert len(old_brackets)==len(COMPARISONS)'),
        ('assert len(ax.collections)==12 and len(ax.patches)==12 and len(fig.texts)==7',
         'assert len(ax.collections)==len(MODELS) and len(ax.patches)==len(MODELS) and len(fig.texts)==6'),
        ('assert len(ax.texts)==11','assert len(ax.texts)==len(COMPARISONS)')]
    for before,after in replacements:
        assert source.count(before)==1,before
        source=source.replace(before,after)
        ADAPTATIONS.append({'function':'numeric main (in-memory drawing section)','before':before,'after':after})
    ns=dict(vars(reference_numeric));ns.update(common)
    ns.update(make_figure=custom,print=lambda *args,**kwargs:None)
    exec(compile(source,'<reference-numeric-drawing-section>','exec'),ns)
    fig,stats=ns['build_numeric'](data)
    ax=fig.axes[0]
    assert [t.get_text() for t in ax.get_xticklabels()]==DISPLAY
    assert len(ax.collections)==len(ax.patches)==9 and len(ax.texts)==7
    assert all(np.allclose(body.get_facecolor()[0],to_rgba(PALETTE[i%7])) for i,body in enumerate(ax.collections))
    assert all(t.get_fontsize()==10 for t in ax.get_xticklabels())
    assert all(t.get_fontsize()==12 for t in ax.get_yticklabels())
    assert all(t.get_fontsize()==8.5 for t in fig.texts)
    assert all(t.get_fontsize()==8.5 for t in ax.texts)
    assert ax.yaxis.label.get_fontsize()==14
    assert len(fig.texts)==6 and len(set(t.get_position()[1] for t in fig.texts))==2
    assert not ax.spines['top'].get_visible() and not ax.spines['right'].get_visible()
    expected=[r for r in all_comparisons if float(r['adjusted_p'])<.05]
    for actual,wanted in zip(stats,expected):
        assert (actual['base'],actual['candidate'])==(wanted['base_code'],wanted['candidate_code'])
        assert abs(actual['adjusted_p']-float(wanted['adjusted_p']))<1e-14
    check_text(fig)
    return fig,{'displayed_codes':DISPLAY,'palette':PALETTE,'comparisons':expected,
        'comparison_family':14,'original_comparisons':all_comparisons,
        'filter':'Keep endpoints of at least one corrected two-sided p<.05 contrast; not new model selection',
        'fonts':{'x_ticks':10,'y_ticks':12,'axis_label':14,'key_and_p':8.5},
        'data_fingerprint':packet['fingerprint']}

def build_scatter(packet):
    capture={}
    def inherited_figure(**kwargs):
        kwargs.setdefault('figsize',CANVAS)
        return new_figure(**kwargs)
    def capture_save(fig,path,title):
        for text in fig.texts:
            if text.get_text()=='Fixed BQIP: exploratory OOF':
                text.set_text('Fixed B + ChIP-seq + eCLIP: exploratory OOF')
        capture['fig']=fig
    reference_scatter.new_figure=inherited_figure
    reference_scatter.CANVAS=CANVAS
    reference_scatter.save_one=capture_save
    source={'oof':{'observed':packet['observed'],'fixed_prediction':packet['predicted']}}
    # Exact final referenced-chat function, no density/colour/style rewriting.
    reference_scatter.draw_oof(source,HERE,local=False)
    fig=capture['fig'];ax=fig.axes[0];image=ax.collections[0]
    assert int(np.sum(image.get_array()))==15351
    assert image.get_rasterized() is True
    from matplotlib.colors import LinearSegmentedColormap
    wanted=LinearSegmentedColormap.from_list('expected_reference',COLORBAR)
    assert np.allclose(image.cmap(np.linspace(0,1,257)),wanted(np.linspace(0,1,257)))
    r=float(pearsonr(packet['observed'],packet['predicted']).statistic)
    rho=float(spearmanr(packet['observed'],packet['predicted']).statistic)
    assert abs(r-packet['metrics']['pearson_r'])<1e-12
    assert abs(rho-packet['metrics']['spearman_rho'])<1e-12
    assert any(f'r = {r:.4f}' in t.get_text() and f'rho = {rho:.4f}' in t.get_text() for t in ax.texts)
    reference_scatter.check_layout(fig,'regulatory exact-reference scatter')
    check_text(fig)
    return fig,{'model':packet['model'],'features':packet['feature_count'],'metrics':packet['metrics'],
        'source_oof_sha256':packet['oof_sha256'],'data_fingerprint':packet['fingerprint'],
        'scope':packet['scope'],'method':'Original draw_oof: hexbin, gridsize65, mincnt1, bins=log, no count/point substitution',
        'palette':COLORBAR,'hexagon_gene_count_sum':int(np.sum(image.get_array())),
        'hexagons_rasterized_as_original':True,'all_text_editable':True}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--directory',type=Path,required=True)
    parser.add_argument('--check-only',action='store_true');args=parser.parse_args()
    for filename,digest in PINNED.items():assert sha(ORIGINAL/filename)==digest,filename
    all_source_hashes={f.name:sha(f) for f in ORIGINAL.glob('*.py')}
    violin_packet=json.loads((ROOT/'reference_violin_packet.json').read_text(encoding='utf-8-sig'))
    scatter_packet=json.loads((ROOT/'fig3b_final_style/selected_scatter_packet.json').read_text(encoding='utf-8-sig'))
    assert violin_packet['fingerprint']==scatter_packet['fingerprint']=='c13b03947d1c64bd21ad4525860c8a1528cb98ff6d9d711936451f360cbaf0a4'
    with (ROOT/'reference_violin_significant_groups_preview/test_paired_comparisons.tsv').open(encoding='utf-8') as f:
        comparisons=list(csv.DictReader(f,delimiter='\t'))
    results=[('fig3b_violin_REG5_exact_chat_code_9p4x5p3.pdf',build_violin(violin_packet,comparisons)),
             ('fig3b_scatter_REG5_exact_chat_code_9p4x5p3.pdf',build_scatter(scatter_packet))]
    outputs=[]
    assert args.directory.is_dir()
    for filename,(fig,audit) in results:
        target=args.directory/filename
        if not args.check_only and target.exists():raise FileExistsError(target)
        metadata={'Title':filename,'Subject':json.dumps({'reference_thread':SOURCE_THREAD,
            'original_source_sha256':all_source_hashes,'canvas_inches':CANVAS,
            'user_override':'Preserve current9.4x5.3 canvas; retained nine-group/significant-only display',
            'data':audit,'data_schema_adaptations':ADAPTATIONS},ensure_ascii=True),
            'Keywords':'Original referenced functions; REG5 actual data; PDF only; no original files changed'}
        raw=fig.savefig
        def with_metadata(path,_save=raw,_meta=metadata,**kwargs):
            return _save(path,metadata=_meta,**kwargs)
        fig.savefig=with_metadata
        if args.check_only:
            stream=io.BytesIO();fig.savefig(stream,format='pdf')
            outputs.append({'filename':filename,'memory_pdf_bytes':len(stream.getvalue()),'audit':audit})
        else:
            save_pdf(fig,target)
            outputs.append({'path':str(target),'bytes':target.stat().st_size,'sha256':sha(target),'audit':audit})
        plt.close(fig)
    assert all_source_hashes=={f.name:sha(f) for f in ORIGINAL.glob('*.py')}
    print(json.dumps({'status':'PASS','check_only':args.check_only,'reference_source_hashes_verified':True,
        'original_sources_unchanged':True,'canvas_inches':CANVAS,'outputs':outputs},ensure_ascii=True))

if __name__=='__main__':main()
