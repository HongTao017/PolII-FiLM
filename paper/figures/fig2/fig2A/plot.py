"""Offline entry point for a single archived Figure 2 panel."""
from pathlib import Path
import sys,argparse,json,hashlib
ROOT=Path(__file__).resolve().parent
CFG=json.loads((ROOT/'package.json').read_text(encoding='utf-8'))
parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=ROOT/CFG['filename'])
args=parser.parse_args();out=args.output
if not out.is_absolute():out=ROOT/out
out=out.resolve();out.parent.mkdir(parents=True,exist_ok=True)
assert out.suffix.lower()=='.pdf' and (out==ROOT.resolve() or ROOT.resolve() in out.parents)
import subprocess
_original_subprocess_run=subprocess.run
def _offline_run(command,*a,**kw):
    name=str(command[0] if isinstance(command,(list,tuple)) else command).split()[0].lower()
    if Path(name).name in {'ssh','ssh.exe','scp','scp.exe','curl','wget'}:
        raise RuntimeError('Offline figure entry point forbids network acquisition')
    return _original_subprocess_run(command,*a,**kw)
subprocess.run=_offline_run
letter=CFG['panel'];sys.argv=[sys.argv[0]]

if letter in 'AB':
    from pypdf import PdfReader,PdfWriter,Transformation
    from pypdf.generic import RectangleObject
    cfg=json.loads((ROOT/'data/diagram_definition.json').read_text(encoding='utf-8'))
    r=PdfReader(ROOT/'data/source_fig2.ai')
    w=PdfWriter();page=w.add_page(r.pages[0]);page.transfer_rotation_to_content()
    left,top,right,bottom=cfg['crop_top_left_pt'];source_h=float(page.mediabox.height)
    low=source_h-bottom
    page.add_transformation(Transformation().translate(-left,-low))
    page.mediabox=RectangleObject([0,0,right-left,bottom-top])
    page.cropbox=RectangleObject([0,0,right-left,bottom-top])
    w.add_metadata({'/Title':'Figure 2'+letter+'; archived vector panel','/Subject':'Vector extraction from original Illustrator composition; manual source edits preserved.'})
    with out.open('wb') as f:w.write(f)
else:
    sys.path.insert(0,str(ROOT/'code'))
    import numpy as np,pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    for f in (ROOT/'data/fonts').glob('*.ttf'):font_manager.fontManager.addfont(str(f))
    if letter=='C':
        import plot_test_tf_binding_site as original
        final=original.original.custom.approved
        data=pd.read_csv(ROOT/'data/fold_metrics.tsv',sep='\t')
        final.retrieve=lambda *_:data.copy()
        class PanelPages(original.original.NALabelPdfPages):
            def savefig(self,figure=None,**kwargs):
                figure.axes[0].set_ylabel('Pearson r',fontsize=14)
                return super().savefig(figure,**kwargs)
        original.original.custom._pdf_pages=PanelPages
        final.make_figure=original.make_figure
        final.DESTINATION=out
        final.PdfPages=original.pdf_pages
        final.main()
    elif letter=='D':
        import redraw_figures as original
        data=pd.read_csv(ROOT/'data/oof_predictions.tsv',sep='\t')
        packet={'oof':{'observed':data.Y_PC1.tolist(),'fixed_prediction':data.test_prediction.tolist()}}
        from scipy.stats import pearsonr,spearmanr
        metrics={'n':len(data),'pearson_r':float(pearsonr(data.Y_PC1,data.test_prediction).statistic),'spearman_rho':float(spearmanr(data.Y_PC1,data.test_prediction).statistic)}
        assert len(data)==15351 and round(metrics['pearson_r'],4)==.7276 and round(metrics['spearman_rho'],4)==.7495
        def save(fig,*_):
            for t in list(fig.texts):t.remove()
            fig.axes[0].set_title('BQIP',fontsize=14,pad=12)
            for t in fig.axes[0].texts:t.set_text(t.get_text().split('\nn =')[0])
            fig.axes[1].set_ylabel('')
            fig.savefig(out)
            plt.close(fig)
        original.save_one=save
        original.draw_oof(packet,ROOT,local=False)
        (ROOT/'data/recomputed_metrics.json').write_text(json.dumps(metrics,indent=2),encoding='utf-8')
    elif letter=='E':
        import plot_rf_importance_reference as original
        packet=json.loads((ROOT/'data/plot_packet.json').read_text())
        fig,qa=original.build_figure(packet)
        for t in list(fig.texts):t.remove()
        fig.savefig(out);plt.close(fig)
    elif letter=='F':
        import redraw_current_correlation_heatmap as original
        original.LEGACY_CODE=ROOT/'code/legacy_heatmap.py'
        packet=json.loads((ROOT/'data/plot_packet.json').read_text())
        fig=original.make_figure(packet)
        if fig._suptitle is not None:fig._suptitle.remove()
        fig.savefig(out);plt.close(fig)
    assert matplotlib.rcParams['pdf.fonttype']==42

record={'panel':letter,'output':str(out.relative_to(ROOT)),'sha256':hashlib.sha256(out.read_bytes()).hexdigest(),
        'offline_redraw':True,'data_refit_or_retraining':False,'output_format':'PDF only'}
(ROOT/'audit'/('redraw_'+out.stem+'.json')).write_text(json.dumps(record,indent=2),encoding='utf-8')
print('SAVED',str(out))
