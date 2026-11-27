"""Reproduce C-F from packaged frozen inputs; no SSH or original data dependency."""
from pathlib import Path
import json, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def build_figure(here, panel):
    data = here/'data'
    if panel == 'C':
        import plot_test_custom_palette as custom
        import plot_test_pvalues_only as approved
        records = pd.read_csv(data/'fold_metrics.tsv', sep='\t')
        original_make = custom.recolored_make_figure
        capture = {}
        def capture_make(frame):
            fig, stats = original_make(frame)
            capture.update(fig=fig, stats=stats)
            return fig, stats
        approved.make_figure = capture_make
        approved.retrieve = lambda unused: records.copy()
        old_argv = sys.argv
        try:
            sys.argv = ['plot.py', '--layout-only']
            approved.main()
        finally:
            sys.argv = old_argv
        fig, stats = capture['fig'], capture['stats']
        na = [t for t in fig.axes[0].texts if t.get_text() == '1']
        assert len(na) == 2
        for text in na:
            text.set_text('N/A')
        expected = json.loads((data/'reference_statistics.json').read_text())
        for actual, reference in zip(stats, expected):
            assert actual['base'] == reference['base'] and actual['candidate'] == reference['candidate']
            assert np.isclose(actual['adjusted_p'], reference['adjusted_p'], rtol=1e-6, atol=1e-15)
        result = {'comparisons':stats, 'test_folds':10,
                  'NA_display_only':'Two adjusted P=1 values are displayed as N/A for display; they are not missing statistics.'}
    elif panel == 'D':
        import redraw_figures as original
        from scipy.stats import pearsonr, spearmanr
        frame = pd.read_csv(data/'oof_predictions.tsv', sep='\t')
        assert len(frame) == 15351 and frame.gene.is_unique
        assert sorted(frame.fold.unique()) == list(range(1,11))
        packet = {'oof':{'observed':frame.Y_PC1.tolist(),
                         'fixed_prediction':frame.test_prediction.tolist()}}
        captured = {}
        def capture_save(fig, path, title):
            original.check_layout(fig, 'fig4.d.pdf')
            captured['fig'] = fig
        original.save_one = capture_save
        original.draw_oof(packet, here, local=False)
        fig = captured['fig']
        result = {'n':len(frame), 'pearson_r':float(pearsonr(frame.Y_PC1, frame.test_prediction).statistic),
                  'spearman_rho':float(spearmanr(frame.Y_PC1, frame.test_prediction).statistic),
                  'scope':'Fixed BQIP, five-seed ensemble, 10 held-out OOF folds; feature-selection evaluation exploratory.'}
        assert round(result['pearson_r'],4) == .7276 and round(result['spearman_rho'],4) == .7495
    elif panel == 'E':
        import plot_rf_importance_reference as original
        packet = json.loads((data/'importance_packet.json').read_text())
        fig, qa = original.build_figure(packet)
        result = {'qa':qa, 'top30':packet['top30'], 'provenance':packet['provenance']}
    elif panel == 'F':
        import heatmap_legend_left_colorbar_right as original
        import redraw_current_correlation_heatmap as base
        base.LEGACY_CODE = here/'code/legacy_heatmap.py'
        packet = json.loads((data/'correlation_packet.json').read_text())
        fig = original.make_figure(packet)
        result = {'matrix_shape':[len(packet['rows']),len(packet['columns'])],
                  'provenance':packet['provenance'], 'clustered_rows':packet['clustered_rows'],
                  'clustered_columns':packet['clustered_columns'],
                  'layout':'Bottom feature-group legend LEFT, Pearson correlation colorbar RIGHT'}
    else:
        raise ValueError(panel)


    from matplotlib.backends.backend_agg import FigureCanvasAgg
    FigureCanvasAgg(fig)
    assert matplotlib.rcParams['pdf.fonttype'] == 42
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    from matplotlib.text import Text
    for text in fig.findobj(Text):
        if not text.get_visible() or not text.get_text().strip():
            continue
        b = text.get_window_extent(renderer)
        assert b.x0 >= -1 and b.x1 <= fig.bbox.x1+1, text.get_text()
        assert b.y0 >= -1 and b.y1 <= fig.bbox.y1+1, text.get_text()
    return fig, result
