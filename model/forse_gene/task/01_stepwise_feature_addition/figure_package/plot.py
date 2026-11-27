"""Portable, offline figure redraw. Only the explicitly selected PDF is exported."""
from pathlib import Path
import argparse, io, json, os, sys, hashlib
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
os.environ.setdefault('MPLBACKEND', 'Agg')
sys.path.insert(0, str(HERE/'code'))
from matplotlib import font_manager
for font in (HERE/'data/fonts').glob('*.ttf'):
    font_manager.fontManager.addfont(str(font))
from render import build_figure
from final_figure_small import save_pdf


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    parser.add_argument('--verify-only', action='store_true', help='Render PDF in memory; do not overwrite the final PDF')
    args = parser.parse_args()
    meta = json.loads((HERE/'package.json').read_text())
    for name, digest in meta['frozen_input_sha256'].items():
        assert hashlib.sha256((HERE/name).read_bytes()).hexdigest() == digest, name
    fig, stats = build_figure(HERE, meta['panel'])
    fig.canvas.draw()
    if args.verify_only:
        stream = io.BytesIO()
        fig.savefig(stream, format='pdf')
        print(json.dumps({'panel':meta['panel'], 'status':'PASS', 'memory_pdf_bytes':len(stream.getvalue()), 'statistics':stats}))
    else:
        output = (args.output or HERE/meta['filename']).resolve()
        assert output.suffix.lower() == '.pdf', 'PDF only'
        output.parent.mkdir(parents=True, exist_ok=True)
        original_save = fig.savefig
        def with_metadata(path, **kwargs):
            return original_save(path, metadata={'Title':meta['description'], 'Subject':json.dumps(stats),
                'Keywords':'Portable frozen-data redraw; archived artwork retained in initial PDF'}, **kwargs)
        fig.savefig = with_metadata
        save_pdf(fig, output)
    import matplotlib.pyplot as plt
    plt.close(fig)


if __name__ == '__main__':
    main()
