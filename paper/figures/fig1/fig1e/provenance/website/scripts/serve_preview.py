"""Serve only the compiled frontend; source and server credentials are never exposed."""
import argparse
import functools
import json
from urllib.parse import urlsplit, parse_qs
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from gene_pi_tsne_jobs import GenePiTsneJobs
from genome_tracks import GenomeTracks

class PreviewHandler(SimpleHTTPRequestHandler):
    protocol_version = 'HTTP/1.1'

    def do_HEAD(self):
        url = urlsplit(self.path)
        if url.path.startswith('/api/v1/genome/reference/'):
            return self.server.genome.reference(self, url.path.rsplit('/', 1)[-1], head=True)
        return super().do_HEAD()

    def genome_api(self, url):
        if url.path.startswith('/api/v1/genome/reference/'):
            return self.server.genome.reference(self, url.path.rsplit('/', 1)[-1])
        try:
            status, value = self.server.genome.get(url.path.removeprefix('/api/v1/genome/'), parse_qs(url.query))
        except (ValueError, TypeError):
            status, value = 400, {'error': 'Invalid genome request. Use a valid region of at most 1 Mb.'}
        except Exception:
            status, value = 503, {'error': 'Genome data is temporarily unavailable. Please retry.'}
        body = json.dumps(value, separators=(',', ':'), allow_nan=False).encode()
        self.send_response(status)
        self.send_header('Content-Type', 'application/json')
        self.send_header('Content-Length', str(len(body)))
        self.send_header('Cache-Control', 'no-store')
        self.end_headers()
        try: self.wfile.write(body)
        except (BrokenPipeError, ConnectionResetError): pass

    def do_GET(self):
        url = urlsplit(self.path)
        if url.path.startswith('/api/v1/genome/'):
            return self.genome_api(url)
        if url.path != '/api/v1/analysis/gene-pi-tsne':
            return super().do_GET()
        params = parse_qs(url.query)
        try:
            symbol, transcript = params.get('symbol', [''])[0], params.get('transcript', [''])[0]
            if len(symbol) > 100 or len(transcript) > 100:
                raise ValueError('Invalid gene / transcript identifier.')
            status, value = self.server.tsne.get(symbol, transcript)
        except ValueError as error:
            status, value = 400, {'status': 'error', 'error': str(error)}
        except Exception:
            status, value = 503, {'status': 'error', 'error': 'PI t-SNE data is unavailable. Please retry.'}
        body = json.dumps(value, separators=(',', ':'), allow_nan=False).encode()
        self.send_response(status)
        self.send_header('Content-Type', 'application/json')
        self.send_header('Content-Length', str(len(body)))
        self.send_header('Cache-Control', 'no-store')
        self.end_headers()
        self.wfile.write(body)

class PreviewServer(ThreadingHTTPServer):
    request_queue_size = 128
    daemon_threads = True

if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--bind',default='127.0.0.1')
    parser.add_argument('--port',default=8088,type=int)
    parser.add_argument('--directory',default=str(Path(__file__).resolve().parent.parent/'dist'))
    args=parser.parse_args()
    root=Path(args.directory).resolve()
    if not (root/'index.html').is_file():
        parser.error('Build output index.html is missing.')
    handler=functools.partial(PreviewHandler,directory=str(root))
    print(f'Serving NascentAtlas dist at http://{args.bind}:{args.port}',flush=True)
    server = PreviewServer((args.bind,args.port),handler)
    server.tsne = GenePiTsneJobs(Path(__file__).resolve().parent.parent, root / 'data')
    server.genome = GenomeTracks(Path(__file__).resolve().parent.parent)
    server.serve_forever()
