
import json
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[3]
FOLD = ROOT / 'folds' / 'fold_02'
sys.path.insert(0, str(FOLD))
import run_pipeline


def main():
    config = json.loads((FOLD / 'experiment_config.json').read_text())
    scans = FOLD / 'instance_scans'
    no_hit = scans / 'r2_neg_patterns_pattern_6'
    with tempfile.TemporaryDirectory(prefix='fimo_repair_check_') as temporary:
        area = Path(temporary)
        output = area / 'empty.tsv'
        run_pipeline.scan_fimo(config, no_hit / 'target.meme', no_hit / 'seqlets.fa', output)
        lines = output.read_text().splitlines()
        assert len(lines) == 1 and 'sequence_name' in lines[0]
        hit = scans / 'r0_pos_patterns_pattern_0'
        output = area / 'hit.tsv'
        run_pipeline.scan_fimo(config, hit / 'target.meme', hit / 'seqlets.fa', output)
        assert output.read_bytes() == (hit / 'fimo.tsv').read_bytes()
        fake = area / 'false_success.py'
        fake.write_text('#!/usr/bin/env python3\nimport sys\nsys.stderr.write("Error: motif parse failure\\n")\n')
        fake.chmod(0o755)
        bad = dict(config, fimo=str(fake))
        try:
            run_pipeline.scan_fimo(bad, no_hit / 'target.meme', no_hit / 'seqlets.fa', area / 'bad.tsv')
        except AssertionError:
            pass
        else:
            raise AssertionError('Suspicious empty output was accepted')
        fake.write_text('#!/usr/bin/env python3\nimport sys\nsys.exit(3)\n')
        try:
            run_pipeline.scan_fimo(bad, no_hit / 'target.meme', no_hit / 'seqlets.fa', area / 'exit3.tsv')
        except subprocess.CalledProcessError:
            pass
        else:
            raise AssertionError('Nonzero FIMO exit was accepted')
    print(json.dumps(dict(status='PASS', verified=['zero_hit_header', 'positive_hit_identical',
        'suspicious_diagnostic_rejected', 'nonzero_exit_rejected'])))


if __name__ == '__main__':
    main()
