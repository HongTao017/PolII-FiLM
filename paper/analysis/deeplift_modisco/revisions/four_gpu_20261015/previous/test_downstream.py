"""Tool integration smoke tests; all synthetic outputs isolated under tests/."""
import inspect
import json
from pathlib import Path
import subprocess

import h5py
import numpy as np
import pandas as pd
from modiscolite.tfmodisco import TFMoDISco
from modiscolite.io import save_hdf5

from run_pipeline import MEME_HEADER, meme_entry, scan_fimo
from deep_lift_core import write_json

root=Path(__file__).resolve().parent
config=json.loads((root/'experiment_config.json').read_text())
test=root/'tests'
test.mkdir(exist_ok=True)
parameters=inspect.signature(TFMoDISco).parameters
assert all(name in parameters for name in ('one_hot','hypothetical_contribs','max_seqlets_per_metacluster','final_min_cluster_size'))
save_hdf5(str(test/'empty_modisco_format.h5'),[],[],100)
with h5py.File(test/'empty_modisco_format.h5','r') as handle:
    assert 'pos_patterns' in handle and 'neg_patterns' in handle
motifs=json.loads((Path(config['selection_root'])/'motifs.json').read_text())
record=motifs['EGR1']
matrix=np.array([record['pfm'][c] for c in 'ACGT']).T
consensus=''.join('ACGT'[int(x.argmax())] for x in matrix)
(test/'motif.meme').write_text(MEME_HEADER+meme_entry(record['matrix_id'],matrix,100))
(test/'synthetic_consensus.fa').write_text(f'>synthetic_consensus_NOT_EXPERIMENTAL\n{consensus}\n')
scan_fimo(config,test/'motif.meme',test/'synthetic_consensus.fa',test/'synthetic_fimo.tsv')
fimo=pd.read_csv(test/'synthetic_fimo.tsv',sep='\t',comment='#')
assert len(fimo)>0 and fimo['p-value'].min()<=config['fimo_p_threshold']
with (test/'synthetic_tomtom.log').open('w') as logfile:
    subprocess.run([config['tomtom'],'-oc',str(test/'synthetic_tomtom'),'-min-overlap','5',
        '-thresh','.05',str(test/'motif.meme'),str(test/'motif.meme')],stdout=logfile,stderr=subprocess.STDOUT,check=True)
tomtom=pd.read_csv(test/'synthetic_tomtom'/'tomtom.tsv',sep='\t',comment='#')
assert len(tomtom)>0 and tomtom['q-value'].min()<=.05
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import logomaker
fig,ax=plt.subplots(figsize=(6,2))
logomaker.Logo(pd.DataFrame(matrix/matrix.sum(1,keepdims=True),columns=list('ACGT')),ax=ax)
ax.set_title('SYNTHETIC TOOL TEST — NOT AN EXPERIMENT RESULT')
fig.savefig(test/'synthetic_plot_smoke.png',dpi=100)
plt.close(fig)
write_json(root/'downstream_tests_complete.json',dict(status='PASS',TFMoDISco_signature=True,HDF5_format=True,
    FIMO_consensus_p=float(fimo['p-value'].min()),TOMTOM_identical_matrix_q=float(tomtom['q-value'].min()),
    plot_smoke=True,synthetic_only=True))
print('DOWNSTREAM TESTS PASS')
