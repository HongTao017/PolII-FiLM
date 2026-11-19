"""Verify both native TOMTOM pairing and the selected gene's FIMO site."""
from pathlib import Path
import json,csv,subprocess,hashlib,xml.etree.ElementTree as ET
import numpy as np,h5py
ROOT=Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019')
R=Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motif/folds/fold_03')
S=ROOT/'analysis/all_motif_ISM_screen_20261020/exact_ism/candidate_003'
O=ROOT/'analysis/PTGIS_EGR1_TOMTOM_matched_20261020';O.mkdir(exist_ok=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,obj):p.write_text(json.dumps(obj,ensure_ascii=False,indent=2))
complete=json.loads((R/'complete.json').read_text());sources={}
for name in ['tomtom/tomtom.tsv','tomtom/tomtom.xml','modisco/tss_sequence.h5','patterns.meme','jaspar_125.meme']:
    p=R/name;sources[str(p)]=sha(p);assert sources[str(p)]==complete['outputs_sha256'][name]
plan=json.loads((S/'frozen_plan.json').read_text());assert (plan['gene'],plan['TF'],plan['matrix_id'])==('PTGIS','EGR1','MA0162.5')
inference=json.loads((S/'complete.json').read_text())
for name,expected in inference['outputs_sha256'].items():assert sha(S/name)==expected
rows=list(csv.DictReader((R/'tomtom/tomtom.tsv').open(),delimiter='\t'))
row=next(x for x in rows if x['Query_ID']==plan['pattern'] and x['Target_ID']==plan['matrix_id'])
assert float(row['q-value'])<.05
xml=ET.parse(R/'tomtom/tomtom.xml').getroot();queries=list(xml.find('queries'));targets=list(xml.find('targets'))
qi=next(i for i,n in enumerate(queries) if n.attrib['id']==plan['pattern']);ti=next(i for i,n in enumerate(targets) if n.attrib['id']==plan['matrix_id'])
query=queries[qi];target=targets[ti];native=next(n for q in xml.find('matches') if int(q.attrib['idx'])==qi for n in q if int(n.attrib['idx'])==ti)
assert int(native.attrib['off'])==int(row['Optimal_offset'])
assert (native.attrib['rc']=='y')==(row['Orientation']=='-')
def matrix(node):
    arr=np.array([[float(p.attrib[b]) for b in 'ACGT'] for p in node]);assert np.max(abs(arr.sum(1)-1))<3e-6
    return arr/arr.sum(1,keepdims=True)
qp=matrix(query);tp=matrix(target);oriented=tp[::-1,::-1] if row['Orientation']=='-' else tp
offset=int(row['Optimal_offset']);left=max(0,-offset);right=min(len(qp),len(oriented)-offset);tl=left+offset;tr=right+offset
assert right-left==int(row['Overlap'])==10 and len(tp)==10
cfg=json.loads((R/'experiment_config.json').read_text());M=json.loads((R/'motifs.json').read_text())['EGR1'];counts=np.array([M['pfm'][b] for b in 'ACGT']).T;ppm=(counts+1e-8)/(counts+1e-8).sum(1,keepdims=True)
# TOMTOM's native XML stores probabilities after its default0.1 motif pseudocount.
reference_lines=(R/'jaspar_125.meme').read_text().splitlines()
refidx=next(i for i,line in enumerate(reference_lines) if line.split()[:2]==['MOTIF','MA0162.5'])
raw_reference=np.array([[float(x) for x in line.split()] for line in reference_lines[refidx+2:refidx+12]])
assert np.allclose(ppm,raw_reference,atol=2e-8,rtol=1e-7),'FIMO and original TOMTOM input reference differ'
help_run=subprocess.run([cfg['tomtom'],'-h'],capture_output=True,text=True)
help_text=help_run.stdout+help_run.stderr
(O/'tomtom_help.txt').write_text(help_text)
assert 'motif-pseudo' in help_text and '0.1' in help_text
xml_expected=(raw_reference*float(target.attrib['nsites'])+.1*.25)/(float(target.attrib['nsites'])+.1)
assert np.allclose(xml_expected,tp,atol=1e-6,rtol=1e-5),'Native XML smoothing must be explained by recorded TOMTOM pseudocount'
print(json.dumps(dict(stage='reference_verified',raw_reference_equal=True,XML_default_pseudocount=.1,XML_max_rounding_error=float(abs(xml_expected-tp).max()))),flush=True)
with h5py.File(R/'modisco/tss_sequence.h5','r') as f:
    g=f['pos_patterns/pattern_0'];q=g['seqlets'];k=plan['seqlet'];onehot=q['sequence'][k]
    if plan['seqlet_is_revcomp']:onehot=onehot[::-1,::-1]
    s,e=plan['window_input0_halfopen'];a,b=plan['motif_input0_halfopen'];sa,sb=plan['seqlet_start'],plan['seqlet_end'];seq=plan['window_WT_sequence']
    assert np.array_equal(onehot,np.eye(4)[['ACGT'.index(base) for base in seq[sa-s:sb-s]]])
    cwm=g['contrib_scores'][:];native_ppm=g['sequence'][:];nseqlets=len(q['start'])
    np.savez_compressed(O/'native_pattern_source.npz',pattern_PPM=native_ppm,pattern_CWM=cwm,selected_seqlet_onehot_forward=onehot,selected_seqlet_contrib_stored=q['contrib_scores'][k])
position_scores=cwm[left:right].sum(1)
pair=dict(TOMTOM_target_pseudocount=.1,TOMTOM_target_nsites=100,TOMTOM_target_XML_smoothing_verified=True,FIMO_raw_reference_same_as_TOMTOM_input=True,fold=3,pattern=plan['pattern'],TF='EGR1',matrix_id='MA0162.5',n_seqlets=nseqlets,TOMTOM=row,XML_match=native.attrib,query_full_PPM=qp.tolist(),target_forward_PPM=tp.tolist(),target_oriented_PPM=oriented.tolist(),query_core_PPM=qp[left:right].tolist(),target_core_PPM=oriented[tl:tr].tolist(),query_core_1based_inclusive=[left+1,right],target_core_oriented_1based_inclusive=[tl+1,tr],query_core_consensus=''.join('ACGT'[i] for i in qp[left:right].argmax(1)),target_core_consensus=''.join('ACGT'[i] for i in oriented[tl:tr].argmax(1)),CWM_core_statistics=dict(sum=float(position_scores.sum()),max=float(position_scores.max()),min=float(position_scores.min())),source_sha256=sources,all_pattern_matches=plan['pattern_matches'],method='Original frozen full-pattern TOMTOM result; displayed logo is native10bp alignment overlap, not a new comparison on the cropped motif. Gene-specific FIMO positions are a separate coordinate system.')
dump(O/'TOMTOM_pair.json',pair)
meme=O/'EGR1_target.meme';meme.write_text('MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\nMOTIF MA0162.5 EGR1\nletter-probability matrix: alength= 4 w= 10 nsites= 100\n'+''.join(' '.join(f'{v:.8f}' for v in r)+'\n' for r in ppm)+'\n')
records=[f'>WT\n{seq}\n']
for i,v in enumerate(plan['variants']):
    t=v['window_offset'];records.append(f'>v{i}\n{seq[:t]+v["alternative"]+seq[t+1:]}\n')
fa=O/'WT_mutants.fa';fa.write_text(''.join(records))
scan=subprocess.run([cfg['fimo'],'--text','--thresh',str(cfg['fimo_p_threshold']),str(meme),str(fa)],capture_output=True,text=True,check=True)
assert 'truncat' not in scan.stderr.lower();(O/'WT_mutants_FIMO.tsv').write_text(scan.stdout)
hits=list(csv.DictReader([x for x in scan.stdout.splitlines() if x and not x.startswith('#')],delimiter='\t'))
def same(name):return [h for h in hits if h['sequence_name']==name and int(h['start'])-1==a-s and int(h['stop'])==b-s and h['strand']==plan['motif_strand_relative_to_input']]
WT=same('WT');assert WT and np.isclose(float(WT[0]['p-value']),plan['FIMO_p'])
h=json.loads((S/'highlight.json').read_text());vi=next(i for i,v in enumerate(plan['variants']) if v['position']==h['position'] and v['alternative']==h['alternative'])
h.update(FIMO_WT=WT,FIMO_same_site_loss=not bool(same(f'v{vi}')),FIMO_same_site_mutant=same(f'v{vi}'),overlapping_target_match_remaining=any(x['sequence_name']==f'v{vi}' and int(x['start'])-1<b-s and int(x['stop'])>a-s for x in hits))
pwm=np.log2(((counts+.5)/(counts+.5).sum(1,keepdims=True))/.25)
if plan['motif_strand_relative_to_input']=='-':pwm=pwm[::-1,::-1]
h['delta_PWM_same_site_bits']=float(pwm[h['position']-a,'ACGT'.index(h['alternative'])]-pwm[h['position']-a,'ACGT'.index(h['reference'])])
dump(O/'verified_highlight.json',h)
selection=dict(gene='PTGIS',fold=3,pattern=plan['pattern'],TF='EGR1',matrix_id='MA0162.5',candidate_id='candidate_003',method='Among the48 exactly perturbed shortlist, keep only real FIMO-hit candidates with a corresponding original TOMTOM pair q<.05 and the101bp maximum inside the site; choose highest actual core/flank peak contrast. User explicitly allows nonzero flanks and any TF.',pair=pair,highlight=h,contrast=json.loads((S/'contrast.json').read_text()),site_sequence_gene_oriented=seq[a-s:b-s],source_seqlet_50bp_region_input0_halfopen=[sa,sb],FIMO_motif_10bp_region_input0_halfopen=[a,b],display_window101bp_region_input0_halfopen=[s,e],scope='Known motif comparison is at aggregate-pattern level; FIMO confirms the selected gene sequence instance. Neither means that all50bp are a known TF binding core.',limitations='EGR1 is a significant alternative match; KLF16 has a lower TOMTOM q for this pattern. This supports an EGR1-like sequence comparison, not exclusive or experimentally validated EGR1 binding. Selected illustration, not a population estimate.')
dump(O/'verified_case.json',selection)
print(json.dumps(dict(status='verified',TOMTOM=row,FIMO_WT=WT,highlight=h,contrast=selection['contrast'],query_core=pair['query_core_1based_inclusive'],consensuses=[pair['query_core_consensus'],pair['target_core_consensus']],n_seqlets=nseqlets),indent=2))
