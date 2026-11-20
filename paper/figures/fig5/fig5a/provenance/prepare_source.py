from pathlib import Path
import sys,json,csv,hashlib,shutil,subprocess,os
import numpy as np,h5py
ROOT=Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019')
DEST=Path('/mnt/HDD8TB/hongtao/end/pig/fig 5')
ISM=DEST/'SP4_MA0685_2_F10_ZBTB7A_ISM_MATCHED_101BP_SMOOTH_SIGMA1p1'
R=Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motif/folds/fold_10')
plan=json.loads((ISM/'data/perturbation/frozen_plan.json').read_text())
cfg=json.loads((R/'experiment_config.json').read_text());man=json.loads((R/'run_manifest.json').read_text());completion=json.loads((R/'complete.json').read_text())
used={}
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for block in iter(lambda:f.read(2**20),b''):h.update(block)
    return h.hexdigest()
def check(p,expected):
    assert sha(p)==expected,str(p)
    used[str(p)]=expected
for name in ['selected_genes.tsv','modisco_input_order.npy','modisco/tss_sequence.h5']:check(R/name,completion['outputs_sha256'][name])
source=Path(cfg['source_root'])/'data/sequence_prepared/tss_sequence.npy'
metadata=Path(cfg['model_root'])/'data/genes_repartitioned.tsv'
for p in [source,metadata]:check(p,man['original_inputs_sha256'][str(p)])
selected=list(csv.DictReader((R/'selected_genes.tsv').open(),delimiter='\t'))
rows=list(csv.DictReader(metadata.open(),delimiter='\t'))
order=np.load(R/'modisco_input_order.npy')
with h5py.File(R/'modisco/tss_sequence.h5','r') as h5:
    group=h5['pos_patterns/pattern_0'];q=group['seqlets'];k=plan['seqlet']
    i=int(order[int(q['example_idx'][k])]);index=int(selected[i]['gene_index']);s,e=int(q['start'][k]),int(q['end'][k]);rev=bool(q['is_revcomp'][k])
    assert (i,index,s,e,rev)==(64,12871,55,105,False) and selected[i]['gene']=='ZBTB7A'
    stored_dna=q['sequence'][k].copy();stored_cwm=q['contrib_scores'][k].copy()
dna=np.array(np.load(source,mmap_mode='r')[index])
attrs=[]
for seed in cfg['seeds']:
    p=R/'attributions'/f'seed_{seed}/hypothetical.npy';info=json.loads(p.with_name('complete.json').read_text());check(p,info['output_sha256']['hypothetical.npy'])
    attrs.append(np.load(p,mmap_mode='r')[i,:,:600].copy())
mean=np.mean(attrs,axis=0).T;actual=dna*mean
assert np.array_equal(dna[s:e],stored_dna)
assert np.allclose(actual[s:e],stored_cwm,rtol=2e-5,atol=1e-7),float(abs(actual[s:e]-stored_cwm).max())
seq=''.join('ACGT'[int(v.argmax())] if v.sum()==1 else 'N' for v in dna)
assert seq[50:151]==plan['window_WT_sequence']
m=rows[index];left=int(m['tss_sequence_requested_start0']);right=int(m['tss_sequence_requested_end0']);tss=int(m['tss_bed0'])
coords=np.arange(left,right) if m['strand']=='+' else np.arange(right-1,left-1,-1)
tssidx=tss-left if m['strand']=='+' else right-1-tss
assert len(coords)==600
info=dict(gene='ZBTB7A',gene_index=index,selected_row=i,fold=10,pattern=plan['pattern'],seqlet=k,seqlet_input0_halfopen=[s,e],motif_input0_halfopen=plan['motif_input0_halfopen'],seqlet_is_revcomp=rev,source_seqlet_sequence=seq[s:e],sequence=seq,tss_input_index=tssidx,metadata=m,seeds=cfg['seeds'],method=cfg['method'],target=cfg['attribution_target'],source_sha256=used,observed_base_definition='onehot * mean hypothetical attribution; summed over four bases at each position',single_gene_seqlet_matches_modisco=True,maximum_seqlet_difference=float(abs(actual[s:e]-stored_cwm).max()),same_instance_as_final_ISM=True,coordinate_convention='0-based input positions; intervals half-open')
reports=[]
for suffix,panel in [('TSS600_IMPORTANCE','full_tss'),('SEQLET50_CONTRIBUTION','seqlet50')]:
    name='SP4_MA0685_2_F10_ZBTB7A_'+suffix;folder=DEST/name
    for sub in ['data/fonts','code','audit','provenance']:(folder/sub).mkdir(parents=True,exist_ok=True)
    shutil.copy2(ROOT/'draw_zbtb7a_attribution.py',folder/'plot.py')
    shutil.copy2(ROOT/'prepare_zbtb7a_attribution.py',folder/'provenance/prepare_source.py')
    shutil.copy2(ROOT/'final_figure_small.py',folder/'code/final_figure_small.py')
    for f in (ISM/'data/fonts').glob('*.ttf'):shutil.copy2(f,folder/'data/fonts'/f.name)
    for src,dst in [(ISM/'data/figure_selection.json','ISM_figure_selection.json'),(ISM/'data/TOMTOM_pair.json','TOMTOM_pair.json')]:shutil.copy2(src,folder/'data'/dst)
    np.savez_compressed(folder/'data/selected_real_attribution.npz',onehot=dna,hypothetical=mean,actual=actual,genomic_position0=coords,seqlet_contrib_scores=stored_cwm,seed_hypothetical=np.stack(attrs).transpose(0,2,1))
    (folder/'data/selection.json').write_text(json.dumps(info,indent=2))
    np.savetxt(folder/'data/attribution_600bp.tsv',np.column_stack([np.arange(600),coords+1,list(seq),actual.sum(1)]),fmt='%s',delimiter='\t',header='input_position0\tgenomic_position1\tbase\tsigned_observed_contribution',comments='')
    package=dict(filename=name+'.pdf',panel=panel,description="ZBTB7A/SP4 complete TSS attribution" if panel=='full_tss' else "ZBTB7A/SP4 50-bp source-seqlet contributions",gene='ZBTB7A',fold=10,pattern=plan['pattern'],seqlet=75)
    (folder/'package.json').write_text(json.dumps(package,ensure_ascii=False,indent=2))
    (folder/'requirements.txt').write_text('numpy==2.2.6\nmatplotlib==3.10.9\nPyMuPDF==1.28.2\n')
    (folder/'README.md').write_text("# ZBTB7A / SP4 attribution\n\nRun `python plot.py` or provide an output PDF with --output.\n\n"+package['description']+". The figure and final 101-bp mutagenesis example use the same gene and fold-10 seqlet. Values are signed PC1 contributions. Complete input positions are 0-599; the source seqlet spans 55-104 and the SP4-like match spans 96-104. Contributions are observed-base one-hot values multiplied by five-seed mean hypothetical attributions, without smoothing or normalization. This is one seqlet, not the 98-seqlet average pattern. Frozen data, coordinates and source hashes are retained.\n")
    trials=[folder/'audit/trial1.pdf',folder/'audit/trial2.pdf']
    import pymupdf
    logs=[]
    for trial in trials:
        run=subprocess.run([sys.executable,str(folder/'plot.py'),'--output',str(trial)],text=True,capture_output=True)
        logs.append(run.stdout+'\n'+run.stderr);assert run.returncode==0,logs[-1]
    with pymupdf.open(trials[0]) as d1,pymupdf.open(trials[1]) as d2:assert d1[0].get_pixmap().samples==d2[0].get_pixmap().samples
    output=folder/package['filename'];os.replace(trials[0],output);trials[1].unlink()
    shutil.copy2(output,ROOT/'figures'/output.name)
    qa_path=folder/'audit/redraw_QA.json';qa=json.loads(qa_path.read_text());qa.update(output=str(output),PDF_sha256=sha(output));qa_path.write_text(json.dumps(qa,indent=2))
    (folder/'audit/redraw.log').write_text('\n'.join(logs))
    (folder/'REPRODUCIBILITY.json').write_text(json.dumps(dict(status='PASS',independent_redraw_pixel_identical=True,PDF_sha256=sha(output),source_files_hash_verified=True,single_gene_seqlet_verified_against_native_modisco=True,same_instance_as_final_ISM=True,data_paths_relative=True),indent=2))
    files=sorted(p for p in folder.rglob('*') if p.is_file() and p.name!='MANIFEST_SHA256.txt' and '__pycache__' not in p.parts)
    (folder/'MANIFEST_SHA256.txt').write_text(''.join(f'{sha(p)}  {p.relative_to(folder).as_posix()}\n' for p in files))
    files=[p for p in folder.rglob('*') if p.is_file() and '__pycache__' not in p.parts]
    record=dict(folder=name,description=package['description'],filename=output.name,source_sha256=sha(output),files=len(files),bytes=sum(p.stat().st_size for p in files),pixel_identical=True)
    reports.append(record);print(json.dumps(record),flush=True)
for filename in ["figure_index.tsv"]:
    path=DEST/filename
    with path.open(newline='') as f:
        reader=csv.DictReader(f,delimiter='\t');fields=reader.fieldnames;rows=list(reader)
    rows=[row for row in rows if row['folder'] not in {r['folder'] for r in reports}]+[{k:r[k] for k in fields} for r in reports]
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fields,delimiter='\t');writer.writeheader();writer.writerows(rows)
path=DEST/'PACKAGING_REPORT.json';report=json.loads(path.read_text());report['packages']=[r for r in report['packages'] if r['folder'] not in {p['folder'] for p in reports}]+reports;report['count']=len(report['packages']);path.write_text(json.dumps(report,ensure_ascii=False,indent=2))
(DEST/'ZBTB7A_ATTRIBUTION_PANELS_20261021.json').write_text(json.dumps(dict(status='PASS',panels=reports,TSS_index0=tssidx,source_attribution_range=[float(actual.sum(1).min()),float(actual.sum(1).max())]),ensure_ascii=False,indent=2))
