from pathlib import Path
import os, sys, json, hashlib, subprocess, shutil, csv, io
import pymupdf

DEST = Path('/mnt/HDD8TB/hongtao/end/pig/fig 5')
SOURCE = Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019/figures')
NAMES = ['EGR1_MA0162_5_F03_PTGIS_ISM_MATCHED_101BP_SMOOTH_SIGMA1p1', 'SP4_MA0685_2_F10_ZBTB7A_ISM_MATCHED_101BP_SMOOTH_SIGMA1p1']
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def write_json(p, obj): p.write_text(json.dumps(obj, ensure_ascii=False, indent=2)+'\n')
def pixels(p):
    with pymupdf.open(p) as doc: return doc[0].get_pixmap().samples
records = []
for name in NAMES:
    folder = DEST/name
    pdf = folder/(name+'.pdf')
    before_pdf = sha(pdf)
    data_hashes = {str(p.relative_to(folder)): sha(p) for p in (folder/'data').rglob('*') if p.is_file()}
    before_curves = sha(folder/'audit/curve_summary.tsv')
    code = folder/'code/render.py'
    text = code.read_text()
    old = ".spines[['left','bottom']].set_linewidth(.65)"
    assert text.count(old) == 2, name
    text = text.replace(old, ".spines[['left','bottom']].set_linewidth(1.3)")
    text = text.replace('qa=dict(output=', 'qa=dict(axis_spine_linewidth_pt=1.3,output=')
    text = text.replace('fig.canvas.draw();renderer=', "assert all(ax.spines[side].get_linewidth()==1.3 for ax in [curve,heat_ax] for side in ['left','bottom'])\nfig.canvas.draw();renderer=")
    code.write_text(text)
    trial = folder/'audit/axis_update_candidate.pdf'
    repro = folder/'audit/axis_update_reproduction.pdf'
    log = []
    try:
        for target in [trial,repro]:
            run = subprocess.run([sys.executable,str(folder/'plot.py'),'--output',str(target)],capture_output=True,text=True)
            log.append(run.stdout+'\n'+run.stderr)
            assert run.returncode==0,log[-1]
        assert pixels(trial)==pixels(repro), 'Redraw not identical'
        assert pixels(trial)!=pixels(pdf), 'Axis appearance unchanged'
        assert before_curves==sha(folder/'audit/curve_summary.tsv')
        assert all(sha(folder/rel)==value for rel,value in data_hashes.items())
        os.replace(trial,pdf)
    finally:
        trial.unlink(missing_ok=True)
        repro.unlink(missing_ok=True)
    shutil.copy2(pdf,SOURCE/pdf.name)
    (folder/'audit/axis_update_render.log').write_text('\n'.join(log))
    qa_path=folder/'audit/redraw_QA.json'
    qa=json.loads(qa_path.read_text()); qa.update(output=str(pdf),sha256=sha(pdf),bytes=pdf.stat().st_size)
    write_json(qa_path,qa)
    for filename in ['package.json','REPRODUCIBILITY.json']:
        p=folder/filename; meta=json.loads(p.read_text())
        meta['axis_spine_linewidth_pt']=1.3
        meta['axis_update_data_and_curves_unchanged']=True
        if filename=='REPRODUCIBILITY.json':
            for key in list(meta):
                if 'sha256' in key and 'pdf' in key.lower(): meta[key]=sha(pdf)
            meta['axis_update_independent_redraw_pixel_identical']=True
        write_json(p,meta)
    readme=folder/'README.md'
    readme.write_text(readme.read_text()+"\nAxis lines were updated to 1.3 pt; data, gain/loss values, smoothing and color scales are unchanged.\n")
    record=dict(folder=name,filename=pdf.name,before_pdf_sha256=before_pdf,source_sha256=sha(pdf),axis_spine_linewidth_pt=1.3,data_files_unchanged=True,curve_summary_unchanged=True,pixel_identical=True)
    write_json(folder/'audit/axis_update_QA.json',record)
    files=sorted(p for p in folder.rglob('*') if p.is_file() and p.name!='MANIFEST_SHA256.txt' and '__pycache__' not in p.parts)
    (folder/'MANIFEST_SHA256.txt').write_text(''.join(f'{sha(p)}  {p.relative_to(folder).as_posix()}\n' for p in files))
    files=[p for p in folder.rglob('*') if p.is_file() and '__pycache__' not in p.parts]
    record.update(files=len(files),bytes=sum(p.stat().st_size for p in files))
    records.append(record)
    print(json.dumps(record),flush=True)
updates={r['folder']:r for r in records}
for filename in ["figure_index.tsv",'FINAL_ISM_EXAMPLES.tsv']:
    path=DEST/filename
    with path.open(newline='') as f:
        reader=csv.DictReader(f,delimiter='\t'); fields=reader.fieldnames; rows=list(reader)
    for row in rows:
        if row['folder'] in updates:
            record=updates[row['folder']]
            for key in ['source_sha256','files','bytes','pixel_identical']:
                if key in row: row[key]=record[key]
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fields,delimiter='\t');writer.writeheader();writer.writerows(rows)
path=DEST/'PACKAGING_REPORT.json';report=json.loads(path.read_text())
for row in report['packages']:
    if row['folder'] in updates:
        row.update({k:updates[row['folder']][k] for k in ['source_sha256','files','bytes','pixel_identical']})
report['last_axis_update']='2026-10-21'
write_json(path,report)
write_json(DEST/'SELECTED_ISM_AXIS_UPDATE_20261021.json',dict(status='PASS',packages=records))
