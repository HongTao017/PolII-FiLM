"""Reference-inspired sequence names; retain all native topology/partitions.

The paper's IUPAC/repeat notation is borrowed, not its manual cluster decisions.
Homopolymer feature classes include reverse complements. A minimum repeat count
is measured from the raw, full-precision Cavener calls surrounding each frozen
core; '+' is used only when observed run lengths actually vary across members.
"""
from pathlib import Path
from collections import defaultdict,Counter
import sys,csv,json,math,hashlib
import numpy as np
ROOT=(Path(__file__).resolve().parents[1]/'data')
sys.path.insert(0,str(ROOT))
import revise_core_partitions as r
import add_two_single_page_overviews as single
from name_and_plot_matched44 import consensus,compact,reverse_complement,pdf_qa
from final_figure_small import save_pdf
import matplotlib.pyplot as plt

def run_bounds(sequence,start,end,letter):
 assert sequence[start:end]==letter*(end-start),(sequence[start:end],letter)
 lo,hi=start,end
 while lo>0 and sequence[lo-1]==letter:lo-=1
 while hi<len(sequence) and sequence[hi]==letter:hi+=1
 return lo,hi

def repeated_name(code):
 # Whole repeated multi-base units get count+unit, as in the paper's 2CT.
 # Internal single-letter runs use the already established unambiguous codes.
 for width in range(1,len(code)//2+1):
  if len(code)%width==0 and code[:width]*(len(code)//width)==code:
   return str(len(code)//width)+code[:width]
 return compact(code)

def orientation(i,med):
 if i==med:return False
 a,b=sorted([i,med]);return r.e.pairs[(a,b)]['direction']=='Reverse'

def prepare():
 old=json.loads((ROOT/'analysis/shared_core_partition_QA.json').read_text());assert old['status']=='PASS'
 groups=old['named_hit_bearing_multi_pattern_clades'];singles=old['independent_hits']
 oldnames={x['pattern_id']:x for x in csv.DictReader((ROOT/'analysis/matched44_sequence_names.tsv').open(),delimiter='\t')}
 source={x['key']:x for x in r.e.src['patterns']};evidence=[];classes=defaultdict(list)
 for u in groups:
  u['old_name']=u['name'];u['name']=repeated_name(u['core']['code'])
  u['scope']='shared_core_clade';med=r.e.idx[u['medoid']];arr=r.e.alignment(u['leaves'],med)
  start=u['core']['start']+50;end=u['core']['end']+50;code=u['core']['code']
  pure=len(set(code))==1 and code[0] in 'ACGT'
  u['homopolymer_class']=('T' if code[0] in 'AT' else 'G') if pure else ''
  for j,i in enumerate(u['leaves']):
   original=r.e.ps[i];raw=r.e.pp[i];rev=orientation(i,med);observed=np.where(np.isfinite(arr[j,:,0]))[0]
   shift=int(observed[0]);a=start-shift;b=end-shift
   raw_a,raw_b=(50-b,50-a) if rev else (a,b)
   seq=consensus(raw);run=0;runstart=runend=None
   if pure:
    oriented=consensus(np.nan_to_num(arr[j]));lo,hi=run_bounds(oriented,start,end,code[0]);run=hi-lo
    al,bl=lo-shift,hi-shift;runstart,runend=(50-bl,50-al) if rev else (al,bl)
    classes[u['homopolymer_class']].append(run)
   s=source[original['key']]
   evidence.append(dict(unit_id=u['group_id'],scope=u['scope'],pattern_id=original['id'],original_pattern=f"Fold {int(s['fold']):02d} / {s['pattern_id']}",pattern_key=s['key'],seqlets=s['seqlets'],old_unit_name=u['old_name'],reverse_complemented_for_clade_alignment=int(rev),shared_core_original_orientation=seq[raw_a:raw_b],shared_core_aligned_orientation=consensus(arr[j,start:end]),raw_core_start1=raw_a+1,raw_core_end1=raw_b,run_length=run,raw_run_start1=runstart+1 if runstart is not None else '',raw_run_end1=runend if runend is not None else '',source_h5=s['source_h5'],source_group=s['source_group']))
 for u in singles:
  i=u['node'];p=r.e.ps[i];s=source[p['key']];old=oldnames[p['id']];code=old['core_IUPAC_6bp'];u['old_name']=old['sequence_name']+'*'
  u['scope']='individual_pattern';u['name']=repeated_name(code)+'*';pure=len(set(code))==1 and code[0] in 'ACGT'
  u['homopolymer_class']=('T' if code[0] in 'AT' else 'G') if pure else ''
  start=int(old['core_start_1based'])-1;end=int(old['core_end_1based']);seq=consensus(r.e.pp[i]);run=0;lo=hi=None
  if pure:
   lo,hi=run_bounds(seq,start,end,code[0]);run=hi-lo;classes[u['homopolymer_class']].append(run)
  evidence.append(dict(unit_id=u['group_id'],scope=u['scope'],pattern_id=p['id'],original_pattern=f"Fold {int(s['fold']):02d} / {s['pattern_id']}",pattern_key=s['key'],seqlets=s['seqlets'],old_unit_name=u['old_name'],reverse_complemented_for_clade_alignment=0,shared_core_original_orientation=code,shared_core_aligned_orientation=code,raw_core_start1=start+1,raw_core_end1=end,run_length=run,raw_run_start1=lo+1 if lo is not None else '',raw_run_end1=hi if hi is not None else '',source_h5=s['source_h5'],source_group=s['source_group']))
 definitions={}
 for letter,lengths in classes.items():
  name=f'{min(lengths)}{letter}'+('+' if max(lengths)>min(lengths) else '')
  definitions[letter]=dict(name=name,min_repeats=min(lengths),max_repeats=max(lengths),members=len(lengths),all_run_lengths=lengths,reverse_complement_symbol='A' if letter=='T' else 'C')
 for u in groups+singles:
  if u['homopolymer_class']:u['name']=definitions[u['homopolymer_class']]['name']+('*' if u['scope']=='individual_pattern' else '')
 # Name each native homopolymer clade from its own observed repeat lengths.
 # A broad reverse-complement feature class does not replace clade-specific
 # minimum counts: old 4T-01 is 5T+, while old 4T-02 is 13T+.
 for u in groups:
  if u['homopolymer_class']:
   lengths=[row['run_length'] for row in evidence if row['unit_id']==u['group_id']]
   assert lengths and min(lengths)>0
   u['name']=f"{min(lengths)}{u['homopolymer_class']}"+('+' if max(lengths)>min(lengths) else '')
   u['branch_repeat_evidence']=dict(min_repeats=min(lengths),max_repeats=max(lengths),all_run_lengths=lengths)
 units=sorted(groups+singles,key=lambda u:u['rank']);byunit={u['group_id']:u for u in units}
 for row in evidence:
  u=byunit[row['unit_id']];row['new_display_name']=u['name'];row['homopolymer_feature_class']=u['homopolymer_class']
  rawbase=row['shared_core_original_orientation'][0]
  row['reverse_complemented_for_feature_name']=int(bool(u['homopolymer_class']) and rawbase in 'AC')
 return old,groups,singles,units,evidence,definitions

def draw(groups,units):
 records=[]
 g=single.balanced_geometry(units);shadows=r.annotate_units(g,units)
 for u in shadows:u['angle']=float(g['theta'][u['node']])
 print('layout_reference_names_29',flush=True)
 scale,rot=single.single_layout(shadows)
 path=ROOT/'figures/motifStack_paper_style_ALL29_SINGLE_PAGE_5x4.pdf'
 fig,qa=r.draw_tree(g,shadows,shadows,scale,rot,'a',1,1);save_pdf(fig,path);plt.close(fig)
 qa.update(pdf_qa(path,[u['name'] for u in shadows]));qa.update(single.verify_pdf(path))
 records.append(dict(path=str(path),annotation_count=29,unique_display_names=len(set(u['name'] for u in shadows)),scale=scale,rotation=rot,QA=qa))
 oldlayout=json.loads((ROOT/'analysis/shared_core_partition_QA.json').read_text())['overview']
 g13=r.geometry(groups);s13=r.annotate_units(g13,units);labels=[u for u in s13 if u['n']>1]
 scale13,rot13=oldlayout['scale'],oldlayout['rotation']
 # Renaming can change label width, so verify before keeping original layout.
 fig,ax=r.setup();texts=[]
 for u in labels:
  t=ax.text(0,0,u['name'],fontsize=12,color='black',rotation_mode='anchor',va='center');r.place(t,u['angle']+math.radians(rot13),scale13);texts.append(t)
 qa13=r.label_check(fig,texts);plt.close(fig)
 if any(qa13.values()):scale13,rot13=single.single_layout(labels)
 path13=ROOT/'figures/motifStack_paper_style_13CLADES_SINGLE_PAGE_5x4.pdf'
 fig,qa=r.draw_tree(g13,s13,labels,scale13,rot13,'a',1,1);save_pdf(fig,path13);plt.close(fig)
 qa.update(pdf_qa(path13,[u['name'] for u in labels]));qa.update(single.verify_pdf(path13))
 records.append(dict(path=str(path13),annotation_count=13,unique_display_names=len(set(u['name'] for u in labels)),scale=scale13,rotation=rot13,QA=qa))
 return records

def main():
 before={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'figures').glob('*.pdf') if not p.name.startswith('motifStack_paper_style_')}
 old,groups,singles,units,evidence,defs=prepare()
 assert len(units)==29 and len(groups)==13 and len(singles)==16
 print('NAMES',json.dumps([dict(id=u['group_id'],old=u['old_name'],new=u['name']) for u in units]),flush=True)
 print('REPEAT_EVIDENCE',json.dumps(defs),flush=True)
 records=draw(groups,units)
 match=list(csv.DictReader((ROOT/'analysis/original_TOMTOM_matches.tsv').open(),delimiter='\t'))
 evidenceby={x['pattern_id']:x for x in evidence};assert len(evidenceby)==65
 for row in evidence:
  hits=[x for x in match if x['query']==row['pattern_id']]
  row.update(matched=int(bool(hits)),TFs=';'.join(x['TF'] for x in hits),JASPAR_ids=';'.join(x['matrix_id'] for x in hits),q_values=';'.join(x['q'] for x in hits))
 r.write_tsv(ROOT/'analysis/paper_style_naming_members65.tsv',evidence)
 r.write_tsv(ROOT/'analysis/paper_style_naming_matched44.tsv',[x for x in evidence if x['matched']])
 r.write_tsv(ROOT/'analysis/paper_style_naming_TOMTOM79.tsv',[dict(new_display_name=evidenceby[x['query']]['new_display_name'],original_pattern=evidenceby[x['query']]['original_pattern'],**x) for x in match])
 mapping=[]
 for u in units:
  members=[x for x in evidence if x['unit_id']==u['group_id']]
  repeat=u.get('branch_repeat_evidence',{})
  mapping.append(dict(unit_id=u['group_id'],old_name=u['old_name'],new_name=u['name'],scope=u['scope'],conserved_IUPAC=u['core']['code'] if u['core'] else members[0]['shared_core_original_orientation'],patterns=u['n'],matched=u['hits'],branch_min_repeats=repeat.get('min_repeats',''),branch_max_repeats=repeat.get('max_repeats',''),branch_run_lengths=';'.join(str(x) for x in repeat.get('all_run_lengths',[])),original_patterns=';'.join(x['original_pattern'] for x in members)))
 r.write_tsv(ROOT/'analysis/paper_style_name_dictionary29.tsv',mapping)
 assert sum(x['matched'] for x in evidence)==44 and len(match)==79
 for path,sha in before.items():assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==sha
 assert hashlib.sha256((ROOT/'analysis/motifStack_hclust.rds').read_bytes()).hexdigest()==r.TREE_SHA
 report=dict(status='PASS',figures=records,repeat_feature_classes=defs,native_tree_sha256=r.TREE_SHA,source_h5_sha256_verified=r.e.source_checks,source_matrices_and_matches_unchanged=True,partitions_unchanged=True,annotation_units=29,clades=13,independent_patterns=16,matched_patterns=44,TOMTOM_records=79,all_old_figures_preserved=before,
 method=dict(notation='IUPAC conserved sequence, minimum repeated unit, optional longer repeats (+); reverse complements represented by a shared feature name',interpretation='Descriptive sequence features without functional family validation',group_cores='unchanged validated aligned 4-12 bp cores from the frozen clade partition',individual_cores='unchanged maximum-information 6bp individual features, explicitly marked *',repeat_definition='maximal exact Cavener single-base run containing the frozen core in each original PPM; clade prefix is the minimum observed within that native clade; independent homopolymer labels retain the broader A/T or C/G pooled feature-class lower bound and *; + only if the relevant observed run lengths vary',mixed_cores='IUPAC order retained; exact repeated whole multi-base units compressed when present, otherwise single-letter repeat notation; no + added without extension evidence',same_name='shared short-sequence feature, not a claim of same native branch or full motif; repeated names retained at distinct tree positions'))
 (ROOT/'analysis/paper_style_naming_QA.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
 note='''# Sequence-based partition labels

The original 314-pattern tree, 13 partitions and 16 independent items are retained. Labels summarize IUPAC cores and repeat lengths without changing tree membership.

## Repeat labels

'''
 for letter,d in defs.items():note+=f"- `{d['name']}`: {d['members']} patterns. Native Cavener consensus has {d['min_repeats']}-{d['max_repeats']} consecutive repeats around the core. {letter} and reverse-complement symbol {d['reverse_complement_symbol']} share a sequence-feature label while retaining separate tree positions.\n"
 note+='''Repeat labels use the shortest observed run in each partition. The two T-rich branches have run lengths 14/5/5 and 18/15/13 and are labeled 5T+ and 13T+. Reverse-complement orientation is used for A-rich members. The plus sign indicates longer observed runs and is not a significance marker. Shared labels do not establish distinct TF families. Mixed cores retain supported contiguous IUPAC features; noncontiguous runs are not joined. Independent-pattern stars identify individual patterns rather than significance. Blue/pink branches indicate negative/positive PC1 attribution, gray sectors indicate partitions and outer-ring purple indicates match status. Partition cores use 4-12-bp descriptive rules; independent patterns use six-base windows.
'''
 for rec in records:note+=f"- `{rec['path']}`: one page, {rec['annotation_count']} annotation positions and {rec['unique_display_names']} distinct display names.\n"
 note+='''- `analysis/paper_style_name_dictionary29.tsv`: partition-label correspondence.
- `analysis/paper_style_naming_matched44.tsv`: 44 matched patterns, original identifiers, annotations and core coordinates.
- `analysis/paper_style_naming_members65.tsv`: 65 partition members.
- `analysis/paper_style_naming_TOMTOM79.tsv`: 79 retained database matches.
- `analysis/paper_style_naming_QA.json`: naming parameters and records.

## Label correspondence

| Position | Previous label | Current label | Range |
|---|---|---|---|
'''
 for m in mapping:note+=f"|{m['unit_id']}|{m['old_name']}|{m['new_name']}|{"Partition" if m['scope']=='shared_core_clade' else "Independent pattern"}|\n"
 (ROOT/'analysis/paper_style_naming_guide.md').write_text(note,encoding='utf-8')
 for path in [ROOT/'analysis/current_figure_guide.md',ROOT/"README_motifStack_redraw.md"]:
  text=path.read_text();marker="\n\n# Sequence-based partition labels";text=text.split(marker)[0];path.write_text(text+'\n\n'+note,encoding='utf-8')
 print(json.dumps(dict(status='PASS',figures=records,classes=defs,dictionary=mapping)),flush=True)

if __name__=='__main__':main()
