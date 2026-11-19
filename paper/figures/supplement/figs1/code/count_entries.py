"""Count finite original log2(PI) transcript entries, before any imputation.

Laboratory grouping uses the existing source-author metadata as its label.
The count is neither a PI sum nor a count of positive log2(PI) values.
"""
from pathlib import Path
import re,json,hashlib
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
META={'gene','transcript','gene_transcript','chrom','strand','tss_bed0','tes_bed0'}
FILES={'GRO-seq':'Gro_seq_human_transcript_log2_pi_matrix.tsv',
       'MNET-seq':'Mnet_seq_human_transcript_log2_pi_matrix.tsv',
       'PRO-seq':'PRO_seq_human_transcript_log2_pi_matrix.tsv'}
ASSAY_ORDER={name:i for i,name in enumerate(FILES)}
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()

def natural(value):
    return tuple((1,int(part)) if part.isdigit() else (0,part.casefold())
                 for part in re.split(r'(\d+)',str(value)) if part)

def count_entries():
    sources=json.loads((ROOT/'provenance/source_files.json').read_text(encoding='utf-8'))
    for record in sources:assert sha(ROOT/record['package_file'])==record['sha256'],record['package_file']
    metadata=pd.read_csv(ROOT/'data/pi_sample_metadata_ordered.tsv',sep='\t')
    # The original heatmap metadata spells the same assay as mNET-seq.
    metadata['assay']=metadata['assay'].replace({'mNET-seq':'MNET-seq'})
    assert len(metadata)==105 and not metadata.duplicated(['assay','sample']).any()
    lookup=metadata.set_index(['assay','sample'])
    union=set();rows=[];matrix_audit=[]
    for assay,filename in FILES.items():
        path=ROOT/'data/raw'/filename
        table=pd.read_csv(path,sep='\t',na_values=['NA','NaN',''],low_memory=False)
        columns=[c for c in table if c not in META]
        keys=(table['gene_transcript'].astype(str) if 'gene_transcript' in table else
              table['gene'].astype(str)+'|'+table['transcript'].astype(str))
        assert keys.notna().all() and keys.is_unique
        union.update(keys)
        numbers=table[columns].apply(pd.to_numeric,errors='raise').to_numpy(dtype=float)
        finite=np.isfinite(numbers);nfinite=finite.sum(axis=0)
        for colidx,sample in enumerate(columns):
            info=lookup.loc[[(assay,sample)]].iloc[0]
            assert info['source_file']==filename
            author=str(info['author']);cell=str(info['cell_line'])
            study=re.search(r'GSE\d+',sample);study=study.group() if study else ''
            rep=re.search(r'(?:^|_)rep(\d+)(?:_|$)',sample,re.IGNORECASE)
            replicate=int(rep.group(1)) if rep else 0
            rest=sample[len(author)+1:] if sample.startswith(author+'_') else sample
            label=f'{author}_{assay}_{rest}'
            rows.append({'sample':sample,'assay':assay,'author':author,'laboratory_group':author,
                         'study':study,'cell_line':cell,'replicate':replicate,
                         'valid_transcript_entries':int(nfinite[colidx]),
                         'missing_or_nonfinite_entries_within_source':int((~finite[:,colidx]).sum()),
                         'zero_log2_PI_entries_counted_as_valid':int((numbers[:,colidx]==0).sum()),
                         'negative_log2_PI_entries_counted_as_valid':int((finite[:,colidx] & (numbers[:,colidx]<0)).sum()),
                         'source_transcript_entries':len(table),'source_file':filename,
                         'display_label':label})
        matrix_audit.append({'assay':assay,'filename':filename,'sha256':sha(path),
                             'transcript_entries':len(table),'samples':len(columns),
                             'minimum_valid_entries':int(nfinite.min()),'maximum_valid_entries':int(nfinite.max())})
    assert len(rows)==105 and len(union)==16696
    assert set((r['assay'],r['sample']) for r in rows)==set(lookup.index)
    rows.sort(key=lambda r:(r['author'].casefold(),natural(r['study']),ASSAY_ORDER[r['assay']],
                            natural(r['cell_line']),r['replicate'],natural(r['sample'])))
    for i,row in enumerate(rows,1):row['display_order']=i
    counts=pd.DataFrame(rows)
    assert counts.valid_transcript_entries.min()==9197 and counts.valid_transcript_entries.max()==14338
    assert counts.author.nunique()==34 and counts.display_label.is_unique
    counts.to_csv(ROOT/'data/sample_valid_transcript_counts.tsv',sep='\t',index=False)
    pd.DataFrame({'transcript_entry':sorted(union)}).to_csv(ROOT/'data/transcript_union.tsv',sep='\t',index=False)
    for group,column in [('laboratory','author'),('assay','assay')]:
        summary=counts.groupby(column,sort=False).agg(samples=('sample','size'),
                min_valid_entries=('valid_transcript_entries','min'),
                max_valid_entries=('valid_transcript_entries','max'),
                mean_valid_entries=('valid_transcript_entries','mean'))
        summary.to_csv(ROOT/f'data/{group}_summary.tsv',sep='\t')
    report={'status':'PASS','sample_count':len(counts),'laboratory_author_groups':int(counts.author.nunique()),
            'assay_counts':{k:int(v) for k,v in counts.assay.value_counts().items()},
            'transcript_union_count':len(union),'minimum_valid_entries':int(counts.valid_transcript_entries.min()),
            'maximum_valid_entries':int(counts.valid_transcript_entries.max()),
            'count_definition':'Number of finite original log2(PI) transcript-entry values, before imputation; zero and negative values included.',
            'entry_identifier':'gene|transcript, or existing gene_transcript key',
            'laboratory_group_definition':'Existing source-author field used as a laboratory/source grouping label; no new institutional affiliation inferred.',
            'ordering':'Source author alphabetically; within author, GSE natural order, assay GRO/MNET/PRO, cell line natural order, replicate number, sample ID.',
            'raw_values_modified':False,'imputation_used':False,'matrices':matrix_audit}
    (ROOT/'audit/counts_QA.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    return counts,report

if __name__=='__main__':
    _,report=count_entries();print(json.dumps(report,ensure_ascii=False,indent=2))
