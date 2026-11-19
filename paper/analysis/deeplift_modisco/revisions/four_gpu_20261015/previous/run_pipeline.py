"""Separate frozen-weight DeepLIFT -> TF-MoDISco -> TOMTOM -> ISM pilot.

One process, exclusive run lock, idle GPU6 only. No training or old output writes.
Any numerical/protocol/annotation error stops; no automatic threshold relaxation.
"""
import argparse
import fcntl
import importlib.metadata as md
import json
import os
from pathlib import Path
import shutil
import subprocess
import time
import traceback

import numpy as np
import pandas as pd
import torch

from deep_lift_core import (ConditionalDNA, REGIONS, LENGTHS, OFFSETS, digest,
    edit_physical, explain_pair, gene_inputs, load_original, physical_map,
    shuffled_references, stable_seed, write_json)

ROOT = Path(__file__).resolve().parent


def state(stage, **fields):
    record = dict(pid=os.getpid(), status='running', stage=stage,
        updated_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), **fields)
    write_json(ROOT/'controller_state.json', record)
    print(json.dumps(record), flush=True)


def wait_gpu(config):
    while True:
        gpu = subprocess.check_output(['nvidia-smi', '--query-gpu=index,uuid,memory.used', '--format=csv,noheader,nounits'], text=True)
        row = next(x.split(',') for x in gpu.splitlines() if int(x.split(',')[0]) == config['physical_gpu'])
        processes = subprocess.check_output(['nvidia-smi', '--query-compute-apps=gpu_uuid,pid', '--format=csv,noheader,nounits'], text=True)
        own_gpu_processes = [int(x.split(',')[1]) for x in processes.splitlines() if x.startswith(row[1].strip()+',')]
        if (own_gpu_processes and set(own_gpu_processes) == {os.getpid()}) or (not own_gpu_processes and int(row[2]) < 64):
            return
        state('waiting_gpu6')
        time.sleep(30)


def checkpoint_hashes(config):
    train = json.loads((Path(config['model_root'])/'training_config.json').read_text())
    base = Path(config['model_root'])/'result'/train['run_name']/f'fold_{config["fold"]:02d}'
    return {str(seed): {name: digest(base/f'seed_{seed}'/name) for name in ('best.pt', 'preprocessing.npz', 'run_spec.json', 'complete.json')} for seed in config['seeds']}


def prepare(config):
    assert not (ROOT/'run_manifest.json').exists(), 'Run already prepared: no implicit resume or overwrite'
    assert json.loads((ROOT/'preflight_complete.json').read_text())['status'] == 'PASS'
    selection = Path(config['selection_root'])/'selected_genes.tsv'
    source_complete = json.loads((Path(config['selection_root'])/'complete.json').read_text())
    assert source_complete['status'] == 'complete'
    assert digest(selection) == source_complete['outputs_sha256']['selected_genes.tsv']
    shutil.copy2(selection, ROOT/'selected_genes.tsv')
    shutil.copy2(Path(config['selection_root'])/'motifs.json', ROOT/'motifs.json')
    selected = pd.read_csv(selection, sep='\t')
    metadata = pd.read_csv(Path(config['model_root'])/'data/genes_repartitioned.tsv', sep='\t')
    assert len(selected) == config['genes'] and selected.gene.is_unique and selected.dna_overlap_group.is_unique
    assert selected.fold.eq(config['fold']).all()
    assert metadata.iloc[selected.gene_index].gene.tolist() == selected.gene.tolist()
    (ROOT/'references').mkdir()
    source = Path(config['source_root'])/'data/bpnet_prepared'
    input_hashes={str(source/name):digest(source/name) for name in
        [*(f'{r}.npy' for r in REGIONS),'regulatory_log10_unscaled.npy','feature_manifest.tsv']}
    input_hashes[str(Path(config['model_root'])/'data/genes_repartitioned.tsv')]=digest(Path(config['model_root'])/'data/genes_repartitioned.tsv')
    (ROOT/'environment.lock.txt').write_text(subprocess.check_output([str(Path(__import__('sys').executable)),'-B','-m','pip','freeze','--all'],text=True))
    arrays = {r: np.load(source/f'{r}.npy', mmap_mode='r') for r in REGIONS}
    reference_hashes = {}
    for i, row in enumerate(selected.itertuples()):
        actual = np.concatenate([arrays[r][row.gene_index] for r in REGIONS], 0).T
        references = shuffled_references(metadata.iloc[row.gene_index], actual, config['references_per_gene'], config['random_seed'])
        name = f'references/gene_{i:03d}.npz'
        np.savez_compressed(ROOT/name, references=references)
        reference_hashes[name] = digest(ROOT/name)
        if (i+1) % 20 == 0:
            state('preparing_physical_consistent_references', genes=i+1, total_genes=len(selected))
    manifest = dict(config_sha256=digest(ROOT/'experiment_config.json'), selected_genes_sha256=digest(ROOT/'selected_genes.tsv'),
        motif_database_sha256=digest(ROOT/'motifs.json'),
        original_checkpoint_files=checkpoint_hashes(config), references_sha256=reference_hashes,
        original_inputs_sha256=input_hashes,
        code_sha256={p.name: digest(p) for p in ROOT.glob('*.py')},
        dependencies={k: md.version(k) for k in ('torch', 'numpy', 'modisco', 'tangermeme', 'h5py', 'leidenalg', 'igraph')},
        reference_policy='Dinucleotide shuffle each unique contiguous physical DNA block; synchronize all overlapping copies. Preserve N and padding. With overlaps, union blocks, not each individual window, preserve dinucleotide composition.',
        selection='Frozen old fold-1 test genes, no Y/ChIP/effect selection', optimizer_updates=0, bn_refits=0, preprocessing_refits=0)
    write_json(ROOT/'run_manifest.json', manifest)
    return selected, metadata, manifest


def attribute(config, selected):
    references_sha256=json.loads((ROOT/'run_manifest.json').read_text())['references_sha256']
    for k, seed in enumerate(config['seeds']):
        wait_gpu(config)
        folder = ROOT/'attributions'/f'seed_{seed}'
        folder.mkdir(parents=True)
        model, processor, identity = load_original(config, seed, 'cuda')
        with np.load(Path(identity['folder'])/'preprocessing.npz') as prep:
            target_mean, target_std = prep['target_mean'].item(), prep['target_std'].item()
        values = np.zeros((len(selected), 4, int(OFFSETS[-1])), dtype=np.float32)
        deltas = np.zeros((len(selected), config['references_per_gene']), dtype=np.float32)
        errors = deltas.copy()
        for i, row in enumerate(selected.itertuples()):
            actual, regulation = gene_inputs(config, row.gene_index, processor, 'cuda')
            wrapper = ConditionalDNA(model, regulation, target_mean, target_std)
            ref_path=f'references/gene_{i:03d}.npz'
            assert digest(ROOT/ref_path)==references_sha256[ref_path]
            refs = np.load(ROOT/ref_path)['references']
            accumulated = np.zeros((4, int(OFFSETS[-1])), dtype=np.float64)
            for j, reference in enumerate(refs):
                reference = torch.from_numpy(reference[None].astype(np.float32)).cuda()
                hypo, observed, delta, error = explain_pair(wrapper, actual, reference)
                if not torch.isfinite(hypo).all() or not torch.isfinite(error).all():
                    raise FloatingPointError((seed, row.gene, j, 'nonfinite DeepLIFT'))
                limit = config['completeness_atol']+config['completeness_rtol']*delta.abs()
                assert (error.abs() <= limit).all(), (seed, row.gene, j, error.item(), delta.item())
                accumulated += hypo[0].cpu().numpy()
                deltas[i,j], errors[i,j] = delta.item(), error.item()
            values[i] = accumulated / len(refs)
            if (i+1) % 5 == 0 or i == 0:
                state('DeepLIFT', complete_seeds=k, seed=seed, gene=row.gene, genes=i+1,
                    total_genes=len(selected), references_per_gene=len(refs), max_completeness_error=float(np.max(np.abs(errors[:i+1]))))
        np.save(folder/'hypothetical.npy', values)
        np.savez_compressed(folder/'completeness.npz', prediction_delta=deltas, error=errors)
        write_json(folder/'complete.json', dict(status='complete', seed=seed, genes=len(selected), **identity,
            method=config['method'], output_sha256={name:digest(folder/name) for name in ('hypothetical.npy', 'completeness.npz')}))
        del model, wrapper
        torch.cuda.empty_cache()


MEME_HEADER = 'MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n'


def meme_entry(name, matrix, nsites=20):
    p = np.asarray(matrix, dtype=float)
    assert p.shape[1] == 4 and np.isfinite(p).all() and (p >= 0).all()
    p = (p+1e-8)/(p+1e-8).sum(1, keepdims=True)
    return f'MOTIF {name}\nletter-probability matrix: alength= 4 w= {len(p)} nsites= {nsites}\n'+''.join(' '.join(f'{x:.8f}' for x in row)+'\n' for row in p)+'\n'


def discover(config, selected):
    import h5py
    from modiscolite.tfmodisco import TFMoDISco
    from modiscolite.io import save_hdf5
    source = Path(config['source_root'])/'data/bpnet_prepared'
    mean = np.mean([np.load(ROOT/'attributions'/f'seed_{s}'/'hypothetical.npy') for s in config['seeds']], 0)
    order = np.array(sorted(range(len(selected)), key=lambda i: stable_seed(config['random_seed'], selected.iloc[i].gene)))
    np.save(ROOT/'modisco_input_order.npy', order)
    (ROOT/'modisco').mkdir()
    patterns, entries = [], []
    for j, region in enumerate(REGIONS):
        state('TF_MoDISco', region=region, complete_seeds=5)
        onehot = np.array(np.load(source/f'{region}.npy', mmap_mode='r')[selected.gene_index.to_numpy(dtype=int)])
        hypothetical = mean[:, :, OFFSETS[j]:OFFSETS[j+1]].transpose(0,2,1).copy()
        # Remove padding/unknown bases and repeated physical copies from motif
        # discovery without altering the saved original DeepLIFT attribution.
        metadata = pd.read_csv(Path(config['model_root'])/'data/genes_repartitioned.tsv', sep='\t')
        for i, row in enumerate(selected.itertuples()):
            a = np.concatenate([np.load(source/f'{r}.npy', mmap_mode='r')[row.gene_index] for r in REGIONS], 0).T
            mapping = physical_map(metadata.iloc[row.gene_index], a)
            for item in mapping.values():
                for index in item['indices'][1:]:
                    if OFFSETS[j] <= index < OFFSETS[j+1]:
                        hypothetical[i, index-OFFSETS[j]] = 0
            hypothetical[i, onehot[i].sum(1) != 1] = 0
        pos, neg = TFMoDISco(one_hot=onehot[order].astype(np.float32), hypothetical_contribs=hypothetical[order],
            max_seqlets_per_metacluster=config['modisco_max_seqlets_per_region'],
            final_min_cluster_size=config['modisco_min_pattern_seqlets'], verbose=True)
        filename = ROOT/'modisco'/f'{region}.h5'
        save_hdf5(str(filename), pos, neg, LENGTHS[j])
        with h5py.File(filename, 'r') as handle:
            for sign in ('pos_patterns', 'neg_patterns'):
                if sign not in handle:
                    continue
                for name in handle[sign]:
                    group = handle[sign][name]
                    pattern_id = f'r{j}_{sign}_{name}'
                    matrix = group['sequence'][:]
                    count = int(group['seqlets/n_seqlets'][0])
                    entries.append(meme_entry(pattern_id, matrix, count))
                    patterns.append(dict(pattern_id=pattern_id, region=region, region_index=j,
                        sign=sign, h5_path=str(filename.relative_to(ROOT)), h5_group=f'{sign}/{name}', seqlets=count))
    write_json(ROOT/'patterns.json', patterns)
    if not patterns:
        write_json(ROOT/'no_annotated_patterns.json', {'reason':'TF-MoDISco found no supported patterns; no forced motif assignment or perturbation'})
        return []
    (ROOT/'patterns.meme').write_text(MEME_HEADER+''.join(entries))
    motifs = json.loads((ROOT/'motifs.json').read_text())
    database = MEME_HEADER+''.join(meme_entry(v['matrix_id'], np.array([v['pfm'][c] for c in 'ACGT']).T, 100) for v in motifs.values())
    (ROOT/'jaspar_125.meme').write_text(database)
    state('TOMTOM', patterns=len(patterns))
    command = [config['tomtom'], '-oc', str(ROOT/'tomtom'), '-min-overlap', str(config['tomtom_min_overlap']),
        '-thresh', str(config['tomtom_q_threshold']), str(ROOT/'patterns.meme'), str(ROOT/'jaspar_125.meme')]
    with (ROOT/'tomtom.log').open('w') as handle:
        subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, check=True)
    table = pd.read_csv(ROOT/'tomtom'/'tomtom.tsv', sep='\t', comment='#')
    id_to_tf = {v['matrix_id']:k for k,v in motifs.items()}
    annotation = []
    for pattern in patterns:
        matches = table.loc[table['Query_ID'].eq(pattern['pattern_id']) & table['q-value'].le(config['tomtom_q_threshold']) & table.Overlap.ge(config['tomtom_min_overlap'])].sort_values(['q-value','p-value','Target_ID'])
        annotation.append({**pattern, 'matches':[{'TF':id_to_tf[r['Target_ID']], 'matrix_id':r['Target_ID'],
            'p_value':float(r['p-value']), 'q_value':float(r['q-value']), 'orientation':r['Orientation']} for _,r in matches.iterrows()]})
    write_json(ROOT/'annotated_patterns.json', annotation)
    return annotation


def pwm_for(motif, pseudocount):
    counts = np.array([motif['pfm'][c] for c in 'ACGT'], dtype=float).T
    ppm = (counts+pseudocount)/(counts+ pseudocount).sum(1, keepdims=True)
    return np.log2(ppm/.25)


def best_footprint(dna, pwm, start, end):
    # This is a PWM score, NOT a binding probability or FIMO p-value.
    choices = []
    for position in range(max(0,start), min(dna.shape[-1]-len(pwm),end-len(pwm))+1):
        local = dna[:,position:position+len(pwm)].T
        if not (local.sum(1) == 1).all():
            continue
        for reverse in (False,True):
            matrix = pwm[::-1,::-1] if reverse else pwm
            choices.append((float((local*matrix).sum()), position, reverse))
    return max(choices, default=None, key=lambda x:(x[0],-x[1],-x[2]))


def freeze_perturbations(config, selected, metadata, annotations):
    import h5py
    motifs = json.loads((ROOT/'motifs.json').read_text())
    order = np.load(ROOT/'modisco_input_order.npy')
    source = Path(config['source_root'])/'data/bpnet_prepared'
    eligible = sorted((x for x in annotations if x['matches']), key=lambda x:(-x['seqlets'],x['pattern_id']))[:config['max_patterns_to_perturb']]
    plan, skipped, seen = [], [], set()
    for pattern in eligible:
        j = pattern['region_index']
        target = pattern['matches'][0]  # All alternate TF matches retained in plan.
        pwm = pwm_for(motifs[target['TF']], config['pwm_pseudocount'])
        candidates = []
        scan_folder=ROOT/'instance_scans'/pattern['pattern_id']
        scan_folder.mkdir(parents=True)
        counts=np.array([motifs[target['TF']]['pfm'][c] for c in 'ACGT']).T
        (scan_folder/'target.meme').write_text(MEME_HEADER+meme_entry(target['matrix_id'],counts,100))
        with h5py.File(ROOT/pattern['h5_path'], 'r') as handle:
            seqlets = handle[pattern['h5_group']+'/seqlets']
            records=[]
            for number,(example,start,end) in enumerate(zip(seqlets['example_idx'][:],seqlets['start'][:],seqlets['end'][:])):
                i=int(order[int(example)])
                dna=np.array(np.load(source/f'{pattern["region"]}.npy',mmap_mode='r')[selected.iloc[i].gene_index])
                seq=''.join('ACGT'[int(x.argmax())] if x.sum()==1 else 'N' for x in dna[int(start):int(end)])
                records.append(f'>s{number}\n{seq}\n')
            (scan_folder/'seqlets.fa').write_text(''.join(records))
            scan_fimo(config,scan_folder/'target.meme',scan_folder/'seqlets.fa',scan_folder/'fimo.tsv')
            fimo=pd.read_csv(scan_folder/'fimo.tsv',sep='\t',comment='#')
            for number,(example,start,end) in enumerate(zip(seqlets['example_idx'][:],seqlets['start'][:],seqlets['end'][:])):
                i = int(order[int(example)])
                row = selected.iloc[i]
                array = np.array(np.load(source/f'{pattern["region"]}.npy', mmap_mode='r')[row.gene_index]).T
                # HDF5 start/end are forward-input coordinates, even when the
                # extracted seqlet sequence is stored reverse complemented.
                hits=fimo.loc[fimo.sequence_name.eq(f's{number}') & fimo['p-value'].le(config['fimo_p_threshold'])].sort_values(['p-value','start','strand'])
                if hits.empty:
                    skipped.append(dict(pattern_id=pattern['pattern_id'],gene=row.gene,reason='No target JASPAR FIMO match at fixed p<=1e-4 inside seqlet'))
                    continue
                hit=hits.iloc[0]
                position=int(start)+int(hit.start)-1
                reverse=hit.strand=='-'
                oriented_pwm=pwm[::-1,::-1] if reverse else pwm
                score=float((array[:,position:position+len(pwm)].T*oriented_pwm).sum())
                # Freeze the illustrative SNP using sequence score alone.
                choices=[]
                for t in range(len(pwm)):
                    ref=int(array[:,position+t].argmax())
                    choices.extend((float(oriented_pwm[t,a]-oriented_pwm[t,ref]),position+t,a,ref)
                        for a in range(4) if a!=ref)
                highlight=min(choices)
                center = position+len(pwm)//2
                window_start = max(0,min(center-config['ism_window_bp']//2,LENGTHS[j]-config['ism_window_bp']))
                window_end = window_start+config['ism_window_bp']
                if not (array[:,window_start:window_end].sum(0) == 1).all():
                    skipped.append(dict(pattern_id=pattern['pattern_id'], gene=row.gene, reason='N/padding in full ISM window'))
                    continue
                meta = metadata.iloc[int(row.gene_index)]
                physical_left = int(meta[f'{pattern["region"]}_requested_start0']+position if meta.strand=='+' else meta[f'{pattern["region"]}_requested_end0']-position-len(pwm))
                key=(meta.chrom, physical_left, physical_left+len(pwm), target['matrix_id'])
                if key in seen:
                    continue
                candidates.append(dict(pattern_id=pattern['pattern_id'], region=pattern['region'], region_index=j,
                    gene=row.gene, gene_row=i, gene_index=int(row.gene_index), chrom=meta.chrom, strand=meta.strand,
                    dna_group=str(row.dna_overlap_group), TF=target['TF'], matrix_id=target['matrix_id'],
                    all_TF_matches=pattern['matches'], pattern_sign=pattern['sign'], motif_start=position,
                    motif_end=position+len(pwm), motif_reverse=bool(reverse), WT_pwm_score_bits=score,
                    WT_fimo_p=float(hit['p-value']),
                    highlight_position=highlight[1], highlight_alternative='ACGT'[highlight[2]],
                    highlight_reference='ACGT'[highlight[3]], highlight_delta_PWM_bits=highlight[0],
                    window_start=window_start, window_end=window_end, physical_motif_start0=physical_left,
                    selection_hash=stable_seed(config['random_seed'],row.gene,pattern['pattern_id'],position)))
        # At most one instance per gene/pattern; order never uses mutation effect.
        candidates.sort(key=lambda x:(x['selection_hash'],x['gene']))
        genes = set()
        retained = 0
        for candidate in candidates:
            if candidate['gene'] in genes:
                continue
            key=(candidate['chrom'],candidate['physical_motif_start0'],candidate['physical_motif_start0']+len(pwm),candidate['matrix_id'])
            if key in seen:
                continue
            candidate['instance_id'] = f'instance_{len(plan):03d}'
            plan.append(candidate)
            seen.add(key)
            genes.add(candidate['gene'])
            retained += 1
            if retained >= config['instances_per_pattern']:
                break
    write_json(ROOT/'frozen_perturbation_plan.json', dict(instances=plan, skipped=skipped,
        selection='up to ten most-supported annotated patterns; five SHA256-ordered distinct genes per pattern; never use ISM effects to select',
        WT_hit_policy='FIMO p<=1e-4 target hit inside a real TF-MoDISco seqlet; lowest p then coordinate. Not measured binding.',
        plan_frozen_before_mutant_inference=True))
    return plan


def scan_fimo(config,meme,fasta,output):
    command=[config['fimo'],'--text','--thresh',str(config['fimo_p_threshold']),str(meme),str(fasta)]
    with Path(output).open('w') as out, Path(str(output)+'.stderr').open('w') as err:
        subprocess.run(command,stdout=out,stderr=err,check=True)
    stderr=Path(str(output)+'.stderr').read_text().lower()
    assert 'truncat' not in stderr and 'max-stored' not in stderr, 'FIMO truncation evidence'
    header=Path(output).read_text().splitlines()
    assert header and 'sequence_name' in header[0], 'Unexpected FIMO text output, not silently treated as no matches'


def perturb(config, selected, metadata, plan):
    motifs = json.loads((ROOT/'motifs.json').read_text())
    (ROOT/'perturbations').mkdir()
    wait_gpu(config)
    results = {x['instance_id']:[] for x in plan}
    for seed in config['seeds']:
        wait_gpu(config)
        model, processor, identity = load_original(config, seed, 'cuda')
        with np.load(Path(identity['folder'])/'preprocessing.npz') as p:
            mean, std = p['target_mean'].item(), p['target_std'].item()
        for k, item in enumerate(plan):
            actual, regulation = gene_inputs(config, item['gene_index'], processor, 'cuda')
            wrapper = ConditionalDNA(model, regulation, mean, std)
            with torch.no_grad():
                wt = float(wrapper(actual).item())
            mapping = physical_map(metadata.iloc[item['gene_index']], actual.cpu().numpy())
            coordinates, variants = [], []
            matrix = np.zeros((4, item['window_end']-item['window_start']), dtype=np.float32)
            for position in range(item['window_start'], item['window_end']):
                global_i = int(OFFSETS[item['region_index']]+position)
                coordinate = next(p for p,x in mapping.items() if global_i in x['indices'])
                ref = mapping[coordinate]['allele']
                for alternative in range(4):
                    if alternative != ref:
                        coordinates.append((coordinate,alternative,position))
            for start in range(0,len(coordinates),config['inference_batch_size']):
                batch = coordinates[start:start+config['inference_batch_size']]
                mutants = torch.cat([edit_physical(actual,mapping,c,a) for c,a,_ in batch])
                with torch.no_grad():
                    prediction = wrapper(mutants).cpu().numpy()
                assert np.isfinite(prediction).all()
                for (coordinate,alternative,position), value in zip(batch,prediction):
                    delta = float(value-wt)
                    matrix[alternative,position-item['window_start']] = delta
                    variants.append(dict(position=position, coordinate0=coordinate, alternative='ACGT'[alternative],
                        delta_signed_Y_PC1=delta, synchronized_copies=len(mapping[coordinate]['indices'])))
            folder = ROOT/'perturbations'/item['instance_id']
            folder.mkdir(exist_ok=True)
            np.save(folder/f'seed_{seed}_delta.npy', matrix)
            write_json(folder/f'seed_{seed}_variants.json', dict(seed=seed, WT_prediction_signed_Y_PC1=wt, variants=variants, checkpoint=identity))
            results[item['instance_id']].append(matrix)
            state('saturation_mutagenesis', seed=seed, instance=k+1, total_instances=len(plan), gene=item['gene'], TF=item['TF'])
        del model, wrapper
        torch.cuda.empty_cache()
    statistics = []
    for item in plan:
        folder = ROOT/'perturbations'/item['instance_id']
        stack = np.stack(results[item['instance_id']])
        np.save(folder/'mean_delta.npy', stack.mean(0))
        np.save(folder/'seed_sd_delta.npy', stack.std(0,ddof=1))
        source=Path(config['source_root'])/'data/bpnet_prepared'
        dna=np.array(np.load(source/f'{item["region"]}.npy', mmap_mode='r')[item['gene_index']]).T
        pwm=pwm_for(motifs[item['TF']], config['pwm_pseudocount'])
        if item['motif_reverse']:
            pwm=pwm[::-1,::-1]
        rows=[]
        for position in range(item['window_start'],item['window_end']):
            ref=int(dna[:,position].argmax())
            motif_pos=position-item['motif_start']
            inside=0 <= motif_pos < len(pwm)
            for alternative in range(4):
                if alternative == ref:
                    continue
                dscore=float(pwm[motif_pos,alternative]-pwm[motif_pos,ref]) if inside else 0.0
                offset=position-item['window_start']
                effects=stack[:,alternative,offset]
                rows.append(dict(instance_id=item['instance_id'], gene=item['gene'], TF=item['TF'], region=item['region'],
                    position=position, reference='ACGT'[ref], alternative='ACGT'[alternative], inside_motif=inside,
                    chrom=item['chrom'], strand=item['strand'],
                    genomic_position1=int(metadata.iloc[item['gene_index']][f'{item["region"]}_requested_start0']+position+1
                        if item['strand']=='+' else metadata.iloc[item['gene_index']][f'{item["region"]}_requested_end0']-position),
                    genomic_reference='ACGT'[ref] if item['strand']=='+' else 'TGCA'[ref],
                    genomic_alternative='ACGT'[alternative] if item['strand']=='+' else 'TGCA'[alternative],
                    delta_PWM_same_site_bits=dscore, delta_signed_Y_PC1=float(effects.mean()),
                    seed_SD=float(effects.std(ddof=1)), positive_seeds=int((effects>0).sum()), negative_seeds=int((effects<0).sum()),
                    motif_annotation='score_decreased' if dscore <= -config['motif_score_change_bits'] else ('score_increased' if dscore >= config['motif_score_change_bits'] else 'small_score_change')))
        table=pd.DataFrame(rows)
        # Re-scan WT and every actual mutant sequence: motif labels must be
        # observed in sequence scoring, not inferred from the model heatmap.
        counts=np.array([motifs[item['TF']]['pfm'][c] for c in 'ACGT']).T
        (folder/'target.meme').write_text(MEME_HEADER+meme_entry(item['matrix_id'],counts,100))
        wtseq=''.join('ACGT'[int(x.argmax())] for x in dna[:,item['window_start']:item['window_end']].T)
        records=[f'>WT\n{wtseq}\n']
        for i,row in table.iterrows():
            position=int(row.position)-item['window_start']
            mutant=wtseq[:position]+row.alternative+wtseq[position+1:]
            records.append(f'>v{i}\n{mutant}\n')
        (folder/'WT_mutants.fa').write_text(''.join(records))
        scan_fimo(config,folder/'target.meme',folder/'WT_mutants.fa',folder/'WT_mutants_fimo.tsv')
        fimo=pd.read_csv(folder/'WT_mutants_fimo.tsv',sep='\t',comment='#')
        start1=item['motif_start']-item['window_start']+1
        end1=item['motif_end']-item['window_start']
        strand='-' if item['motif_reverse'] else '+'
        def hitset(name):
            return {(int(r.start),int(r.stop),r.strand) for r in fimo.loc[fimo.sequence_name.eq(name)].itertuples()}
        wt_hits=hitset('WT')
        assert (start1,end1,strand) in wt_hits, 'Frozen WT FIMO hit not reproduced in plotting window'
        losses,overlaps,gains=[],[],[]
        for i in table.index:
            hits=hitset(f'v{i}')
            losses.append((start1,end1,strand) not in hits)
            overlaps.append(any(s<=end1 and e>=start1 for s,e,_ in hits))
            gains.append(len(hits-wt_hits))
        table['same_site_FIMO_loss']=losses
        table['overlapping_target_match_remaining']=overlaps
        table['new_target_FIMO_hits']=gains
        # Highlight the strongest pre-inference sequence-score decrease, never
        # the mutation with the largest prediction effect.
        highlight=table.loc[table.inside_motif].sort_values(['delta_PWM_same_site_bits','position','alternative']).iloc[0]
        assert (int(highlight.position),highlight.alternative)==(item['highlight_position'],item['highlight_alternative'])
        write_json(folder/'highlight.json',dict(position=int(highlight.position),reference=highlight.reference,
            alternative=highlight.alternative,delta_PWM_bits=float(highlight.delta_PWM_same_site_bits),
            FIMO_same_site_loss=bool(highlight.same_site_FIMO_loss),
            overlapping_target_match_remaining=bool(highlight.overlapping_target_match_remaining),
            rule='Minimum fixed PWM score change; independent of predicted mutation effect'))
        table.to_csv(folder/'variant_effects.tsv',sep='\t',index=False)
        from scipy.stats import pearsonr, spearmanr
        local=table.loc[table.inside_motif]
        r=rho=None
        if len(local)>=3 and local.delta_PWM_same_site_bits.max()>local.delta_PWM_same_site_bits.min() and local.delta_signed_Y_PC1.max()>local.delta_signed_Y_PC1.min():
            r=float(pearsonr(local.delta_PWM_same_site_bits,local.delta_signed_Y_PC1).statistic)
            rho=float(spearmanr(local.delta_PWM_same_site_bits,local.delta_signed_Y_PC1).statistic)
            r=r if np.isfinite(r) else None
            rho=rho if np.isfinite(rho) else None
        statistics.append(dict(instance_id=item['instance_id'], gene=item['gene'], TF=item['TF'],
            pattern_sign=item['pattern_sign'], SNPs_inside_motif=len(local),
            PWM_score_change_vs_prediction_r=r, PWM_score_change_vs_prediction_rho=rho,
            mean_abs_effect_motif=float(local.delta_signed_Y_PC1.abs().mean()),
            mean_abs_effect_flanks=float(table.loc[~table.inside_motif,'delta_signed_Y_PC1'].abs().mean()),
            strict_FIMO_loss_SNPs=int((table.inside_motif & table.same_site_FIMO_loss & ~table.overlapping_target_match_remaining).sum()),
            note='Substitutions are dependent; descriptive within-site comparison, not independent biological tests. Motif score loss does not necessarily lower model Y_PC1.'))
        write_json(folder/'complete.json',dict(status='complete', **item, seeds=config['seeds'],
            outputs_sha256={p.name:digest(p) for p in folder.iterdir() if p.is_file() and p.name!='complete.json'}))
    write_json(ROOT/'perturbation_summary.json',statistics)


def main():
    config=json.loads((ROOT/'experiment_config.json').read_text())
    os.environ['CUDA_VISIBLE_DEVICES']=str(config['physical_gpu'])
    torch.set_num_threads(4)
    np.random.seed(config['random_seed'])
    torch.manual_seed(config['random_seed'])
    torch.backends.cuda.matmul.allow_tf32=False
    torch.backends.cudnn.allow_tf32=False
    lock=(ROOT/'run.lock').open('a')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    selected,metadata,manifest=prepare(config)
    attribute(config,selected)
    annotations=discover(config,selected)
    plan=freeze_perturbations(config,selected,metadata,annotations)
    if plan:
        perturb(config,selected,metadata,plan)
        from plot_perturbations import make_plots
        make_plots(ROOT,config,plan)
    assert checkpoint_hashes(config)==manifest['original_checkpoint_files']
    assert all(digest(name)==expected for name,expected in manifest['original_inputs_sha256'].items())
    for name,expected in manifest['code_sha256'].items():
        assert digest(ROOT/name)==expected, 'Code changed during run'
    write_json(ROOT/'complete.json',dict(status='complete',genes=len(selected),seeds=config['seeds'],
        annotated_patterns=sum(bool(x['matches']) for x in annotations),perturbed_instances=len(plan),
        no_motif_is_not_zero_biological_effect=not bool(plan),
        original_checkpoint_files_unchanged=True, optimizer_updates=0,scope=config['scope'],
        outputs_sha256={str(p.relative_to(ROOT)):digest(p) for p in ROOT.rglob('*') if p.is_file() and not any(x in p.parts for x in ('analysis_venv','__pycache__')) and p!=ROOT/'complete.json' and p.name not in ('controller.log','controller_state.json','run.lock')}))
    write_json(ROOT/'controller_state.json',dict(status='complete',pid=os.getpid(),stage='complete',perturbed_instances=len(plan)))


if __name__=='__main__':
    try:
        main()
    except Exception as e:
        write_json(ROOT/'controller_state.json',dict(status='failed',pid=os.getpid(),error_type=type(e).__name__,error=str(e),traceback=traceback.format_exc()))
        traceback.print_exc()
        raise
