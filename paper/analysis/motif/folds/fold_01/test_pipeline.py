"""Deterministic unit tests for physical mapping/editing and DeepLIFT rules."""
import unittest
import numpy as np
import pandas as pd
import torch
from deep_lift_core import (LENGTHS, OFFSETS, REGIONS, GroupNormRescale, Rescale,
    dinuc_counts, edit_physical, physical_map, shuffled_references, digest, write_json)
from preflight import unit_tests


def example(strand='+'):
    meta=dict(gene='synthetic_not_biological',strand=strand)
    # Deliberately overlap TSS/donor/acceptor; last window is disjoint.
    for region,start,length in zip(REGIONS,(1000,1200,2000,10000),LENGTHS):
        meta[f'{region}_requested_start0']=start
        meta[f'{region}_requested_end0']=start+length
    total=int(OFFSETS[-1])
    array=np.zeros((4,total),dtype=np.int8)
    for j,region in enumerate(REGIONS):
        left,right=meta[f'{region}_requested_start0'],meta[f'{region}_requested_end0']
        for i in range(LENGTHS[j]):
            coordinate=left+i if strand=='+' else right-1-i
            # Same oriented allele for all copies; N at fixed physical positions.
            if coordinate%97 != 0:
                array[(coordinate*7+coordinate//11)%4,OFFSETS[j]+i]=1
    return pd.Series(meta),array


class PipelineTests(unittest.TestCase):
    def test_groupnorm_rescale_conservation(self):
        self.assertEqual(len(unit_tests()),4)

    def test_relu_gelu(self):
        for kind in ('relu','gelu'):
            x=torch.tensor([[-2.,-.1,0.,.1,2.],[1.,.1,0.,-.1,-1.]],requires_grad=True)
            y=Rescale.apply(x,kind)
            g=torch.autograd.grad(y[0].sum(),x)[0][0]
            self.assertTrue(torch.allclose((g*(x[0]-x[1])).sum(),(y[0]-y[1]).sum(),atol=1e-6))

    def test_physical_mapping_and_editing_both_strands(self):
        for strand in ('+','-'):
            meta,array=example(strand)
            mapping=physical_map(meta,array)
            coordinate=next(p for p,x in mapping.items() if len(x['indices'])>=2)
            original=torch.from_numpy(array[None].astype(np.float32))
            alt=(mapping[coordinate]['allele']+1)%4
            edited=edit_physical(original,mapping,coordinate,alt)
            changed=(edited!=original).any(1)[0].nonzero().flatten().tolist()
            self.assertEqual(sorted(changed),sorted(mapping[coordinate]['indices']))
            physical_map(meta,edited.numpy())

    def test_overlapping_disagreement_fails(self):
        meta,array=example()
        mapping=physical_map(meta,array)
        item=next(x for x in mapping.values() if len(x['indices'])>=2)
        bad=array.copy()
        i=item['indices'][1]
        bad[:,i]=0
        bad[(item['allele']+1)%4,i]=1
        with self.assertRaises(AssertionError):
            physical_map(meta,bad)

    def test_reference_reproducibility_mask_and_overlap(self):
        for strand in ('+','-'):
            meta,array=example(strand)
            first=shuffled_references(meta,array,2,20261015)
            second=shuffled_references(meta,array,2,20261015)
            self.assertTrue(np.array_equal(first,second))
            self.assertTrue(np.array_equal(first.sum(1),np.repeat(array.sum(0)[None],2,0)))
            for ref in first:
                physical_map(meta,ref)


if __name__=='__main__':
    from pathlib import Path
    torch.set_num_threads(4)
    result=unittest.main(verbosity=2,exit=False).result
    root=Path(__file__).resolve().parent
    write_json(root/'tests_complete.json',dict(status='PASS' if result.wasSuccessful() else 'FAIL',
        tests=result.testsRun,failures=len(result.failures),errors=len(result.errors),
        code_sha256={p.name:digest(p) for p in root.glob('*.py')}))
    raise SystemExit(0 if result.wasSuccessful() else 1)
