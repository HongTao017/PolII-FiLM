"""Synthetic scientific checks, followed separately by real-checkpoint preflight."""
import sys
sys.dont_write_bytecode=True
import numpy as np
import torch
from torch import nn
from core import Readout, attribute
from crosscell import build_reference


def main():
    groups=[{'indices':[0,1],'features':['A_gene_body_count','A_pause_count'],'category':'TF','tf':'A','cell':'c1'},
            {'indices':[2,3],'features':['A_gene_body_count','A_pause_count'],'category':'TF','tf':'A','cell':'c2'},
            {'indices':[4,5],'features':['B_gene_body_count','B_pause_count'],'category':'TF','tf':'B','cell':'c1'}]
    raw=np.zeros((2,2554));raw[0,:7]=[2,4,6,8,9,10,27];raw[1,:7]=[3,4,3,4,5,6,28]
    reference,z,eligible,varying,_=build_reference(raw,groups)
    assert np.allclose(reference[0,:7],[4,6,4,6,9,10,27])
    assert np.allclose(z[0,:2],[-1,1]) and np.isnan(z[:,2]).all()
    assert np.isnan(z[1,:2]).all() and eligible.tolist()==[True,True,False]
    assert np.array_equal(reference[:,4:],raw[:,4:])
    class Toy(nn.Module):
        def forward(self,batch):
            x=batch['regulatory_features']
            return x[:,0]*x[:,2]+0.2*x[:,1]**2+0.01*x[:,6]+batch['tss_sequence'].sum((1,2))*0.001
    model=Toy().eval()
    actual=torch.tensor(raw[:1],dtype=torch.float32)
    baseline=torch.tensor(reference[:1],dtype=torch.float32)
    dna=tuple(torch.ones((1,3,4)) for _ in range(4))
    cfg={'absolute_tolerance_in_target_sd':1e-4,'relative_tolerance':.01,
         'step_stability_relative_tolerance':.02,'ig_steps':[16,32,64], 'internal_batch_size':4}
    for absolute in (True,False):
        wrapper=Readout(model,10,2,absolute)
        a,gradient,qc=attribute(wrapper,actual,baseline,dna,np.arange(6),cfg,2.)
        assert qc['passed'] and np.isfinite(gradient).all()
        assert abs(a.sum()-qc['endpoint_difference'])<1e-4
        assert a[4]==a[5]==0
    print('PASS: same-target reference, singleton handling, fixed other features, z-score, IG completeness')


if __name__=='__main__':main()
