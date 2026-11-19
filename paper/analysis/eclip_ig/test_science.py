"""Synthetic scientific checks, followed separately by real-checkpoint preflight."""
import sys
sys.dont_write_bytecode=True
import numpy as np
import torch
from torch import nn
from core import Readout
from integrate import integrate


def main():
    raw=np.zeros((1,2554));raw[0,:7]=[2,4,6,8,0,0,27]
    reference=raw.copy();reference[:,:6]=0
    assert np.array_equal(reference[:,6:],raw[:,6:])
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
        a,gradient,qc=integrate(wrapper,actual,baseline,dna,np.arange(6),cfg,2.)
        assert qc['passed'] and np.isfinite(gradient).all()
        assert abs(a.sum()-qc['endpoint_difference'])<1e-4
        assert a[4]==a[5]==0
    print('PASS: eCLIP zero-count reference, fixed other features, and IG completeness')


if __name__=='__main__':main()
