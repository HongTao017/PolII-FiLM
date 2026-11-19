"""Same IG path and tolerances, composite quadrature fallback for kinks."""
import numpy as np
import torch
from core import attribute, require


def integrate(wrapper, actual, baseline, dna, chip, cfg, target_sd):
    result, gradient, report = attribute(wrapper, actual, baseline, dna, chip, cfg, target_sd)
    if report['passed']:
        report['method'] = 'Captum_GaussLegendre'
        return result, gradient, report
    direction = (actual-baseline).detach()
    nodes, weights = np.polynomial.legendre.leggauss(16)
    previous = None
    previous_pass = False
    atol = cfg['absolute_tolerance_in_target_sd']*target_sd
    effect = report['endpoint_difference']
    other = np.setdiff1d(np.arange(2554), chip)
    for intervals in cfg['composite_intervals']:
        alpha = ((np.arange(intervals)[:,None]+(nodes[None,:]+1)/2)/intervals).ravel()
        weight = np.tile(weights/(2*intervals), intervals)
        total = torch.zeros_like(actual, dtype=torch.float64)
        size = cfg['internal_batch_size']
        for start in range(0,len(alpha),size):
            end = min(start+size,len(alpha))
            a = torch.tensor(alpha[start:end],device=actual.device,dtype=actual.dtype)[:,None]
            points = (baseline.detach()+a*direction).requires_grad_(True)
            fixed = tuple(d.expand(end-start,*d.shape[1:]) for d in dna)
            output = wrapper(points,*fixed)
            grad = torch.autograd.grad(output.sum(),points)[0]
            require(bool(torch.isfinite(grad).all()), 'Nonfinite gradient in composite quadrature')
            w = torch.tensor(weight[start:end],device=actual.device,dtype=torch.float64)[:,None]
            total += (w*grad.to(torch.float64)).sum(0,keepdim=True)
        vector = (total*direction).detach().cpu().numpy()[0]
        require(np.isfinite(vector).all() and np.all(vector[other]==0), 'Invalid conditional IG')
        residual = float(vector.sum()-effect)
        limit = atol+cfg['relative_tolerance']*abs(effect)
        change = None if previous is None else float(np.abs(vector-previous).sum())
        change_limit = atol+cfg['step_stability_relative_tolerance']*float(np.abs(vector).sum())
        good = abs(residual)<=limit
        accepted = good and previous_pass and change is not None and change<=change_limit
        report['history'].append({'method':'composite_GL16','steps':int(16*intervals),
                                  'residual':residual,'completeness_limit':limit,
                                  'attribution_L1_change':change,'stability_limit':change_limit,
                                  'accepted':bool(accepted)})
        if accepted:
            report.update(passed=True,method='composite_GL16')
            return vector[chip],total.detach().cpu().numpy()[0,chip],report
        previous, previous_pass = vector,good
    return None,None,report
