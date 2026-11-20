"""Conditional ChIP IG. No training and no cached FiLM-dependent embeddings."""
import hashlib
import json
from pathlib import Path

import numpy as np
import torch
from torch import nn
from captum.attr import IntegratedGradients

REGIONS = ('tss_sequence', 'first_donor_sequence',
           'first_acceptor_sequence', 'pre_tes_sequence')
LENGTHS = (600, 2000, 2000, 1000)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def save_json(path, obj):
    # Output directories are exclusively created by the runner.
    Path(path).write_text(json.dumps(obj, indent=2, ensure_ascii=False,
                                   allow_nan=False) + '\n', encoding='utf-8')


class Readout(nn.Module):
    def __init__(self, model, target_mean, target_std, absolute):
        super().__init__()
        self.model = model
        self.mean = float(target_mean)
        self.std = float(target_std)
        self.absolute = absolute

    def forward(self, regulatory, *dna):
        require(len(dna) == 4, 'Four actual DNA tensors are required')
        batch = dict(zip(REGIONS, dna))
        batch['regulatory_features'] = regulatory
        delta = self.model(batch).reshape(-1) * self.std
        return delta.abs() if self.absolute else delta + self.mean


def attribute(wrapper, actual, baseline, dna, chip, cfg, target_sd):
    """One gene. All ChIP coordinates move jointly; others have zero delta.

    Captum without the input multiplier gives the path-integrated gradient.
    We explicitly multiply it by actual-baseline and check completeness.
    This also preserves average sensitivity when an input difference is zero.
    """
    require(actual.shape == (1, 2554), 'Expected exactly one 2554-column input')
    other = np.setdiff1d(np.arange(2554), chip)
    require(torch.equal(actual[:, other], baseline[:, other]), 'Other inputs changed')
    with torch.no_grad():
        start = float(wrapper(baseline, *dna).item())
        end = float(wrapper(actual, *dna).item())
    delta = actual - baseline
    previous = None
    history = []
    ig = IntegratedGradients(wrapper, multiply_by_inputs=False)
    atol = cfg['absolute_tolerance_in_target_sd'] * target_sd
    for steps in cfg['ig_steps']:
        mean_gradient = ig.attribute(
            actual, baselines=baseline, additional_forward_args=dna,
            n_steps=int(steps), method='gausslegendre',
            internal_batch_size=int(cfg['internal_batch_size']))
        values = (mean_gradient * delta).detach().cpu().numpy()[0]
        require(np.isfinite(values).all(), 'Nonfinite attribution')
        require(np.all(values[other] == 0), 'Non-ChIP attribution must be exactly zero')
        residual = float(values.astype(np.float64).sum() - (end-start))
        threshold = atol + cfg['relative_tolerance'] * abs(end-start)
        change = None if previous is None else float(np.abs(values-previous).sum())
        stability_limit = atol + cfg['step_stability_relative_tolerance'] * float(np.abs(values).sum())
        accepted = abs(residual) <= threshold and change is not None and change <= stability_limit
        history.append(dict(steps=int(steps), residual=residual,
                            completeness_limit=threshold, attribution_L1_change=change,
                            stability_limit=stability_limit, accepted=bool(accepted)))
        if accepted:
            return values[chip], mean_gradient.detach().cpu().numpy()[0, chip], {
                'passed': True, 'baseline_readout': start, 'actual_readout': end,
                'endpoint_difference': end-start, 'history': history}
        previous = values
    return None, None, {'passed': False, 'baseline_readout': start,
                       'actual_readout': end, 'endpoint_difference': end-start,
                       'history': history}


def load_model(cfg, fold, seed, device):
    # Imports resolve to the original, subsequently hash-checked implementation.
    from region_film.model import PolIIFiLMRegressor
    from pausing_utils.preprocessing import FoldPreprocessor
    root = Path(cfg['model_root'])
    training = read_json(root/'training_config.json')
    directory = root/'result'/training['run_name']/f'fold_{fold:02d}'/f'seed_{seed}'
    files = {name: sha(directory/name) for name in
             ('best.pt', 'preprocessing.npz', 'run_spec.json', 'complete.json', 'test_predictions.tsv')}
    spec = read_json(directory/'run_spec.json')
    done = read_json(directory/'complete.json')
    ckpt = torch.load(directory/'best.pt', map_location='cpu', weights_only=True)
    require(spec['fold'] == fold and spec['seed'] == seed, 'Wrong fold/seed')
    require(ckpt['fingerprint'] == spec['fingerprint'] == done['fingerprint'], 'Checkpoint mismatch')
    for key, origin in [('implementation_sha256', root),
                        ('reused_source_sha256', Path(cfg['source_root']))]:
        for name, expected in spec[key].items():
            require(sha(origin/name) == expected, f'Implementation changed: {origin/name}')
    require(ckpt['architecture_config']['region_lengths'] == dict(zip(REGIONS, LENGTHS)),
            'Not the original four-window A model')
    model = PolIIFiLMRegressor(ckpt['architecture_config'], regulatory_dim=2554)
    model.load_state_dict(ckpt['model'], strict=True)
    model.eval().to(device).requires_grad_(False)
    require(hasattr(model, 'film_mlp_conv') and not hasattr(model, 'regulatory_mlp'), 'Wrong model family')
    norm = model.input_normalization_summary()
    require(norm['fitted'] and not norm['batch_statistics_used'], 'BN is not fixed')
    require(norm['statistics_sha256'] == done['input_normalization']['statistics_sha256'], 'BN hash mismatch')
    processor = FoldPreprocessor.load(directory/'preprocessing.npz')
    require(float(processor.target_std) > 0, 'Invalid target scale')
    return model, processor, directory, done, {
        'fold': fold, 'seed': seed, 'files_sha256': files, 'normalization': norm}
