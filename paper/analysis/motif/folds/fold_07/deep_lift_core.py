"""Conditional-DNA DeepLIFT: Rescale nonlinearities + symmetric product GN.

This is an explicit, custom rule implementation, NOT vanilla Captum DeepLift.
Forward uses the ORIGINAL model's parameters and original GroupNorm operation.
GroupNorm backward decomposes center -> square -> mean -> rsqrt -> multiply.
Square and rsqrt use Rescale; the product uses symmetric (two-order) splitting.
Completeness alone is necessary, not proof of biological attribution validity.
"""
import hashlib
import json
from pathlib import Path

import numpy as np
import torch
from torch import nn
from torch.nn import functional as F

REGIONS = ('tss_sequence', 'first_donor_sequence', 'first_acceptor_sequence', 'pre_tes_sequence')
LENGTHS = (600, 2000, 2000, 1000)
OFFSETS = np.cumsum((0,) + LENGTHS)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False))
    temporary.replace(path)


def stable_seed(*parts):
    return int.from_bytes(hashlib.sha256('|'.join(map(str, parts)).encode()).digest()[:4], 'little')


class Rescale(torch.autograd.Function):
    @staticmethod
    def forward(ctx, x, kind):
        y = F.relu(x) if kind == 'relu' else F.gelu(x, approximate='tanh')
        ctx.save_for_backward(x, y)
        ctx.kind = kind
        return y

    @staticmethod
    def backward(ctx, grad):
        x, y = ctx.saved_tensors
        a, b = x.chunk(2)
        ya, yb = y.chunk(2)
        dx = a - b
        # Zero delta contributes zero, but a finite local derivative is needed.
        if ctx.kind == 'relu':
            local = ((a + b) > 0).to(x.dtype)
        else:
            middle = (a + b) / 2
            z = 0.7978845608028654 * (middle + .044715 * middle**3)
            t = torch.tanh(z)
            local = .5 * (1 + t) + .5 * middle * (1 - t*t) * .7978845608028654 * (1 + .134145 * middle**2)
        safe = torch.where(dx == 0, torch.ones_like(dx), dx)
        slope = torch.where(dx == 0, local, (ya-yb)/safe)
        return grad * torch.cat((slope, slope)), None


class GroupNormRescale(torch.autograd.Function):
    @staticmethod
    def forward(ctx, x, weight, bias, groups, eps):
        ctx.save_for_backward(x, weight)
        ctx.groups, ctx.eps = groups, eps
        return F.group_norm(x, groups, weight, bias, eps)

    @staticmethod
    def backward(ctx, grad):
        x, weight = ctx.saved_tensors
        n, channels, length = x.shape
        a, b = x.reshape(n, ctx.groups, -1).chunk(2)
        ca = a - a.mean(-1, keepdim=True)
        cb = b - b.mean(-1, keepdim=True)
        va = ca.square().mean(-1, keepdim=True) + ctx.eps
        vb = cb.square().mean(-1, keepdim=True) + ctx.eps
        sa, sb = va.rsqrt(), vb.rsqrt()
        # Stable exact divided difference of reciprocal square root.
        ra, rb = va.sqrt(), vb.sqrt()
        secant_rsqrt = -1 / (ra * rb * (ra + rb))
        direct = (sa + sb) / 2
        factor = (ca + cb) / 2
        g = (grad * weight[None, :, None]).reshape(n, ctx.groups, -1)
        direct = torch.cat((direct, direct))
        factor = torch.cat((factor, factor))
        secant_rsqrt = torch.cat((secant_rsqrt, secant_rsqrt))
        square_slope = torch.cat((ca+cb, ca+cb))
        centered_grad = g * direct + (g*factor).sum(-1, keepdim=True) * secant_rsqrt * square_slope / ca.shape[-1]
        result = centered_grad - centered_grad.mean(-1, keepdim=True)
        return result.reshape_as(x), None, None, None, None


class ConditionalDNA(nn.Module):
    """Original forward, fixed actual regulation; optional paired DeepLIFT rules.

    [B,4,5600] is a storage layout only: each region is split before encoding.
    No convolution, seqlet or motif ever bridges an artificial region boundary.
    """
    def __init__(self, model, regulation, target_mean, target_std):
        super().__init__()
        self.model = model
        self.register_buffer('regulation', regulation.detach().clone())
        self.target_mean, self.target_std = float(target_mean), float(target_std)
        with torch.no_grad():
            params = model.film_mlp_conv(self.regulation)
            scales, shifts = model._split_film_params(params, model.config['film_dimensions_conv'])
        self.scales = {k: v.detach() for k, v in scales.items()}
        self.shifts = {k: v.detach() for k, v in shifts.items()}

    def forward(self, dna, paired=False):
        if paired and dna.shape[0] % 2:
            raise ValueError('DeepLIFT actual/reference batch must be paired')
        vectors = []
        for j, name in enumerate(REGIONS):
            encoder = self.model.region_encoders[name]
            x = encoder.motif_conv(dna[:, :, OFFSETS[j]:OFFSETS[j+1]])
            x = Rescale.apply(x, 'relu') if paired else F.relu(x)
            for i, block in enumerate(encoder.blocks):
                residual = x
                gn = block.normalization
                x = GroupNormRescale.apply(x, gn.weight, gn.bias, gn.num_groups, gn.eps) if paired else gn(x)
                site = self.model.region_film_offsets[name] + i
                x = block.film(x, self.scales[site].expand(len(dna), -1), self.shifts[site].expand(len(dna), -1))
                x = block.conv(x)
                x = Rescale.apply(x, 'relu') if paired else F.relu(x)
                x = residual + x
            x = self.model.crop(self.model.channel_adapters[name](x)).mean(-1)
            vectors.append(x)
        head = self.model.fusion_head
        x = head[0](torch.cat(vectors, 1))
        x = Rescale.apply(x, 'gelu') if paired else head[1](x)
        x = head[3](head[2](x)).squeeze(-1)
        return x * self.target_std + self.target_mean


def explain_pair(wrapper, actual, reference):
    x = torch.cat((actual, reference)).detach().requires_grad_(True)
    output = wrapper(x, paired=True)
    multipliers = torch.autograd.grad(output[:len(actual)].sum(), x)[0][:len(actual)]
    observed = multipliers * (actual-reference)
    # One-hot hypothetical contributions, not the effect of an actual mutation.
    hypothetical = multipliers - (multipliers*reference).sum(1, keepdim=True)
    projected = hypothetical * actual
    delta = output[:len(actual)] - output[len(actual):]
    error = observed.flatten(1).sum(1) - delta
    return hypothetical.detach(), projected.detach(), delta.detach(), error.detach()


def load_original(config, seed, device='cpu'):
    import sys
    for root in (config['source_root'], config['model_root']):
        if root not in sys.path:
            sys.path.insert(0, root)
    from region_film.model import PolIIFiLMRegressor
    from pausing_utils.preprocessing import FoldPreprocessor
    root = Path(config['model_root'])
    train = json.loads((root/'training_config.json').read_text())
    folder = root/'result'/train['run_name']/f'fold_{config["fold"]:02d}'/f'seed_{seed}'
    spec = json.loads((folder/'run_spec.json').read_text())
    complete = json.loads((folder/'complete.json').read_text())
    ckpt = torch.load(folder/'best.pt', map_location='cpu', weights_only=True)
    assert ckpt['fingerprint'] == spec['fingerprint'] == complete['fingerprint']
    assert spec['fold'] == config['fold'] and spec['seed'] == seed
    for key, base in (('implementation_sha256', config['model_root']), ('reused_source_sha256', config['source_root'])):
        for name, expected in spec[key].items():
            assert digest(Path(base)/name) == expected, (name, 'source fingerprint mismatch')
    model = PolIIFiLMRegressor(ckpt['architecture_config'], regulatory_dim=2554)
    model.load_state_dict(ckpt['model'], strict=True)
    model.eval().requires_grad_(False).to(device)
    bn = model.input_normalization_summary()
    assert bn['statistics_sha256'] == complete['input_normalization']['statistics_sha256']
    assert bn['fitted'] and not bn['batch_statistics_used']
    assert tuple(model.config['region_lengths'][r] for r in REGIONS) == LENGTHS
    processor = FoldPreprocessor.load(folder/'preprocessing.npz')
    return model, processor, {'folder': str(folder), 'checkpoint_sha256': digest(folder/'best.pt'),
        'preprocessing_sha256': digest(folder/'preprocessing.npz'), 'fingerprint': ckpt['fingerprint']}


def gene_inputs(config, gene_index, processor, device='cpu'):
    root = Path(config['source_root'])/'data/sequence_prepared'
    pieces = [np.array(np.load(root/f'{r}.npy', mmap_mode='r')[gene_index], dtype=np.float32) for r in REGIONS]
    raw = np.array(np.load(root/'regulatory_log10_unscaled.npy', mmap_mode='r')[[gene_index]], dtype=np.float32)
    features = processor.transform(raw)
    dna = torch.from_numpy(np.concatenate(pieces).T[None]).to(device)
    regulation = torch.from_numpy(np.asarray(features, dtype=np.float32)).to(device)
    return dna, regulation


def physical_map(meta, dna):
    """Map one physical coordinate to all oriented window copies (N excluded)."""
    array = np.asarray(dna).reshape(4, -1)
    mapping = {}
    for j, region in enumerate(REGIONS):
        left, right = meta[f'{region}_requested_start0'], meta[f'{region}_requested_end0']
        if not np.isfinite(left):
            assert array[:, OFFSETS[j]:OFFSETS[j+1]].sum() == 0
            continue
        for i in range(LENGTHS[j]):
            column = array[:, OFFSETS[j]+i]
            if column.sum() == 0:
                continue
            assert column.sum() == 1 and np.isin(column, (0, 1)).all()
            coordinate = int(left+i if meta.strand == '+' else right-1-i)
            allele = int(column.argmax())
            if coordinate in mapping:
                assert mapping[coordinate]['allele'] == allele, 'Overlapping DNA copies disagree'
            else:
                mapping[coordinate] = {'allele': allele, 'indices': []}
            mapping[coordinate]['indices'].append(int(OFFSETS[j]+i))
    return mapping


def dinuc_counts(array):
    letters = array.argmax(0)
    return np.bincount(letters[:-1]*4+letters[1:], minlength=16)


def shuffled_references(meta, actual, count, random_seed):
    from tangermeme.ersatz import dinucleotide_shuffle
    original = np.asarray(actual).reshape(4, -1)
    mapping = physical_map(meta, original)
    coordinates = sorted(mapping, reverse=meta.strand == '-')
    blocks = []
    for coordinate in coordinates:
        if not blocks or abs(coordinate-blocks[-1][-1]) != 1:
            blocks.append([])
        blocks[-1].append(coordinate)
    references = np.repeat(original[None], count, 0).astype(np.int8)
    for block in blocks:
        template = np.eye(4, dtype=np.int8)[[mapping[p]['allele'] for p in block]].T
        if len(block) >= 3:
            shuffled = dinucleotide_shuffle(torch.from_numpy(template[None]).float(), n=count,
                random_state=stable_seed(random_seed, meta.gene, block[0], block[-1])).numpy()[0]
        else:
            shuffled = np.repeat(template[None], count, 0)
        for k in range(count):
            assert np.array_equal(dinuc_counts(template), dinuc_counts(shuffled[k]))
            assert np.array_equal(template.sum(1), shuffled[k].sum(1))
            for t, coordinate in enumerate(block):
                for index in mapping[coordinate]['indices']:
                    references[k, :, index] = shuffled[k, :, t]
    assert np.array_equal(references.sum(1), np.repeat(original.sum(0)[None], count, 0))
    for reference in references:
        physical_map(meta, reference)  # Verify all overlapping copies again.
    return references


def edit_physical(actual, mapping, coordinate, alternative):
    result = actual.clone()
    indices = mapping[int(coordinate)]['indices']
    result[:, :, indices] = 0
    result[:, alternative, indices] = 1
    # No changes at any unrelated position.
    mask = torch.ones(actual.shape[-1], dtype=torch.bool, device=actual.device)
    mask[indices] = False
    assert torch.equal(result[:, :, mask], actual[:, :, mask])
    return result
