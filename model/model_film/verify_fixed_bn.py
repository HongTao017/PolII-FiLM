"""Regression tests for frozen, fold-fitted FiLM input normalization."""
import io
import json
from pathlib import Path

import torch
from torch import nn
from region_film.modules import FixedTrainingBatchNorm1d


def tests():
    torch.set_num_threads(2)
    torch.manual_seed(20260918)
    bn = FixedTrainingBatchNorm1d(4)
    try:
        bn(torch.zeros(2, 4))
    except RuntimeError:
        pass
    else:
        raise AssertionError('Unfitted normalization accepted inference')
    train = torch.randn(100, 4)
    train[:, 0] = 0
    train[:, 1] = 0
    train[0, 1] = 10  # A rare feature: most minibatches have zero variance.
    bn.fit_training_statistics(train)
    expected_mean = train.double().mean(0).float()
    expected_var = train.double().var(0, unbiased=True).float()
    expected_var[0] = 1
    torch.testing.assert_close(bn.running_mean, expected_mean)
    torch.testing.assert_close(bn.running_var, expected_var)
    assert bn.constant_training_features.tolist() == [True, False, False, False]
    assert int(bn.fit_rows) == 100
    fingerprint = bn.statistics_fingerprint()
    optimizer = torch.optim.AdamW(bn.parameters(), lr=1e-3)
    before_weight = bn.weight.detach().clone()
    before_bias = bn.bias.detach().clone()
    for _ in range(6):
        bn.train()
        optimizer.zero_grad()
        values = torch.ones(8, 4)
        output = bn(values)
        output.square().mean().backward()
        assert bn.weight.grad is not None and bn.bias.grad is not None
        assert torch.isfinite(bn.weight.grad).all() and torch.isfinite(bn.bias.grad).all()
        optimizer.step()
        assert bn.statistics_fingerprint() == fingerprint
        assert int(bn.num_batches_tracked) == 0
    assert not torch.equal(before_weight, bn.weight)
    assert not torch.equal(before_bias, bn.bias)
    # Identical affine normalization in training and evaluation modes.
    probe = torch.tensor([[1., 10., 2., -3.], [2., 0., -1., 1.]])
    bn.train()
    training_output = bn(probe)
    bn.eval()
    torch.testing.assert_close(bn(probe), training_output, rtol=0, atol=0)
    torch.testing.assert_close(bn(probe[:1]), training_output[:1], rtol=0, atol=0)
    # A constant training column that is nonzero at validation is not amplified.
    assert abs(float(bn(probe)[0, 0])) < 2
    try:
        bn.fit_training_statistics(torch.randn(10, 4))
    except RuntimeError:
        pass
    else:
        raise AssertionError('Refitting on validation/other rows was accepted')
    storage = io.BytesIO()
    torch.save(bn.state_dict(), storage)
    storage.seek(0)
    restored = FixedTrainingBatchNorm1d(4)
    restored.load_state_dict(torch.load(storage, weights_only=True))
    assert restored.statistics_fingerprint() == fingerprint
    torch.testing.assert_close(restored(probe), bn(probe), rtol=0, atol=0)
    # Match ordinary eval-mode BN when using the same safe training statistics.
    legacy = nn.BatchNorm1d(4).eval()
    with torch.no_grad():
        legacy.weight.copy_(bn.weight)
        legacy.bias.copy_(bn.bias)
        legacy.running_mean.copy_(bn.running_mean)
        legacy.running_var.copy_(bn.running_var)
    torch.testing.assert_close(legacy(probe), bn(probe), rtol=0, atol=0)
    # Parent.train() propagates its flag but must not re-enable batch statistics.
    parent = nn.Sequential(restored)
    parent.train()
    parent(torch.zeros(64, 4))
    assert restored.statistics_fingerprint() == fingerprint
    assert int(restored.num_batches_tracked) == 0
    return {'status': 'PASS', 'train_eval_same_statistics': True,
            'minibatch_composition_independent': True, 'affine_weights_learn': True,
            'parent_train_does_not_update_statistics': True,
            'constant_columns_safe_variance': 1.0, 'checkpoint_roundtrip': True,
            'unfitted_forward_rejected': True, 'refitting_rejected': True}


if __name__ == '__main__':
    report = tests()
    (Path(__file__).resolve().parent / 'fixed_bn_tests.json').write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))
