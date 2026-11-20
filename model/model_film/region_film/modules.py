# Copyright 2026 Ekin Deniz Aksu (adapted FiLM routines).
# SPDX-License-Identifier: Apache-2.0

# No attention/Transformer definitions, imports or dependencies.
import torch
from torch import nn
from torch.nn import functional as F
import hashlib


class FixedTrainingBatchNorm1d(nn.BatchNorm1d):
    """BN affine transform with statistics fitted once on the training fold.

    Parent model.train() never enables minibatch statistics. No validation or
    test observations are used to fit the buffers. Training-constant features
    receive variance=1 instead of sqrt(eps) amplification; no input is removed.
    """

    def __init__(self, num_features, eps=1e-5):
        super().__init__(num_features, eps=eps, affine=True, track_running_stats=True)
        self.register_buffer('statistics_fitted', torch.tensor(False))
        self.register_buffer('fit_rows', torch.tensor(0, dtype=torch.int64))
        self.register_buffer('constant_training_features', torch.zeros(num_features, dtype=torch.bool))

    @torch.no_grad()
    def fit_training_statistics(self, training_features):
        if bool(self.statistics_fitted):
            raise RuntimeError('Fixed BN statistics already fitted; do not refit during training')
        values = torch.as_tensor(training_features, dtype=torch.float64, device='cpu')
        if values.ndim != 2 or values.shape[1] != self.num_features or values.shape[0] < 2:
            raise ValueError('Fit requires complete training rows with the configured feature count')
        if not torch.isfinite(values).all():
            raise ValueError('Fixed BN requires finite preprocessed training features')
        mean = values.mean(dim=0)
        variance = values.var(dim=0, unbiased=True)
        constant = variance <= 1e-8
        safe_variance = torch.where(constant, torch.ones_like(variance), variance)
        self.running_mean.copy_(mean.to(self.running_mean))
        self.running_var.copy_(safe_variance.to(self.running_var))
        self.constant_training_features.copy_(constant.to(self.constant_training_features.device))
        self.num_batches_tracked.zero_()
        self.fit_rows.fill_(values.shape[0])
        self.statistics_fitted.fill_(True)

    def statistics_fingerprint(self):
        if not bool(self.statistics_fitted):
            raise RuntimeError('Fixed BN statistics have not been fitted')
        digest = hashlib.sha256()
        for name in ('running_mean', 'running_var', 'constant_training_features', 'fit_rows'):
            digest.update(getattr(self, name).detach().cpu().contiguous().numpy().tobytes())
        digest.update(str(self.eps).encode())
        return digest.hexdigest()

    def statistics_summary(self):
        return {'kind': 'fixed_training_batchnorm', 'fitted': bool(self.statistics_fitted),
                'fit_rows': int(self.fit_rows), 'constant_features': int(self.constant_training_features.sum()),
                'constant_variance': 1.0, 'eps': self.eps, 'affine_trainable': True,
                'batch_statistics_used': False, 'running_statistics_updated': False,
                'variance_min': float(self.running_var.min()),
                'variance_max': float(self.running_var.max()),
                'statistics_sha256': self.statistics_fingerprint()}

    def forward(self, x):
        if not bool(self.statistics_fitted):
            raise RuntimeError('Fit input BN on the training fold before forward or evaluation')
        self._check_input_dim(x)
        # Deliberately independent of self.training: affine weights still learn.
        return F.batch_norm(x, self.running_mean, self.running_var, self.weight,
                            self.bias, training=False, momentum=0.0, eps=self.eps)


class FiLM(nn.Module):
    def __init__(self, num_channels):
        super().__init__()
        self.num_channels = num_channels

    def forward(self, x, scale, shift):
        return x * scale.unsqueeze(-1) + shift.unsqueeze(-1)


class FiLM_MLP(nn.Module):
    def __init__(self, film_dim, input_dim, hidden_dims, dropout):
        super().__init__()
        self.film_dimensions = film_dim
        self.input_dim = input_dim
        self.hidden_dims = hidden_dims
        self.dropout = dropout
        self.output_dim = 2 * sum(film_dim)
        layers = [FixedTrainingBatchNorm1d(input_dim)]
        prev_dim = input_dim
        for hd in hidden_dims:
            layers.extend([nn.Linear(prev_dim, hd), nn.GELU(approximate="tanh"),
                           nn.Dropout(dropout)])
            prev_dim = hd
        self.film_network = nn.Sequential(*layers)
        self.final_layer = nn.Linear(prev_dim, self.output_dim)

    def forward(self, x):
        return torch.tanh(self.final_layer(self.film_network(x)) / 3) * 3.0


class TargetLengthCrop(nn.Module):
    def __init__(self, target_length):
        super().__init__()
        self.target_length = target_length

    def forward(self, x):
        if self.target_length == -1:
            return x
        if x.shape[-1] < self.target_length:
            raise ValueError("Input sequence is shorter than the requested crop")
        offset = (x.shape[-1] - self.target_length) // 2
        return x[..., offset:offset + self.target_length]
