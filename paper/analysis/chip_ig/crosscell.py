"""Cross-source-cell reference in log10(count+1) space, before train scaling."""
from collections import defaultdict
import numpy as np


def build_reference(raw_rows, groups):
    raw_rows = np.asarray(raw_rows, dtype=np.float64)
    reference = raw_rows.copy()
    pairs = np.asarray([g['indices'] for g in groups], dtype=int)
    values = raw_rows[:, pairs]
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError('ChIP inputs must be finite, nonnegative log counts')
    families = defaultdict(list)
    for k, group in enumerate(groups):
        if len(group['features']) != 2 or not ('gene_body' in group['features'][0] and 'pause' in group['features'][1]):
            raise ValueError('Expected gene-body then pause columns')
        families[(group['category'], group['tf'])].append(k)
    z = np.full_like(values, np.nan)
    eligible = np.zeros(len(groups), dtype=bool)
    varying = np.zeros(values.shape[:2], dtype=bool)
    family_rows = []
    for (category, target), ids in sorted(families.items()):
        cells = [groups[k]['cell'] for k in ids]
        if len(set(cells)) != len(cells):
            raise ValueError(f'Duplicate cell-target groups: {category}, {target}')
        count = len(cells)
        family_rows.append({'category': category, 'target': target, 'n_cells': count,
                            'source_cells': cells, 'eligible': count >= 2})
        if count < 2:
            # No cross-cell comparator: retain actual input; not evidence of zero importance.
            continue
        mean = values[:, ids].mean(axis=1)
        sd = values[:, ids].std(axis=1, ddof=0)
        region_varying = sd > 1e-12
        for k in ids:
            reference[:, pairs[k]] = mean
            eligible[k] = True
            varying[:, k] = region_varying.any(axis=-1)
            z[:, k] = np.divide(values[:, k]-mean, sd,
                               out=np.zeros_like(mean), where=region_varying)
    # If neither region varies, cross-cell z is undefined: retain NA in plots.
    z_pair = z.mean(axis=-1)
    z_pair[~varying] = np.nan
    return reference, z_pair, eligible, varying, family_rows
