"""Regression checks for strand mapping and physical-consistent motif disruption."""
import numpy as np
import pandas as pd
from run_native import BASES, COMP, REGIONS, copies, orient, score, shuffle_core


def main():
    pwm = np.full((6, 4), .02)
    reference = np.array([0, 1, 2, 3, 0, 2], dtype=np.uint8)
    pwm[np.arange(6), reference] = .94
    changed = shuffle_core(reference, pwm, 8)
    assert changed is not None
    mutant, wt, delta = changed
    assert np.array_equal(np.bincount(reference, minlength=4), np.bincount(mutant, minlength=4))
    assert delta < wt and not np.array_equal(reference, mutant)
    assert shuffle_core(np.zeros(6, dtype=np.uint8), pwm, 8) is None
    assert np.array_equal(orient(orient(reference, True), True), reference)
    assert np.allclose(score(np.stack([reference, reference]), pwm), wt)

    for strand in ('+', '-'):
        row = {'gene_index': 0, 'strand': strand}
        coordinates = {'tss_sequence': (1000, 1600), 'first_donor_sequence': (800, 2800),
                       'first_acceptor_sequence': (1500, 3500), 'pre_tes_sequence': (1700, 2700)}
        codes = {}
        for name, (left, right) in coordinates.items():
            row[name + '_requested_start0'], row[name + '_requested_end0'] = left, right
            genome = np.arange(left, right)
            sequence = (genome % 4).astype(np.uint8)
            if strand == '-':
                sequence = COMP[sequence[::-1]]
            codes[name] = sequence[None, :]
        for name in REGIONS:
            mapping, signature = copies(pd.Series(row), name, 260, 269, codes)
            assert len(signature) == 9 and len(mapping) >= 9
            original = codes[name][0, 260:269]
            altered = original[::-1].copy()
            inputs = {key: a.copy() for key, a in codes.items()}
            for region, offset, j in mapping:
                assert inputs[region][0, offset] == original[j]
                inputs[region][0, offset] = altered[j]
            assert np.array_equal(inputs[name][0, 260:269], altered)
            for region, offset, j in mapping:
                assert inputs[region][0, offset] == altered[j]
        missing = {key: a.copy() for key, a in codes.items()}
        missing['first_acceptor_sequence'][:] = 4
        row['first_acceptor_sequence_requested_start0'] = np.nan
        row['first_acceptor_sequence_requested_end0'] = np.nan
        mapping, _ = copies(pd.Series(row), 'tss_sequence', 260, 269, missing)
        assert all(region != 'first_acceptor_sequence' for region, _, _ in mapping)
    print('PASS: composition, immutable cores, PWM score selection, strand round trip, shared physical edits')


if __name__ == '__main__':
    main()
