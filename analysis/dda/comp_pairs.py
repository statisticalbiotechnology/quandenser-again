#!/usr/bin/env python3
"""
Calibration of the complementary-pair statistic on DDA data.

For a peptide of neutral monoisotopic mass M, singly charged b and y ions
satisfy  m(b_i) + m(y_{n-i}) = M + 2*1.00728, independent of i.  The .ms2 Z
line carries (M + 1.00728), so the target sum is simply Z + 1.00728.

For every spectrum we count peak pairs summing to the target within a
tolerance, and repeat at 20 decoy sums obtained by shifting the target by an
integer number of Da.  Integer shifts preserve the nominal-mass phase of the
pair-sum distribution, which is periodic at 1 Da because peptide fragment
masses follow the averagine line; a non-integer shift would compare against a
different phase and manufacture a signal.

The per-spectrum permutation p-value is
    p = (#decoys with count >= target + 1) / 21
which is uniform on {1/21, ..., 21/21} under the null.
"""
import sys, collections
import numpy as np

PROTON = 1.00727646
MS2 = sys.argv[1] if len(sys.argv) > 1 else \
    '/home/user/csf/results/quandenser/consensus_spectra/Quandenser.consensus.part1.ms2'
SHIFTS = np.array([-43,-41,-37,-31,-29,-23,-19,-13,-9,-6,
                     6,  9, 13, 19, 23, 29, 31, 37, 41, 43], float)
TOLS = [0.005, 0.01, 0.02, 0.05, 0.10, 0.30]


def read_ms2(path):
    z = mh = None
    mz = []
    for line in open(path):
        c = line[0]
        if c == 'S':
            if z is not None and mz:
                yield z, mh, np.array(mz)
            z = mh = None
            mz = []
        elif c == 'Z':
            f = line.split()
            z, mh = int(f[1]), float(f[2])
        elif c in 'HID':
            continue
        else:
            mz.append(float(line.split(None, 1)[0]))
    if z is not None and mz:
        yield z, mh, np.array(mz)


def count_pairs(sorted_mz, target, tol):
    """Unordered peak pairs (i<j) with mz_i + mz_j within tol of target."""
    lo = np.searchsorted(sorted_mz, target - tol - sorted_mz, side='left')
    hi = np.searchsorted(sorted_mz, target + tol - sorted_mz, side='right')
    total = int((hi - lo).sum())                       # ordered, incl. self
    half = target / 2.0
    self_hits = int(np.sum(np.abs(2 * sorted_mz - target) <= tol))
    return (total - self_hits) // 2


def main():
    spectra = [(z, mh, np.sort(mz)) for z, mh, mz in read_ms2(MS2)]
    print(f'spectra: {len(spectra)}')
    print(f'decoy shifts: {len(SHIFTS)} integer Da offsets in '
          f'[{SHIFTS.min():.0f}, {SHIFTS.max():.0f}]\n')

    for tol in TOLS:
        tgt = np.zeros(len(spectra), int)
        dec = np.zeros((len(spectra), len(SHIFTS)), int)
        zs = np.zeros(len(spectra), int)
        for i, (z, mh, mz) in enumerate(spectra):
            zs[i] = z
            s0 = mh + PROTON                      # = M + 2*proton
            tgt[i] = count_pairs(mz, s0, tol)
            for k, d in enumerate(SHIFTS):
                dec[i, k] = count_pairs(mz, s0 + d, tol)

        beats = (dec >= tgt[:, None]).sum(1)
        pval = (beats + 1) / (len(SHIFTS) + 1)
        top = pval <= 1.0 / (len(SHIFTS) + 1)     # target strictly beats all decoys
        null_rate = 1.0 / (len(SHIFTS) + 1)

        print(f'--- tolerance +/- {tol:g} Da ---')
        print(f'  mean count   target {tgt.mean():7.3f}   decoy {dec.mean():7.3f}'
              f'   ratio {tgt.mean()/max(dec.mean(),1e-9):6.2f}')
        print(f'  spectra where target beats all 20 decoys: '
              f'{top.sum():5d} / {len(spectra)}  = {100*top.mean():5.2f}%'
              f'   (null {100*null_rate:.2f}%)')
        excess = top.mean() - null_rate
        print(f'  excess over null: {100*excess:5.2f} percentage points  '
              f'= {int(round(excess*len(spectra)))} spectra')
        print('  by precursor charge:')
        for z in sorted(set(zs.tolist())):
            m = zs == z
            if m.sum() < 50:
                continue
            print(f'    z={z}  n={m.sum():5d}  target {tgt[m].mean():6.2f}'
                  f'  decoy {dec[m].mean():6.2f}'
                  f'  beats-all {100*top[m].mean():5.2f}%')
        print()

    # p-value histogram at the tolerance that looked best
    return spectra


if __name__ == '__main__':
    main()
