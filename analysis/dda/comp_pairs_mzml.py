#!/usr/bin/env python3
"""
Same calibration, but on raw DDA MS2 spectra rather than MaRaCluster consensus
spectra, which are truncated to 160 peaks and therefore cannot carry the
high-m/z partner of a complementary pair for higher-charge precursors.

Generalised pair relation.  For fragments of charge c1 and c2,
    c1*mz_1 + c2*mz_2 = M + (c1+c2)*proton
so (1,1) gives mz_1 + mz_2 = M + 2p, and (1,2) gives mz_1 + 2*mz_2 = M + 3p.
Adding (1,2) matters for 3+ and 4+ precursors, whose fragments are often 2+.
The decoy shifts see exactly the same relation, so the null stays valid.
"""
import sys, collections
import numpy as np
from pyteomics import mzml

PROTON = 1.00727646
SHIFTS = np.array([-43,-41,-37,-31,-29,-23,-19,-13,-9,-6,
                     6,  9, 13, 19, 23, 29, 31, 37, 41, 43], float)
TOL = 0.01
TOPN = int(sys.argv[2]) if len(sys.argv) > 2 else 0     # 0 = keep all peaks


def count_11(mz, target, tol):
    lo = np.searchsorted(mz, target - tol - mz, 'left')
    hi = np.searchsorted(mz, target + tol - mz, 'right')
    total = int((hi - lo).sum())
    self_hits = int(np.sum(np.abs(2 * mz - target) <= tol))
    return (total - self_hits) // 2


def count_12(mz, target, tol):
    """mz_i + 2*mz_j = target, i != j.  Ordered, so no halving."""
    twice = 2.0 * mz
    lo = np.searchsorted(mz, target - tol - twice, 'left')
    hi = np.searchsorted(mz, target + tol - twice, 'right')
    total = int((hi - lo).sum())
    self_hits = int(np.sum(np.abs(3 * mz - target) <= tol))
    return total - self_hits


def load(path):
    out = []
    for s in mzml.read(path):
        if s.get('ms level') != 2:
            continue
        pre = s['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]
        z = pre.get('charge state')
        if z is None:
            continue
        mzp = float(pre['selected ion m/z'])
        mz = np.asarray(s['m/z array'], float)
        inten = np.asarray(s['intensity array'], float)
        if mz.size < 4:
            continue
        if TOPN and mz.size > TOPN:
            keep = np.argpartition(inten, -TOPN)[-TOPN:]
            mz = mz[keep]
        out.append((int(z), (mzp - PROTON) * int(z), np.sort(mz)))
    return out


def report(spectra, label, counter, name):
    n = len(spectra)
    tgt = np.zeros(n, int)
    dec = np.zeros((n, len(SHIFTS)), int)
    zs = np.zeros(n, int)
    for i, (z, M, mz) in enumerate(spectra):
        zs[i] = z
        s0 = M + (3 if name == '(1,2)' else 2) * PROTON
        tgt[i] = counter(mz, s0, TOL)
        for k, d in enumerate(SHIFTS):
            dec[i, k] = counter(mz, s0 + d, TOL)
    beats = (dec >= tgt[:, None]).sum(1)
    top = beats == 0
    null = 1.0 / (len(SHIFTS) + 1)
    print(f'  {name} relation, {label}')
    print(f'    {"z":<3}{"n":>7}{"target":>9}{"decoy":>9}{"ratio":>8}'
          f'{"beats-all":>11}{"excess pp":>11}')
    for z in sorted(set(zs.tolist())):
        m = zs == z
        if m.sum() < 100:
            continue
        t, d = tgt[m].mean(), dec[m].mean()
        print(f'    {z:<3}{m.sum():>7}{t:>9.2f}{d:>9.2f}'
              f'{t/max(d,1e-9):>8.1f}{100*top[m].mean():>10.2f}%'
              f'{100*(top[m].mean()-null):>10.2f}')
    t, d = tgt.mean(), dec.mean()
    print(f'    {"all":<3}{n:>7}{t:>9.2f}{d:>9.2f}{t/max(d,1e-9):>8.1f}'
          f'{100*top.mean():>10.2f}%{100*(top.mean()-null):>10.2f}')
    return tgt, dec, zs


if __name__ == '__main__':
    path = sys.argv[1]
    spectra = load(path)
    lab = f'top-{TOPN} peaks' if TOPN else 'all peaks'
    print(f'{path.split("/")[-1]}: {len(spectra)} MS2 spectra with a charge state, '
          f'{lab}, tol +/- {TOL} Da, null beats-all rate 4.76%')
    report(spectra, lab, count_11, '(1,1)')
    print()
    report(spectra, lab, count_12, '(1,2)')
