#!/usr/bin/env python3
"""
Fair test for channel 2 on real DIA.

The merged window spectra carry ~3000 peaks, most of them noise, while the CSF
test that showed the statistic working used denoised 65-peak spectra.  Judging
channel 2 on the raw peak list would be blaming co-isolation for something
denoising fixes, so the statistic is rerun on the top-N most intense peaks.
"""
import pickle
import numpy as np

PROTON = 1.00727646
SHIFTS = np.array([s * d for d in range(6, 61) for s in (-1, 1)], float)
NDEC = len(SHIFTS)
NULL = 1.0 / (NDEC + 1)


def count_pairs(mz, target, tol):
    lo = np.searchsorted(mz, target - tol - mz, 'left')
    hi = np.searchsorted(mz, target + tol - mz, 'right')
    tot = int((hi - lo).sum())
    self_hits = int(np.sum(np.abs(2 * mz - target) <= tol))
    return (tot - self_hits) // 2


peps = pickle.load(open('dia_traces.pkl', 'rb'))
npk = np.array([len(p['spec_mz']) for p in peps])
print(f'{len(peps)} peptides, merged spectrum peaks: median {np.median(npk):.0f}, '
      f'q10 {np.percentile(npk,10):.0f}, q90 {np.percentile(npk,90):.0f}')
print(f'null beats-all rate {100*NULL:.2f}%\n')
print(f'{"top-N peaks":>12}{"tol":>7}{"used":>7}{"target":>9}{"decoy":>9}'
      f'{"ratio":>8}{"beats-all":>11}{"enrich":>8}')
for topn in (50, 100, 200, 500, 0):
    for tol in (0.02, 0.05):
        t_all, d_all, beats, used = [], [], [], []
        for p in peps:
            mz, it = p['spec_mz'], p['spec_int']
            if len(mz) < 10:
                continue
            if topn and len(mz) > topn:
                k = np.argpartition(it, -topn)[-topn:]
                mz = mz[k]
            mz = np.sort(mz)
            used.append(len(mz))
            M = p['mz'] * p['z'] - p['z'] * PROTON
            s0 = M + 2 * PROTON
            t = count_pairs(mz, s0, tol)
            d = np.fromiter((count_pairs(mz, s0 + s, tol) for s in SHIFTS), int, NDEC)
            t_all.append(t); d_all.append(d.mean()); beats.append((d >= t).sum() == 0)
        t_all = np.array(t_all); d_all = np.array(d_all); beats = np.array(beats)
        lab = f'{topn}' if topn else 'all'
        print(f'{lab:>12}{tol:>7.2f}{int(np.median(used)):>7}{t_all.mean():>9.2f}'
              f'{d_all.mean():>9.2f}{t_all.mean()/max(d_all.mean(),1e-9):>8.2f}'
              f'{100*beats.mean():>10.1f}%{beats.mean()/NULL:>8.1f}x')
