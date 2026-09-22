#!/usr/bin/env python3
"""
Independence, using the denoised channel-2 score where it has signal.

If channel 2's hits fall on the peptides channel 1 already gets right, the two
are redundant and combining them buys nothing.  If they fall preferentially
where channel 1 is weak, the combination is worth building.
"""
import sys, pickle
import numpy as np
from scipy.stats import spearmanr

TRACES = sys.argv[1] if len(sys.argv) > 1 else 'dia_traces.pkl'

PROTON = 1.00727646
SHIFTS = np.array([s * d for d in range(6, 61) for s in (-1, 1)], float)
NDEC = len(SHIFTS); NULL = 1.0 / (NDEC + 1)
TOPN, TOL = 100, 0.02


def pearson(a, b):
    a = a - a.mean(); b = b - b.mean()
    da = np.sqrt((a * a).sum()); db = np.sqrt((b * b).sum())
    return np.nan if da <= 0 or db <= 0 else float((a * b).sum() / (da * db))


def count_pairs(mz, target, tol=TOL):
    lo = np.searchsorted(mz, target - tol - mz, 'left')
    hi = np.searchsorted(mz, target + tol - mz, 'right')
    tot = int((hi - lo).sum())
    return (tot - int(np.sum(np.abs(2 * mz - target) <= tol))) // 2


peps = pickle.load(open(TRACES, 'rb'))
c1_own = np.full(len(peps), np.nan)
c1_osw = np.array([p['xcorr'] for p in peps])
c2_p = np.full(len(peps), np.nan)
c2_excess = np.full(len(peps), np.nan)

for i, p in enumerate(peps):
    ms1, T = p['ms1'], p['tr_true']
    if ms1.size >= 5 and T.size:
        n = min(len(ms1), len(T))
        w = p['hi'] - p['lo']; half = max(1.5 * w, 15.0)
        keep = (p['rts'][:n] >= p['rt'] - half) & (p['rts'][:n] <= p['rt'] + half)
        if keep.sum() >= 5:
            ref = ms1[:n][keep].astype(float)
            Tk = T[:n][keep]
            rs = [pearson(ref, Tk[:, k].astype(float)) for k in range(Tk.shape[1])]
            rs = [r for r in rs if not np.isnan(r)]
            if rs:
                c1_own[i] = np.median(rs)
    mz, it = p['spec_mz'], p['spec_int']
    if len(mz) >= 10:
        if len(mz) > TOPN:
            mz = mz[np.argpartition(it, -TOPN)[-TOPN:]]
        mz = np.sort(mz)
        M = p['mz'] * p['z'] - p['z'] * PROTON
        s0 = M + 2 * PROTON
        t = count_pairs(mz, s0)
        d = np.fromiter((count_pairs(mz, s0 + s) for s in SHIFTS), int, NDEC)
        c2_p[i] = (np.sum(d >= t) + 1) / (NDEC + 1)
        c2_excess[i] = t - d.mean()

m = ~np.isnan(c1_own) & ~np.isnan(c2_p) & ~np.isnan(c1_osw)
print(f'{int(m.sum())} peptides with both channels scored\n')
hit = c2_p[m] <= NULL
print(f'channel-2 beats-all rate overall: {100*hit.mean():.1f}%  '
      f'(null {100*NULL:.2f}%, enrichment {hit.mean()/NULL:.1f}x)\n')

print('rank correlations (higher = more redundant):')
print(f'  OpenSWATH xcorr_shape  vs channel-2 evidence : '
      f'rho {spearmanr(c1_osw[m], -np.log10(c2_p[m])).statistic:+.3f}')
print(f'  own co-elution score   vs channel-2 evidence : '
      f'rho {spearmanr(c1_own[m], -np.log10(c2_p[m])).statistic:+.3f}')
print(f'  own co-elution score   vs channel-2 excess   : '
      f'rho {spearmanr(c1_own[m], c2_excess[m]).statistic:+.3f}')
print(f'  (reference) OpenSWATH xcorr vs own score     : '
      f'rho {spearmanr(c1_osw[m], c1_own[m]).statistic:+.3f}\n')

for label, score in (('OpenSWATH xcorr_shape', c1_osw[m]),
                     ('own co-elution score ', c1_own[m])):
    q = np.percentile(score, [33, 67])
    print(f'channel-2 hit rate split by {label}:')
    for name, sel in (('weak  ', score <= q[0]),
                      ('middle', (score > q[0]) & (score <= q[1])),
                      ('strong', score > q[1])):
        h = c2_p[m][sel] <= NULL
        print(f'    {name}  n={int(sel.sum()):5d}   beats-all {100*h.mean():5.1f}%'
              f'   enrichment {h.mean()/NULL:5.1f}x')
    print()
