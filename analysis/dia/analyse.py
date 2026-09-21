#!/usr/bin/env python3
"""Both channels on real DIA data, and the independence question."""
import pickle
import numpy as np
from scipy.stats import spearmanr

PROTON = 1.00727646
SHIFTS = np.array([s * d for d in range(6, 61) for s in (-1, 1)], float)
NDEC = len(SHIFTS)
TOL = 0.05


def pearson(a, b):
    a = a - a.mean(); b = b - b.mean()
    da = np.sqrt((a * a).sum()); db = np.sqrt((b * b).sum())
    return np.nan if da <= 0 or db <= 0 else float((a * b).sum() / (da * db))


def auc(pos, neg):
    pos = np.asarray(pos); neg = np.asarray(neg)
    pos = pos[~np.isnan(pos)]; neg = neg[~np.isnan(neg)]
    if not len(pos) or not len(neg):
        return np.nan
    allv = np.concatenate([pos, neg])
    o = np.argsort(allv, kind='mergesort'); sv = allv[o]
    rr = np.empty(len(allv)); i = 0
    while i < len(sv):
        j = i
        while j + 1 < len(sv) and sv[j + 1] == sv[i]:
            j += 1
        rr[o[i:j + 1]] = (i + j) / 2.0 + 1.0
        i = j + 1
    return (rr[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2.0) / (len(pos) * len(neg))


def count_pairs(mz, target, tol=TOL):
    lo = np.searchsorted(mz, target - tol - mz, 'left')
    hi = np.searchsorted(mz, target + tol - mz, 'right')
    tot = int((hi - lo).sum())
    self_hits = int(np.sum(np.abs(2 * mz - target) <= tol))
    return (tot - self_hits) // 2


peps = pickle.load(open('dia_traces.pkl', 'rb'))
print(f'{len(peps)} peptides\n')

# ---- channel 1: fragment trace vs MS1 precursor trace ----
tp, dp, per_pep = [], [], []
for p in peps:
    ms1 = p['ms1']
    T, D = p['tr_true'], p['tr_decoy']
    if ms1.size < 5 or T.size == 0:
        per_pep.append(np.nan); continue
    n = min(len(ms1), len(T))
    if n < 5:
        per_pep.append(np.nan); continue
    # Correlate over twice the elution width, not the padded +/-30 s window.
    # A window that is mostly baseline inflates the correlation of any two
    # traces, true or decoy, because they agree on being zero.
    # OpenSWATH's integration bounds are narrow here (~7 s against a ~3 s
    # cycle), so the correlation window is widened to give enough points while
    # staying well short of the padded collection range.
    w = p['hi'] - p['lo']
    half = max(1.5 * w, 15.0)
    rts = p['rts'][:n]
    keep = (rts >= p['rt'] - half) & (rts <= p['rt'] + half)
    if keep.sum() < 5:
        per_pep.append(np.nan); continue
    T = T[:n][keep]
    D = D[:n][keep] if D.size else D
    ref = ms1[:n][keep].astype(float)
    n = int(keep.sum())
    rs = [pearson(ref, T[:n, k].astype(float)) for k in range(T.shape[1])]
    tp += [r for r in rs if not np.isnan(r)]
    per_pep.append(np.nanmedian(rs) if len(rs) else np.nan)
    if D.size:
        m = min(n, len(D))
        dp += [r for r in (pearson(ref[:m], D[:m, k].astype(float))
                           for k in range(D.shape[1])) if not np.isnan(r)]
per_pep = np.array(per_pep)
print('CHANNEL 1  fragment trace vs MS1 precursor trace')
print(f'  true fragments  n={len(tp):6d}  median r {np.median(tp):+.3f}')
print(f'  decoy fragments n={len(dp):6d}  median r {np.median(dp):+.3f}')
a1 = auc(tp, dp)
thr = np.percentile(dp, 95)
print(f'  AUC {a1:.3f}   TPR@5%FPR {100*np.mean(np.array(tp) >= thr):.1f}%'
      f'   (threshold r >= {thr:.3f})')
thr1 = np.percentile(dp, 99)
print(f'                  TPR@1%FPR {100*np.mean(np.array(tp) >= thr1):.1f}%'
      f'   (threshold r >= {thr1:.3f})\n')

# ---- channel 2: complementary pairs on the merged window spectrum ----
c2, npk = np.full(len(peps), np.nan), np.zeros(len(peps))
tgt = np.full(len(peps), np.nan); dec = np.full(len(peps), np.nan)
for i, p in enumerate(peps):
    mz = np.sort(p['spec_mz'])
    npk[i] = len(mz)
    if len(mz) < 10:
        continue
    M = p['mz'] * p['z'] - p['z'] * PROTON
    s0 = M + 2 * PROTON
    t = count_pairs(mz, s0)
    d = np.fromiter((count_pairs(mz, s0 + s) for s in SHIFTS), int, NDEC)
    tgt[i] = t; dec[i] = d.mean()
    c2[i] = (np.sum(d >= t) + 1) / (NDEC + 1)
ok = ~np.isnan(c2)
print('CHANNEL 2  complementary pairs, merged window spectrum')
print(f'  peptides scored {int(ok.sum())}, median peaks per merged spectrum '
      f'{np.median(npk[ok]):.0f}')
print(f'  target count mean {np.nanmean(tgt):.2f}   decoy mean {np.nanmean(dec):.2f}'
      f'   ratio {np.nanmean(tgt)/max(np.nanmean(dec),1e-9):.2f}')
beats = c2[ok] <= 1.0 / (NDEC + 1)
print(f'  beats all {NDEC} decoy masses: {100*beats.mean():.1f}%'
      f'   (null {100/(NDEC+1):.2f}%)   enrichment {beats.mean()*(NDEC+1):.1f}x\n')

# ---- independence ----
xc = np.array([p['xcorr'] for p in peps])
m = ok & ~np.isnan(xc)
rho = spearmanr(xc[m], -np.log10(c2[m])).statistic
print('INDEPENDENCE')
print(f'  rank correlation, OpenSWATH xcorr_shape vs channel-2 evidence: '
      f'rho = {rho:+.3f}  (n={int(m.sum())})')
m2 = ok & ~np.isnan(per_pep)
rho2 = spearmanr(per_pep[m2], -np.log10(c2[m2])).statistic
print(f'  rank correlation, own channel-1 score vs channel-2 evidence:   '
      f'rho = {rho2:+.3f}  (n={int(m2.sum())})')
q = np.percentile(xc[m], [33, 67])
print('\n  channel-2 hit rate split by channel-1 strength:')
for name, sel in (('weak co-elution   ', xc[m] <= q[0]),
                  ('middle            ', (xc[m] > q[0]) & (xc[m] <= q[1])),
                  ('strong co-elution ', xc[m] > q[1])):
    b = c2[m][sel] <= 1.0 / (NDEC + 1)
    print(f'    {name} n={int(sel.sum()):5d}   beats-all {100*b.mean():5.1f}%')
