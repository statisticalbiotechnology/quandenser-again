#!/usr/bin/env python3
"""
The payoff number: how many features that no search engine identified can be
linked across runs with calibrated confidence, and therefore quantified.

The threshold is set on the false-link distribution measured in linkvalue.py
(query against a same-m/z feature more than 120 s away), so it is a real 5%
false-link rate, not a guess.  Yield is reported by intensity decile because
most of the 97k detected features are noise and an aggregate would hide that.
"""
import pickle
import numpy as np

R1 = '18484_REP3_1ug_Ecoli_NewStock2_SWATH_1'
R2 = '18486_REP3_1ug_Ecoli_NewStock2_SWATH_2'
MZTOL, PPM, C13 = 0.05, 10e-6, 1.0033548
RTW = 120.0

exec(open('linkfull.py').read().split('# index run-2 features by window')[0])

o2 = np.argsort(f2[:, 0]); mz2 = f2[o2, 0]
rng = np.random.default_rng(0)

# which run-1 features correspond to an identified peptide
ident_mask = np.zeros(len(f1), bool)
for key in i1:
    mz, rt, _ = i1[key]
    tol = max(mz * PPM * 2, 0.01)
    c = np.where((np.abs(f1[:, 0] - mz) <= tol) & (np.abs(f1[:, 1] - rt) <= 30.0))[0]
    ident_mask[c] = True
have_sig = np.array([s is not None for s in s1])
print(f'run1: {len(f1)} features, {have_sig.sum()} with a signature, '
      f'{int((ident_mask & have_sig).sum())} matching an identified peptide')


def best_link(qi, rtw=RTW, exclude_near=None):
    lo = np.searchsorted(mz2, f1[qi, 0] * (1 - PPM * 2))
    hi = np.searchsorted(mz2, f1[qi, 0] * (1 + PPM * 2))
    best = 0.0
    for k in range(lo, hi):
        c = o2[k]
        if s2[c] is None:
            continue
        dt = abs(f2[c, 1] - f1[qi, 1])
        if exclude_near is None:
            if dt > rtw:
                continue
        else:
            if dt <= exclude_near:
                continue
        v = sim(s1[qi], s2[c])
        if v > best:
            best = v
    return best


idx = np.where(have_sig)[0]
if len(idx) > 30000:
    idx = rng.choice(idx, 30000, replace=False)
print(f'scoring {len(idx)} run-1 features\n')

link = np.array([best_link(i) for i in idx])
# null: the same features linked against same-m/z partners far away in RT
null = np.array([best_link(i, exclude_near=300.0) for i in idx])
thr = np.percentile(null[null > 0], 95) if (null > 0).any() else 1.0
print(f'false-link null: {int((null>0).sum())} non-zero of {len(null)}, '
      f'95th percentile {thr:.3f}')

isid = ident_mask[idx]
inten = f1[idx, 4]
conf = link >= thr
print(f'\nconfidently linked at a 5% false-link rate:')
print(f'  identified   {int((conf & isid).sum()):>6} / {int(isid.sum()):>6}'
      f'  ({100*conf[isid].mean():.1f}%)')
print(f'  unidentified {int((conf & ~isid).sum()):>6} / {int((~isid).sum()):>6}'
      f'  ({100*conf[~isid].mean():.1f}%)')
r = (conf & ~isid).sum() / max((conf & isid).sum(), 1)
print(f'  unidentified per identified: {r:.1f}x')

dec = np.percentile(inten, np.arange(10, 100, 10))
bins = np.digitize(inten, dec)
print(f'\n{"intensity decile":>17}{"n":>8}{"n unid":>8}{"linked (unid)":>15}'
      f'{"linked (id)":>13}')
for d in range(10):
    m = bins == d
    if m.sum() < 20:
        continue
    mu = m & ~isid
    mi = m & isid
    lu = f'{100*conf[mu].mean():.1f}%' if mu.sum() >= 10 else '-'
    li = f'{100*conf[mi].mean():.1f}%' if mi.sum() >= 10 else '-'
    print(f'{d+1:>17}{int(m.sum()):>8}{int(mu.sum()):>8}{lu:>15}{li:>13}')
