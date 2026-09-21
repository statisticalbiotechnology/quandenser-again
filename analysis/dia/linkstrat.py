#!/usr/bin/env python3
"""
The linking test as an actual retrieval task, stratified by intensity.

Retrieval: given a feature in run 1, rank every run-2 feature that shares its
isolation window and falls within +/-30 s by signature similarity, and ask
whether the correct one comes first.  This is the operation a quant-first
pipeline performs, and unlike an AUC it degrades honestly as competitors
multiply.

Intensity strata matter because the population this is meant to serve, the
features no search engine identified, is the weak one.  Numbers from the top
quartile say nothing useful about it.
"""
import pickle, collections
import numpy as np

exec(open('linktest.py').read().split("a = {f['key']")[0])   # reuse build/sim/pearson_cols

a = {f['key']: f for f in pickle.load(open(R1, 'rb'))}
b = {f['key']: f for f in pickle.load(open(R2, 'rb'))}
sa = {k: v for k, v in ((k, build(f)) for k, f in a.items()) if v}
sb = {k: v for k, v in ((k, build(f)) for k, f in b.items()) if v}
shared = sorted(set(sa) & set(sb))
inten = np.array([a[k]['spec_int_sum'] for k in shared])
q = np.percentile(inten, [25, 50, 75])
print(f'{len(shared)} features with a signature in both runs')
print(f'window-summed intensity quartiles: {q[0]:.3g}  {q[1]:.3g}  {q[2]:.3g}\n')

bywin = collections.defaultdict(list)
for k in shared:
    bywin[a[k]['w']].append(k)

ncomp, top1, top1_rand, sims_pos, sims_best_neg = [], [], [], [], []
for k in shared:
    comp = [q2 for q2 in bywin[a[k]['w']]
            if abs(b[q2]['rt'] - a[k]['rt']) <= 30.0]
    if k not in comp:
        comp.append(k)
    ncomp.append(len(comp))
    s = np.array([sim(sa[k], sb[q2]) for q2 in comp])
    order = np.argsort(s)[::-1]
    top1.append(comp[order[0]] == k)
    top1_rand.append(1.0 / len(comp))
    ki = comp.index(k)
    sims_pos.append(s[ki])
    others = np.delete(s, ki)
    sims_best_neg.append(others.max() if others.size else 0.0)

ncomp = np.array(ncomp); top1 = np.array(top1)
top1_rand = np.array(top1_rand)
sims_pos = np.array(sims_pos); sims_best_neg = np.array(sims_best_neg)
print(f'competitors per feature (same window, within 30 s): '
      f'median {np.median(ncomp):.0f}, q90 {np.percentile(ncomp,90):.0f}, '
      f'max {ncomp.max()}\n')

strata = [('weak   (Q1)', inten <= q[0]),
          ('lower mid (Q2)', (inten > q[0]) & (inten <= q[1])),
          ('upper mid (Q3)', (inten > q[1]) & (inten <= q[2])),
          ('strong (Q4)', inten > q[2]),
          ('all', np.ones(len(shared), bool))]
print(f'{"stratum":<16}{"n":>7}{"comp":>7}{"top-1":>9}{"chance":>9}'
      f'{"med sim +":>11}{"med best -":>12}{"margin>0":>10}')
for name, m in strata:
    if m.sum() < 20:
        continue
    marg = np.mean(sims_pos[m] > sims_best_neg[m])
    print(f'{name:<16}{int(m.sum()):>7}{np.median(ncomp[m]):>7.0f}'
          f'{100*top1[m].mean():>8.1f}%{100*top1_rand[m].mean():>8.1f}%'
          f'{np.median(sims_pos[m]):>11.3f}{np.median(sims_best_neg[m]):>12.3f}'
          f'{100*marg:>9.1f}%')

print('\nretrieval restricted to features with at least 3 competitors:')
m3 = ncomp >= 3
for name, m in strata:
    mm = m & m3
    if mm.sum() < 20:
        continue
    print(f'  {name:<16}n={int(mm.sum()):>5}   competitors median '
          f'{np.median(ncomp[mm]):.0f}   top-1 {100*top1[mm].mean():5.1f}%'
          f'   chance {100*top1_rand[mm].mean():4.1f}%')
