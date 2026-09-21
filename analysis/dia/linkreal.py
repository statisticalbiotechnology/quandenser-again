#!/usr/bin/env python3
"""
The linking task as a pipeline actually poses it.

Candidates for a run-1 feature are run-2 features within a few ppm of its m/z
and inside a retention-time window, not everything sharing a 25 m/z isolation
window.  Two questions:

  1. How often is m/z plus retention time already unambiguous, and when it is
     not, does the fragment signature pick the right candidate?
  2. How far can the retention-time window be widened before m/z alone fails,
     and does the signature hold up where it does?  That is what matters for
     match-between-runs across many runs, where the retention-time search has
     to be generous.

Candidates within an isotope spacing of the query are dropped: they are the
same molecule at another m/z, they have the same MS1 trace and therefore the
same signature, and counting them as wrong answers measures nothing.
"""
import pickle, collections
import numpy as np

R1 = '18484_REP3_1ug_Ecoli_NewStock2_SWATH_1'
R2 = '18486_REP3_1ug_Ecoli_NewStock2_SWATH_2'
MZTOL, PPM, C13 = 0.05, 10e-6, 1.0033548
RTWINS = [30.0, 60.0, 120.0, 300.0, 600.0]

exec(open('linkfull.py').read().split('# index run-2 features by window')[0])

o2 = np.argsort(f2[:, 0])
mz2 = f2[o2, 0]


def same_molecule(mq, mc):
    d = abs(mq - mc)
    if d <= max(mq * PPM * 2, 0.01):
        return True
    for z in (1, 2, 3, 4):
        for k in (1, 2, 3):
            if abs(d - k * C13 / z) <= max(mq * PPM * 2, 0.01):
                return True
    return False


pairs = []
for key in shared:
    qi = nearest(f1, s1, i1[key][0], i1[key][1])
    ti = nearest(f2, s2, i2[key][0], i2[key][1])
    if qi is not None and ti is not None:
        pairs.append((qi, ti))
print(f'{len(pairs)} peptides detected with a signature in both runs\n')

print(f'{"RT window":>10}{"cands":>8}{"ambiguous":>11}{"m/z+RT top-1":>14}'
      f'{"+signature":>12}{"chance":>9}   (all)')
for rtw in RTWINS:
    nc, hit_rt, hit_sig, chance, amb = [], [], [], [], []
    for qi, ti in pairs:
        lo = np.searchsorted(mz2, f1[qi, 0] * (1 - PPM * 2))
        hi = np.searchsorted(mz2, f1[qi, 0] * (1 + PPM * 2))
        cand = [o2[k] for k in range(lo, hi)
                if s2[o2[k]] is not None
                and abs(f2[o2[k], 1] - f1[qi, 1]) <= rtw]
        cand = [c for c in cand
                if c == ti or not same_molecule(f1[qi, 0], f2[c, 0])
                or abs(f2[c, 1] - f2[ti, 1]) > 20.0]
        if ti not in cand:
            cand.append(ti)
        nc.append(len(cand))
        amb.append(len(cand) > 1)
        chance.append(1.0 / len(cand))
        # baseline: closest in retention time
        dt = np.array([abs(f2[c, 1] - f1[qi, 1]) for c in cand])
        hit_rt.append(cand[int(dt.argmin())] == ti)
        sc = np.array([sim(s1[qi], s2[c]) for c in cand])
        hit_sig.append(cand[int(sc.argmax())] == ti)
    nc = np.array(nc)
    print(f'{rtw:>9.0f}s{np.median(nc):>8.0f}{100*np.mean(amb):>10.1f}%'
          f'{100*np.mean(hit_rt):>13.1f}%{100*np.mean(hit_sig):>11.1f}%'
          f'{100*np.mean(chance):>8.1f}%')

print(f'\nrestricted to queries where m/z + RT is genuinely ambiguous '
      f'(more than one candidate):')
print(f'{"RT window":>10}{"n":>7}{"cands":>8}{"RT-nearest":>12}'
      f'{"signature":>11}{"chance":>9}')
for rtw in RTWINS:
    nc, hit_rt, hit_sig, chance = [], [], [], []
    for qi, ti in pairs:
        lo = np.searchsorted(mz2, f1[qi, 0] * (1 - PPM * 2))
        hi = np.searchsorted(mz2, f1[qi, 0] * (1 + PPM * 2))
        cand = [o2[k] for k in range(lo, hi)
                if s2[o2[k]] is not None
                and abs(f2[o2[k], 1] - f1[qi, 1]) <= rtw]
        cand = [c for c in cand
                if c == ti or not same_molecule(f1[qi, 0], f2[c, 0])
                or abs(f2[c, 1] - f2[ti, 1]) > 20.0]
        if ti not in cand:
            cand.append(ti)
        if len(cand) < 2:
            continue
        nc.append(len(cand)); chance.append(1.0 / len(cand))
        dt = np.array([abs(f2[c, 1] - f1[qi, 1]) for c in cand])
        hit_rt.append(cand[int(dt.argmin())] == ti)
        sc = np.array([sim(s1[qi], s2[c]) for c in cand])
        hit_sig.append(cand[int(sc.argmax())] == ti)
    if len(nc) < 20:
        continue
    print(f'{rtw:>9.0f}s{len(nc):>7}{np.median(nc):>8.0f}'
          f'{100*np.mean(hit_rt):>11.1f}%{100*np.mean(hit_sig):>10.1f}%'
          f'{100*np.mean(chance):>8.1f}%')
