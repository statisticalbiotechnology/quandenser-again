#!/usr/bin/env python3
"""
Where a library-free signature could actually earn its place in match-between-
runs, given that on back-to-back replicates nearest-retention-time already gets
95% right.

A. Alignment residual.  These two runs differ by 2.5 s.  Real experiments do
   not: the two CSF replicates in this session, run back to back on the same
   instrument, differ by about 290 s varying with retention time, and an
   aligner leaves a residual.  Retention-time evidence degrades with that
   residual; the signature does not.  A perturbation is added to the query's
   retention time and the crossover located.

B. Link validation.  Choosing among candidates is not the only job.  Match-
   between-runs also has to decide whether a proposed link is real at all,
   which is what Quandenser's linkPEP is for.  Positives are true cross-run
   pairs; negatives are features at the same m/z and a different retention
   time, i.e. exactly the false link that m/z alone would accept.  Retention
   time cannot score this without already trusting the alignment; the
   signature can.
"""
import pickle, collections
import numpy as np

R1 = '18484_REP3_1ug_Ecoli_NewStock2_SWATH_1'
R2 = '18486_REP3_1ug_Ecoli_NewStock2_SWATH_2'
MZTOL, PPM, C13 = 0.05, 10e-6, 1.0033548

exec(open('linkfull.py').read().split('# index run-2 features by window')[0])

o2 = np.argsort(f2[:, 0]); mz2 = f2[o2, 0]
pairs = []
for key in shared:
    qi = nearest(f1, s1, i1[key][0], i1[key][1])
    ti = nearest(f2, s2, i2[key][0], i2[key][1])
    if qi is not None and ti is not None:
        pairs.append((qi, ti))
print(f'{len(pairs)} peptides detected with a signature in both runs\n')


def same_mol(mq, mc):
    t = max(mq * PPM * 2, 0.01)
    if abs(mq - mc) <= t:
        return True
    return any(abs(abs(mq - mc) - k * C13 / z) <= t
               for z in (1, 2, 3, 4) for k in (1, 2, 3))


def candidates(qi, ti, rtw, shift=0.0):
    lo = np.searchsorted(mz2, f1[qi, 0] * (1 - PPM * 2))
    hi = np.searchsorted(mz2, f1[qi, 0] * (1 + PPM * 2))
    qrt = f1[qi, 1] + shift
    c = [o2[k] for k in range(lo, hi)
         if s2[o2[k]] is not None and abs(f2[o2[k], 1] - qrt) <= rtw]
    c = [x for x in c if x == ti or not same_mol(f1[qi, 0], f2[x, 0])
         or abs(f2[x, 1] - f2[ti, 1]) > 20.0]
    if ti not in c:
        c.append(ti)
    return c, qrt


print('A. effect of alignment residual (candidate window 120 s)')
print(f'{"residual":>10}{"n amb":>8}{"cands":>7}{"RT-nearest":>12}'
      f'{"signature":>11}{"chance":>9}')
rng = np.random.default_rng(0)
for res_s in (0.0, 10.0, 20.0, 40.0, 60.0):
    hr, hs, ch, nc = [], [], [], []
    for qi, ti in pairs:
        sh = rng.normal(0, res_s) if res_s > 0 else 0.0
        c, qrt = candidates(qi, ti, 120.0, sh)
        if len(c) < 2:
            continue
        nc.append(len(c)); ch.append(1.0 / len(c))
        dt = np.array([abs(f2[x, 1] - qrt) for x in c])
        hr.append(c[int(dt.argmin())] == ti)
        sc = np.array([sim(s1[qi], s2[x]) for x in c])
        hs.append(c[int(sc.argmax())] == ti)
    print(f'{res_s:>9.0f}s{len(nc):>8}{np.median(nc):>7.0f}'
          f'{100*np.mean(hr):>11.1f}%{100*np.mean(hs):>10.1f}%'
          f'{100*np.mean(ch):>8.1f}%')


def auc(p, n):
    p = np.asarray(p); n = np.asarray(n)
    allv = np.concatenate([p, n])
    o = np.argsort(allv, kind='mergesort'); sv = allv[o]
    rr = np.empty(len(allv)); i = 0
    while i < len(sv):
        j = i
        while j + 1 < len(sv) and sv[j + 1] == sv[i]:
            j += 1
        rr[o[i:j + 1]] = (i + j) / 2.0 + 1.0
        i = j + 1
    A = (rr[:len(p)].sum() - len(p) * (len(p) + 1) / 2.0) / (len(p) * len(n))
    return A, np.mean(p >= np.percentile(n, 95)), np.mean(p >= np.percentile(n, 99))


print('\nB. link validation: true cross-run pair vs a same-m/z false link')
pos, neg, posint = [], [], []
for qi, ti in pairs:
    pos.append(sim(s1[qi], s2[ti])); posint.append(f1[qi, 4])
    lo = np.searchsorted(mz2, f1[qi, 0] * (1 - PPM * 2))
    hi = np.searchsorted(mz2, f1[qi, 0] * (1 + PPM * 2))
    wrong = [o2[k] for k in range(lo, hi)
             if s2[o2[k]] is not None
             and abs(f2[o2[k], 1] - f2[ti, 1]) > 120.0]
    if wrong:
        neg.append(sim(s1[qi], s2[wrong[rng.integers(len(wrong))]]))
A, t5, t1 = auc(pos, neg)
print(f'  n pos {len(pos)}, n neg {len(neg)}')
print(f'  median similarity: true {np.median(pos):.3f}, false {np.median(neg):.3f}')
print(f'  AUC {A:.3f}   TPR@5%FPR {100*t5:.1f}%   TPR@1%FPR {100*t1:.1f}%')
posint = np.array(posint); pos = np.array(pos)
q = np.percentile(posint, [25, 75])
print('\n  by feature intensity:')
for name, m in (('weak   (Q1)', posint <= q[0]),
                ('mid           ', (posint > q[0]) & (posint <= q[1])),
                ('strong (Q4)', posint > q[1])):
    thr5 = np.percentile(neg, 95)
    print(f'    {name:<16}n={int(m.sum()):>5}   median sim {np.median(pos[m]):.3f}'
          f'   above 5% FPR threshold {100*np.mean(pos[m] >= thr5):5.1f}%')
