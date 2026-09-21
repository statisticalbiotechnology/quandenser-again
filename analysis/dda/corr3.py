#!/usr/bin/env python3
"""
Same measurement, controlled two more ways.

- All cycle times are evaluated on the SAME feature subset (those with at least
  5 points at the coarsest cycle), so a row is not flattered by the coarse
  cycles keeping only broad, well-behaved peaks.
- Results are split by feature intensity.  A weak trace correlates badly with
  its own isotope for reasons that have nothing to do with the method, and the
  weak population is precisely the one DIA is supposed to rescue, so an average
  over all intensities answers the wrong question.
"""
import sys, pickle, random
import numpy as np

RUN = sys.argv[1] if len(sys.argv) > 1 else '20150513_20_hTau_CSF_Frac9of12_rep1'
WINDOW_MZ, MIN_PTS, MAX_NEG = 12.5, 5, 3
C13 = 1.0033548
CYCLES = [(1, 0.44), (4, 1.8), (7, 3.1), (11, 4.8), (16, 7.0)]


def pearson(a, b):
    a = a - a.mean(); b = b - b.mean()
    da = np.sqrt((a * a).sum()); db = np.sqrt((b * b).sum())
    return np.nan if da <= 0 or db <= 0 else float((a * b).sum() / (da * db))


def auc(pos, neg):
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


def same_molecule(a, b):
    d = abs(a - b)
    if d < 0.02:
        return True
    k = round(d / C13)
    return k <= 3 and abs(d - k * C13) < 0.02


def main():
    d = pickle.load(open(f'/home/user/csf/calib/xic_{RUN}.pkl', 'rb'))
    feats, rts, mono, iso1 = d['feats'], d['rts'], d['mono'], d['iso1']
    mz = np.array([f['mz'] for f in feats])
    apex = np.array([f['apex'] for f in feats])
    fwhm = np.array([f['fwhm'] for f in feats])
    inten = np.array([f['inten'] for f in feats])
    mass = np.array([f['mz'] * f['z'] - f['z'] * 1.00728 for f in feats])
    span = np.maximum(1.5 * fwhm, 3.0)
    i0 = np.searchsorted(rts, apex - span)
    i1 = np.searchsorted(rts, apex + span)
    coarse = CYCLES[-1][0]
    idx = np.where((fwhm > 0) & (((i1 - i0) + coarse - 1) // coarse >= MIN_PTS))[0]
    print(f'{len(idx)} features usable at every cycle time tested')
    q = np.percentile(inten[idx], [25, 50, 75])
    print(f'intensity apex quartiles: {q[0]:.3g}  {q[1]:.3g}  {q[2]:.3g}\n')

    rng = random.Random(0)
    order = np.argsort(mz); mzs = mz[order]; okset = set(idx.tolist())
    pairs = []
    for a in idx:
        l = np.searchsorted(mzs, mz[a] - WINDOW_MZ)
        r = np.searchsorted(mzs, mz[a] + WINDOW_MZ)
        cand = [order[k] for k in range(l, r)
                if order[k] != a and order[k] in okset
                and abs(apex[order[k]] - apex[a]) <= 0.5 * max(fwhm[a], fwhm[order[k]])
                and not same_molecule(mass[a], mass[order[k]])]
        rng.shuffle(cand)
        pairs += [(a, b) for b in cand[:MAX_NEG]]
    print(f'{len(pairs)} interference pairs on that subset\n')

    strata = [('weak   (Q1)', lambda v: v <= q[0]),
              ('mid    (Q2-Q3)', lambda v: (v > q[0]) & (v <= q[2])),
              ('strong (Q4)', lambda v: v > q[2]),
              ('all', lambda v: np.ones_like(v, bool))]

    for step, cyc in CYCLES:
        pos = np.full(len(idx), np.nan)
        pi = {a: k for k, a in enumerate(idx)}
        for k, a in enumerate(idx):
            sl = slice(i0[a], i1[a], step)
            x = mono[a][sl].astype(float)
            if len(x) >= MIN_PTS:
                pos[k] = pearson(x, iso1[a][sl].astype(float))
        neg, negsrc = [], []
        for a, b in pairs:
            sl = slice(i0[a], i1[a], step)
            x = mono[a][sl].astype(float)
            if len(x) >= MIN_PTS:
                neg.append(pearson(x, mono[b][sl].astype(float)))
                negsrc.append(a)
        neg = np.array(neg); negsrc = np.array(negsrc)
        pts = int(np.median([(i1[a] - i0[a] + step - 1) // step for a in idx]))
        print(f'cycle {cyc:.1f} s, {pts} points per peak')
        print(f'    {"stratum":<16}{"n pos":>7}{"n neg":>8}{"pos med r":>11}'
              f'{"neg med r":>11}{"AUC":>8}{"TPR@5%FPR":>11}')
        for name, sel in strata:
            pm = sel(inten[idx])
            nm = sel(inten[negsrc])
            p, n = pos[pm], neg[nm]
            pp = p[~np.isnan(p)]; nn = n[~np.isnan(n)]
            if len(pp) < 20 or len(nn) < 20:
                continue
            thr = np.percentile(nn, 95)
            print(f'    {name:<16}{len(pp):>7}{len(nn):>8}{np.median(pp):>11.3f}'
                  f'{np.median(nn):>11.3f}{auc(p, n):>8.3f}'
                  f'{100*np.mean(pp >= thr):>10.1f}%')
        print()


if __name__ == '__main__':
    main()
