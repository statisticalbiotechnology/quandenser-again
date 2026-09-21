#!/usr/bin/env python3
"""
Cross-run pooling, with the retention-time warp actually estimated.

The two replicates were acquired back to back and still differ by about 290 s,
varying with retention time.  A fixed tolerance matches almost nothing, which
is the concrete reason naive cross-run fragment correlation is reported to work
only on homogeneous datasets.  The warp is estimated here from mass-matched
features: a coarse global mode, then a median refinement per retention-time
bin, then matching within a tight window around the predicted shift.
"""
import pickle, random
import numpy as np

R1 = '20150513_20_hTau_CSF_Frac9of12_rep1'
R2 = '20150513_48_hTau_CSF_Frac9of12_rep2'
WINDOW_MZ, MIN_PTS, MAX_NEG = 12.5, 5, 3
C13 = 1.0033548
HALF = 12.0
CYCLES = [(1, 0.44), (7, 3.1), (16, 7.0)]
BIN = 150.0


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


def load(run):
    d = pickle.load(open(f'/home/user/csf/calib/xic_{run}.pkl', 'rb'))
    f = d['feats']
    for k, fn in (('mz', lambda x: x['mz']), ('apex', lambda x: x['apex']),
                  ('fwhm', lambda x: x['fwhm']), ('inten', lambda x: x['inten']),
                  ('mass', lambda x: x['mz'] * x['z'] - x['z'] * 1.00728)):
        d[k] = np.array([fn(x) for x in f])
    return d


def candidates(a, b):
    ob = np.argsort(b['mass']); mbs = b['mass'][ob]
    out = []
    for i in range(len(a['mass'])):
        m = a['mass'][i]; tol = max(m * 10e-6, 0.01)
        l = np.searchsorted(mbs, m - tol); r = np.searchsorted(mbs, m + tol)
        for k in range(l, r):
            out.append((i, ob[k], b['apex'][ob[k]] - a['apex'][i]))
    return out


def estimate_warp(a, cands):
    dt = np.array([c[2] for c in cands])
    rt = np.array([a['apex'][c[0]] for c in cands])
    h, e = np.histogram(dt, bins=np.arange(-1800, 1801, 20))
    g = (e[h.argmax()] + e[h.argmax() + 1]) / 2.0
    edges = np.arange(rt.min(), rt.max() + BIN, BIN)
    centres, shifts = [], []
    for lo in edges[:-1]:
        m = (rt >= lo) & (rt < lo + BIN) & (np.abs(dt - g) <= 120)
        if m.sum() >= 30:
            centres.append(lo + BIN / 2); shifts.append(np.median(dt[m]))
    print(f'global modal shift {g:+.0f} s; warp estimated in '
          f'{len(centres)} bins, range {min(shifts):+.0f} to {max(shifts):+.0f} s')
    return np.array(centres), np.array(shifts), g


def main():
    a1, a2 = load(R1), load(R2)
    cands = candidates(a1, a2)
    cx, cy, g = estimate_warp(a1, cands)
    pred = lambda t: np.interp(t, cx, cy, left=cy[0], right=cy[-1])
    match = {}
    for i, j, dt in cands:
        r = abs(dt - pred(a1['apex'][i]))
        if r <= 20.0 and (i not in match or r < match[i][1]):
            match[i] = (j, r)
    match = {i: v[0] for i, v in match.items()}
    print(f'{len(match)} of {len(a1["mass"])} rep1 features matched to rep2 '
          f'after warping\n')

    keep = [i for i in match if a1['fwhm'][i] > 0 and a2['fwhm'][match[i]] > 0]
    keepset = set(keep)
    q = np.percentile(a1['inten'][keep], [25, 75])
    rng = random.Random(0)
    oz = np.argsort(a1['mz']); mzs = a1['mz'][oz]
    pairs = []
    for a in keep:
        l = np.searchsorted(mzs, a1['mz'][a] - WINDOW_MZ)
        r = np.searchsorted(mzs, a1['mz'][a] + WINDOW_MZ)
        cand = [oz[k] for k in range(l, r)
                if oz[k] != a and oz[k] in keepset
                and abs(a1['apex'][oz[k]] - a1['apex'][a])
                <= 0.5 * max(a1['fwhm'][a], a1['fwhm'][oz[k]])
                and not same_molecule(a1['mass'][a], a1['mass'][oz[k]])]
        rng.shuffle(cand)
        pairs += [(a, b) for b in cand[:MAX_NEG]]
    print(f'{len(keep)} matched features usable, {len(pairs)} interference pairs\n')

    def grid(d, i, step):
        lo = np.searchsorted(d['rts'], d['apex'][i] - HALF)
        hi = np.searchsorted(d['rts'], d['apex'][i] + HALF)
        return slice(lo, hi, step)

    strata = [('weak   (Q1)', lambda v: v <= q[0]),
              ('strong (Q4)', lambda v: v > q[1]),
              ('all', lambda v: np.ones(len(v), bool))]
    for step, cyc in CYCLES:
        res = {}
        for tag in ('one', 'two'):
            pos = np.full(len(keep), np.nan)
            for k, a in enumerate(keep):
                s1 = grid(a1, a, step)
                x = a1['mono'][a][s1].astype(float)
                y = a1['iso1'][a][s1].astype(float)
                if tag == 'two':
                    j = match[a]; s2 = grid(a2, j, step)
                    n = min(len(x), len(a2['mono'][j][s2]))
                    x = np.concatenate([x[:n], a2['mono'][j][s2][:n].astype(float)])
                    y = np.concatenate([y[:n], a2['iso1'][j][s2][:n].astype(float)])
                if len(x) >= MIN_PTS:
                    pos[k] = pearson(x, y)
            neg, src = [], []
            for a, b in pairs:
                s1 = grid(a1, a, step)
                x = a1['mono'][a][s1].astype(float)
                y = a1['mono'][b][s1].astype(float)
                if tag == 'two':
                    ja, jb = match[a], match[b]
                    sa, sb = grid(a2, ja, step), grid(a2, jb, step)
                    n = min(len(x), len(a2['mono'][ja][sa]), len(a2['mono'][jb][sb]))
                    x = np.concatenate([x[:n], a2['mono'][ja][sa][:n].astype(float)])
                    y = np.concatenate([y[:n], a2['mono'][jb][sb][:n].astype(float)])
                if len(x) >= MIN_PTS:
                    neg.append(pearson(x, y)); src.append(a)
            res[tag] = (pos, np.array(neg), np.array(src, dtype=int))
        print(f'cycle {cyc:.1f} s')
        print(f'    {"stratum":<16}{"n pos":>7}{"n neg":>8}{"AUC 1 run":>11}'
              f'{"AUC 2 runs":>12}{"TPR 1":>8}{"TPR 2":>8}')
        for name, sel in strata:
            out = []
            for tag in ('one', 'two'):
                pos, neg, src = res[tag]
                p = pos[sel(a1['inten'][keep])]
                n = neg[sel(a1['inten'][src])] if len(src) else neg
                pp = p[~np.isnan(p)]; nn = n[~np.isnan(n)]
                if len(nn) < 10 or len(pp) < 10:
                    out.append((np.nan, np.nan, len(pp), len(nn))); continue
                thr = np.percentile(nn, 95)
                out.append((auc(p, n), np.mean(pp >= thr), len(pp), len(nn)))
            print(f'    {name:<16}{out[0][2]:>7}{out[0][3]:>8}{out[0][0]:>11.3f}'
                  f'{out[1][0]:>12.3f}{100*out[0][1]:>7.1f}%{100*out[1][1]:>7.1f}%')
        print()


if __name__ == '__main__':
    main()
