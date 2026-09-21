#!/usr/bin/env python3
"""
Can library-free pseudo-spectra link DIA features across runs?

Positives  : the same peptide's signature in run 1 and run 2.
Hard negatives: signatures of two DIFFERENT features that share an isolation
             window AND co-elute.  These draw their candidate peaks from the
             same merged spectrum, so if the MS1-correlation filter does not
             make a signature specific to its own feature, they will look
             identical and the method collapses to m/z plus retention time.
Easy negatives: different features in the same window at any retention time.

Identity is used only to label pairs, never to build a signature.
"""
import sys, pickle, collections
import numpy as np

R1 = 'pseudo_18484_REP3_1ug_Ecoli_NewStock2_SWATH_1.pkl'
R2 = 'pseudo_18486_REP3_1ug_Ecoli_NewStock2_SWATH_2.pkl'
RCUT = 0.5          # a fragment joins the signature above this correlation
NSIG = 20           # at most this many fragments per signature
MZTOL = 0.05


def pearson_cols(ref, M):
    r = ref - ref.mean()
    d = np.sqrt((r * r).sum())
    if d <= 0:
        return np.full(M.shape[1], np.nan)
    Mc = M - M.mean(0)
    dn = np.sqrt((Mc * Mc).sum(0))
    out = np.full(M.shape[1], np.nan)
    ok = dn > 0
    out[ok] = (r[:, None] * Mc[:, ok]).sum(0) / (d * dn[ok])
    return out


def build(f):
    ms2rt, C = f['ms2rt'], f['ctr']
    if len(ms2rt) < 5 or C.size == 0:
        return None
    w = f['hi'] - f['lo']
    half = max(1.5 * w, 15.0)
    keep = (ms2rt >= f['rt'] - half) & (ms2rt <= f['rt'] + half)
    if keep.sum() < 5:
        return None
    if len(f['rts']) < 2:
        return None
    ref = np.interp(ms2rt[keep], f['rts'], f['ms1']).astype(float)
    M = C[keep].astype(float)
    r = pearson_cols(ref, M)
    inten = M.sum(0)
    sel = np.where((r >= RCUT) & (inten > 0))[0]
    if sel.size == 0:
        return None
    sel = sel[np.argsort(r[sel])[::-1][:NSIG]]
    o = np.argsort(f['cand'][sel])
    return f['cand'][sel][o], inten[sel][o], float(np.median(r[sel]))


def sim(a, b, idf=None):
    """Normalised dot product on sqrt intensities over matched m/z."""
    ma, ia = a[0], a[1]
    mb, ib = b[0], b[1]
    j = np.searchsorted(mb, ma)
    lo = np.clip(j - 1, 0, len(mb) - 1); hi = np.clip(j, 0, len(mb) - 1)
    dl = np.abs(mb[lo] - ma); dh = np.abs(mb[hi] - ma)
    k = np.where(dl < dh, lo, hi); d = np.minimum(dl, dh)
    hit = d <= MZTOL
    if not hit.any():
        return 0.0
    wa = np.sqrt(ia); wb = np.sqrt(ib)
    if idf is not None:
        wgt = idf(ma[hit])
        num = (wa[hit] * wb[k[hit]] * wgt).sum()
        da = np.sqrt((wa ** 2 * idf(ma)).sum()); db = np.sqrt((wb ** 2 * idf(mb)).sum())
    else:
        num = (wa[hit] * wb[k[hit]]).sum()
        da = np.sqrt((wa ** 2).sum()); db = np.sqrt((wb ** 2).sum())
    return float(num / (da * db)) if da > 0 and db > 0 else 0.0


a = {f['key']: f for f in pickle.load(open(R1, 'rb'))}
b = {f['key']: f for f in pickle.load(open(R2, 'rb'))}
sa = {k: v for k, v in ((k, build(f)) for k, f in a.items()) if v}
sb = {k: v for k, v in ((k, build(f)) for k, f in b.items()) if v}
shared = sorted(set(sa) & set(sb))
print(f'signatures built: run1 {len(sa)}/{len(a)}, run2 {len(sb)}/{len(b)}, '
      f'both {len(shared)}')
ns = [len(sa[k][0]) for k in shared]
print(f'fragments per signature: median {np.median(ns):.0f}, '
      f'q10 {np.percentile(ns,10):.0f}, q90 {np.percentile(ns,90):.0f}')

# global RT offset between runs, from the shared set
off = np.median([b[k]['rt'] - a[k]['rt'] for k in shared])
print(f'median RT offset run2 - run1: {off:+.1f} s\n')

# rarity weighting, MaRaCluster style: a fragment m/z seen in many signatures
# carries little evidence
BIN = 0.05
cnt = collections.Counter()
for k in sa:
    for m in sa[k][0]:
        cnt[round(m / BIN)] += 1
N = max(len(sa), 1)
def idf(mzs):
    return np.array([np.log(N / (1.0 + cnt.get(round(m / BIN), 0))) for m in mzs])

bywin = collections.defaultdict(list)
for k in shared:
    bywin[a[k]['w']].append(k)

pos, neg_hard, neg_easy = [], [], []
pos_i, neg_hard_i, neg_easy_i = [], [], []
rng = np.random.default_rng(0)
for k in shared:
    pos.append(sim(sa[k], sb[k])); pos_i.append(sim(sa[k], sb[k], idf))
    members = bywin[a[k]['w']]
    hard = [q for q in members if q != k
            and abs((a[q]['rt']) - a[k]['rt']) <= 30.0]
    easy = [q for q in members if q != k
            and abs((a[q]['rt']) - a[k]['rt']) > 300.0]
    if hard:
        q = hard[rng.integers(len(hard))]
        neg_hard.append(sim(sa[k], sb[q])); neg_hard_i.append(sim(sa[k], sb[q], idf))
    if easy:
        q = easy[rng.integers(len(easy))]
        neg_easy.append(sim(sa[k], sb[q])); neg_easy_i.append(sim(sa[k], sb[q], idf))


def auc(p, n):
    p = np.asarray(p); n = np.asarray(n)
    if not len(p) or not len(n):
        return np.nan, np.nan, np.nan
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


print(f'{"comparison":<44}{"n":>7}{"med pos":>9}{"med neg":>9}'
      f'{"AUC":>8}{"TPR@5%":>9}{"TPR@1%":>9}')
for name, P, Nn in (
        ('same peptide vs co-eluting same-window', pos, neg_hard),
        ('  ... with rarity weighting', pos_i, neg_hard_i),
        ('same peptide vs same-window elsewhere', pos, neg_easy),
        ('  ... with rarity weighting', pos_i, neg_easy_i)):
    A, t5, t1 = auc(P, Nn)
    print(f'{name:<44}{len(Nn):>7}{np.median(P):>9.3f}{np.median(Nn):>9.3f}'
          f'{A:>8.3f}{100*t5:>8.1f}%{100*t1:>8.1f}%')
