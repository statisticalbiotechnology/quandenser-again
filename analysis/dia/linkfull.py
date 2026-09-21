#!/usr/bin/env python3
"""
Quant-first linking across DIA runs, against the full field of detected MS1
features rather than only the identified ones.

Everything is library-free: features come from MS1 peak picking, signatures
from correlating window fragment traces against each feature's own MS1 trace.
Identity enters only to label which run-2 feature is the right answer, and to
split results by whether a search engine found the peptide at all.
"""
import sys, pickle, collections
import numpy as np

R1 = '18484_REP3_1ug_Ecoli_NewStock2_SWATH_1'
R2 = '18486_REP3_1ug_Ecoli_NewStock2_SWATH_2'
RCUT, NSIG, MZTOL = 0.5, 20, 0.05
PPM = 10e-6


def load(run):
    f = pickle.load(open(f'ms1feat_{run}.pkl', 'rb'))
    b = pickle.load(open(f'block_{run}.pkl', 'rb'))
    lo = np.array([a for a, _ in b['windows']])
    hi = np.array([c for _, c in b['windows']])
    return f['feat'], f['traces'], f['rts'], b['cells'], b['block'], lo, hi


def signatures(run):
    feat, traces, rts, cells, block, lo_w, hi_w = load(run)
    wi = np.searchsorted(hi_w, feat[:, 0], 'right')
    ok = (wi < len(lo_w)) & (feat[:, 0] >= lo_w[np.clip(wi, 0, len(lo_w) - 1)])
    sigs = [None] * len(feat)
    for i in np.where(ok)[0]:
        key = (int(wi[i]), int(feat[i, 1] // block))
        cell = cells.get(key)
        if cell is None:
            continue
        cmz, C, crt = cell
        if C.shape[0] < 5:
            continue
        sc, it = traces[i]
        if len(sc) < 3:
            continue
        ftr_rt = rts[sc]
        keep = (crt >= feat[i, 2] - 10) & (crt <= feat[i, 3] + 10)
        if keep.sum() < 5:
            continue
        ref = np.interp(crt[keep], ftr_rt, it.astype(float))
        M = C[keep].astype(float)
        r = ref - ref.mean()
        dr = np.sqrt((r * r).sum())
        if dr <= 0:
            continue
        Mc = M - M.mean(0)
        dn = np.sqrt((Mc * Mc).sum(0))
        good = dn > 0
        rr = np.full(M.shape[1], -1.0)
        rr[good] = (r[:, None] * Mc[:, good]).sum(0) / (dr * dn[good])
        inten = M.sum(0)
        sel = np.where((rr >= RCUT) & (inten > 0))[0]
        if sel.size == 0:
            continue
        sel = sel[np.argsort(rr[sel])[::-1][:NSIG]]
        o = np.argsort(cmz[sel])
        sigs[i] = (cmz[sel][o], inten[sel][o])
    n = sum(1 for s in sigs if s)
    print(f'{run}: {len(feat)} features, {n} with a signature '
          f'({100*n/len(feat):.0f}%)')
    return feat, sigs, int(wi.max()) + 1, wi


def sim(a, b):
    ma, ia = a; mb, ib = b
    j = np.searchsorted(mb, ma)
    lo = np.clip(j - 1, 0, len(mb) - 1); hi = np.clip(j, 0, len(mb) - 1)
    dl = np.abs(mb[lo] - ma); dh = np.abs(mb[hi] - ma)
    k = np.where(dl < dh, lo, hi); d = np.minimum(dl, dh)
    hit = d <= MZTOL
    if not hit.any():
        return 0.0
    wa = np.sqrt(ia); wb = np.sqrt(ib)
    num = (wa[hit] * wb[k[hit]]).sum()
    da = np.sqrt((wa ** 2).sum()); db = np.sqrt((wb ** 2).sum())
    return float(num / (da * db)) if da > 0 and db > 0 else 0.0


f1, s1, nw, w1 = signatures(R1)
f2, s2, _, w2 = signatures(R2)

# identified peptides, to label the right answer
def idents(run):
    rows = open(f'split_{run}_mprophet_all_peakgroups_DDA.xls').read().rstrip('\n').split('\n')
    h = rows[0].split('\t'); ix = {k: i for i, k in enumerate(h)}
    out = {}
    for r in rows[1:]:
        c = r.split('\t')
        if c[ix['decoy']].strip().upper() != 'FALSE':
            continue
        try:
            q = float(c[ix['m_score']])
            if q > 0.01:
                continue
            key = c[ix['Sequence']] + '/' + c[ix['Charge']]
            rec = (float(c[ix['m.z']]), float(c[ix['RT']]), q)
        except (ValueError, IndexError):
            continue
        if key not in out or q < out[key][2]:
            out[key] = rec
    return out


i1, i2 = idents(R1), idents(R2)
shared = sorted(set(i1) & set(i2))
print(f'peptides confident in both runs: {len(shared)}')


def nearest(feat, sigs, mz, rt, dt=20.0):
    tol = max(mz * PPM * 2, 0.01)
    c = np.where((np.abs(feat[:, 0] - mz) <= tol) & (np.abs(feat[:, 1] - rt) <= dt))[0]
    c = [i for i in c if sigs[i] is not None]
    if not c:
        return None
    return c[int(np.argmax(feat[c, 4]))]


# index run-2 features by window for fast competitor lookup
byw2 = collections.defaultdict(list)
for i in range(len(f2)):
    if s2[i] is not None:
        byw2[int(w2[i])].append(i)
for w in byw2:
    byw2[w].sort(key=lambda i: f2[i, 1])

res = []
for key in shared:
    qi = nearest(f1, s1, i1[key][0], i1[key][1])
    ti = nearest(f2, s2, i2[key][0], i2[key][1])
    if qi is None or ti is None:
        continue
    w = int(w1[qi])
    members = byw2.get(w, [])
    rts2 = np.array([f2[i, 1] for i in members])
    a = np.searchsorted(rts2, f1[qi, 1] - 30.0)
    b = np.searchsorted(rts2, f1[qi, 1] + 30.0)
    comp = members[a:b]
    if ti not in comp:
        comp = list(comp) + [ti]
    if len(comp) < 2:
        continue
    sc = np.array([sim(s1[qi], s2[i]) for i in comp])
    order = np.argsort(sc)[::-1]
    rank = int(np.where(np.array(comp)[order] == ti)[0][0]) + 1
    res.append((len(comp), rank, f1[qi, 4], sc[list(comp).index(ti)],
                np.max(np.delete(sc, list(comp).index(ti)))))

res = np.array(res)
print(f'\n{len(res)} peptides detected as MS1 features with signatures in both runs')
print(f'competitors per query: median {np.median(res[:,0]):.0f}, '
      f'q90 {np.percentile(res[:,0],90):.0f}, max {res[:,0].max():.0f}')
print(f'top-1 {100*np.mean(res[:,1]==1):.1f}%   top-3 {100*np.mean(res[:,1]<=3):.1f}%'
      f'   chance {100*np.mean(1/res[:,0]):.1f}%')
q = np.percentile(res[:, 2], [25, 50, 75])
print(f'\n{"stratum":<16}{"n":>7}{"comp":>7}{"top-1":>9}{"top-3":>9}{"chance":>9}'
      f'{"sim +":>9}{"best -":>9}')
for name, m in (('weak   (Q1)', res[:, 2] <= q[0]),
                ('lower mid (Q2)', (res[:, 2] > q[0]) & (res[:, 2] <= q[1])),
                ('upper mid (Q3)', (res[:, 2] > q[1]) & (res[:, 2] <= q[2])),
                ('strong (Q4)', res[:, 2] > q[2]),
                ('all', np.ones(len(res), bool))):
    if m.sum() < 20:
        continue
    print(f'{name:<16}{int(m.sum()):>7}{np.median(res[m,0]):>7.0f}'
          f'{100*np.mean(res[m,1]==1):>8.1f}%{100*np.mean(res[m,1]<=3):>8.1f}%'
          f'{100*np.mean(1/res[m,0]):>8.1f}%{np.median(res[m,3]):>9.3f}'
          f'{np.median(res[m,4]):>9.3f}')
