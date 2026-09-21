#!/usr/bin/env python3
"""
Library-free pseudo-spectra for DIA features, for LINKING rather than
identification.

For each feature the candidate fragments are the most intense peaks of its
merged isolation-window spectrum, with no library involved.  A fragment joins
the feature's signature if its MS2 trace correlates with the feature's own MS1
precursor trace.  The signature will contain interference; that is expected and
does not matter here, because the question is whether the SAME signature comes
back in another run, not whether it is correct.

Pass 2 over the mzXML: pass 1 (diatest.py) produced the merged window spectra
that define the candidate m/z.
"""
import sys, pickle
import numpy as np
from pyteomics import mzxml

RUN = sys.argv[1]
TOPK = 150            # candidate fragments per feature
OUT = f'pseudo_{RUN}.pkl'
FRAG_TOL = 0.05

peps = pickle.load(open(f'traces_{RUN}.pkl', 'rb'))
print(f'{RUN}: {len(peps)} features')

windows = []
seen = set()
n1 = 0
for s in mzxml.read(RUN + '.mzXML'):
    if s['msLevel'] == 1:
        n1 += 1
        if n1 > 1:
            break
        continue
    p = s['precursorMz'][0]
    c = float(p['precursorMz']); w = float(p.get('windowWideness', 25.0))
    if (c, w) not in seen:
        seen.add((c, w)); windows.append((c - w / 2, c + w / 2))
windows.sort()


def win_of(mz):
    for i, (a, b) in enumerate(windows):
        if a <= mz < b:
            return i
    return -1


keep = []
for i, p in enumerate(peps):
    mz, it = p['spec_mz'], p['spec_int']
    if len(mz) < 20:
        p['cand'] = np.array([]); continue
    k = np.argpartition(it, -min(TOPK, len(it)))[-min(TOPK, len(it)):]
    p['cand'] = np.sort(mz[k])
    p['ctr'] = []
    p['ms2rt'] = []
    keep.append(i)
print(f'{len(keep)} features with a usable candidate list')

pad = 45.0
lo = np.array([p['lo'] - pad for p in peps])
hi = np.array([p['hi'] + pad for p in peps])
bywin = {}
for i in keep:
    bywin.setdefault(peps[i]['w'], []).append(i)


def pick(smz, sint, targets, tol):
    if smz.size == 0 or len(targets) == 0:
        return np.zeros(len(targets), np.float32)
    j = np.searchsorted(smz, targets)
    a = np.clip(j - 1, 0, smz.size - 1); b = np.clip(j, 0, smz.size - 1)
    da = np.abs(smz[a] - targets); db = np.abs(smz[b] - targets)
    k = np.where(da < db, a, b); d = np.minimum(da, db)
    return np.where(d <= tol, sint[k], 0.0).astype(np.float32)


n2 = 0
for s in mzxml.read(RUN + '.mzXML'):
    if s['msLevel'] != 2:
        continue
    n2 += 1
    rt = float(s['retentionTime']) * 60.0
    c = float(s['precursorMz'][0]['precursorMz'])
    w = win_of(c)
    if w < 0 or w not in bywin:
        continue
    act = [i for i in bywin[w] if lo[i] <= rt <= hi[i]]
    if not act:
        continue
    smz = np.asarray(s['m/z array'], float)
    sint = np.asarray(s['intensity array'], float)
    for i in act:
        peps[i]['ctr'].append(pick(smz, sint, peps[i]['cand'], FRAG_TOL))
        peps[i]['ms2rt'].append(rt)
print(f'streamed {n2} MS2 scans')

out = []
for i in keep:
    p = peps[i]
    if not p['ctr']:
        continue
    out.append(dict(key=p['key'], mz=p['mz'], z=p['z'], rt=p['rt'],
                    lo=p['lo'], hi=p['hi'], w=p['w'],
                    ms1=p['ms1'], rts=p['rts'],
                    cand=p['cand'], ctr=np.array(p['ctr']),
                    ms2rt=np.array(p['ms2rt']),
                    spec_int_sum=float(p['spec_int'].sum())))
pickle.dump(out, open(OUT, 'wb'), 4)
print(f'saved {OUT} ({len(out)} features)')
