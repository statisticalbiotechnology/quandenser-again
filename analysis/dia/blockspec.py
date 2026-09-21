#!/usr/bin/env python3
"""
Candidate fragment traces per (isolation window, retention-time block).

The candidate pool is a property of the window and the time, not of any one
feature, so every feature in a block is scored against the same pool and all
the specificity has to come from the MS1-correlation filter.  That is the thing
being tested, so it must not be helped by giving each feature its own pool.

Two streaming passes: one to find the most intense peaks of each block's merged
spectrum, one to extract their traces.
"""
import sys, pickle, collections
import numpy as np
from pyteomics import mzxml

RUN = sys.argv[1]
BLOCK = 60.0
TOPK = 150
FRAG_TOL = 0.05

windows, seen, n1 = [], set(), 0
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
lo_w = np.array([a for a, _ in windows]); hi_w = np.array([b for _, b in windows])
print(f'{len(windows)} windows of {hi_w[0]-lo_w[0]:.0f} m/z')


def win_of(c):
    i = np.searchsorted(hi_w, c, 'right')
    return i if i < len(windows) and lo_w[i] <= c < hi_w[i] else -1


# pass 1: merged spectrum per block
acc = collections.defaultdict(lambda: ([], []))
for s in mzxml.read(RUN + '.mzXML'):
    if s['msLevel'] != 2:
        continue
    rt = float(s['retentionTime']) * 60.0
    w = win_of(float(s['precursorMz'][0]['precursorMz']))
    if w < 0:
        continue
    a, b = acc[(w, int(rt // BLOCK))]
    a.append(np.asarray(s['m/z array'], np.float64))
    b.append(np.asarray(s['intensity array'], np.float32))
print(f'{len(acc)} (window, block) cells')

cand = {}
for key, (mzl, itl) in acc.items():
    mz = np.concatenate(mzl); it = np.concatenate(itl)
    o = np.argsort(mz); mz, it = mz[o], it[o]
    grp = np.concatenate([[0], np.cumsum(np.diff(mz) > FRAG_TOL)])
    n = grp[-1] + 1
    wsum = np.bincount(grp, weights=mz * it, minlength=n)
    isum = np.bincount(grp, weights=it, minlength=n)
    k = isum > 0
    mmz, mit = wsum[k] / isum[k], isum[k]
    if len(mmz) > TOPK:
        sel = np.argpartition(mit, -TOPK)[-TOPK:]
        mmz = mmz[sel]
    cand[key] = np.sort(mmz)
del acc
print(f'candidate lists built, median {np.median([len(v) for v in cand.values()]):.0f} m/z each')


def pick(smz, sint, targets):
    j = np.searchsorted(smz, targets)
    a = np.clip(j - 1, 0, smz.size - 1); b = np.clip(j, 0, smz.size - 1)
    da = np.abs(smz[a] - targets); db = np.abs(smz[b] - targets)
    k = np.where(da < db, a, b); d = np.minimum(da, db)
    return np.where(d <= FRAG_TOL, sint[k], 0.0).astype(np.float32)


# pass 2: traces on those candidates
tr = collections.defaultdict(lambda: ([], []))
for s in mzxml.read(RUN + '.mzXML'):
    if s['msLevel'] != 2:
        continue
    rt = float(s['retentionTime']) * 60.0
    w = win_of(float(s['precursorMz'][0]['precursorMz']))
    if w < 0:
        continue
    key = (w, int(rt // BLOCK))
    c = cand.get(key)
    if c is None or c.size == 0:
        continue
    smz = np.asarray(s['m/z array'], np.float64)
    sint = np.asarray(s['intensity array'], np.float32)
    if smz.size == 0:
        continue
    a, b = tr[key]
    a.append(pick(smz, sint, c)); b.append(rt)
out = {k: (cand[k], np.array(v[0]), np.array(v[1])) for k, v in tr.items()}
print(f'{len(out)} cells with traces, median scans per cell '
      f'{np.median([len(v[2]) for v in out.values()]):.0f}')
pickle.dump(dict(cells=out, windows=windows, block=BLOCK),
            open(f'block_{RUN}.pkl', 'wb'), 4)
print(f'saved block_{RUN}.pkl')
