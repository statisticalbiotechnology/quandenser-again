#!/usr/bin/env python3
"""
Minimal MS1 feature detection on a SWATH run, to give the linking test a
realistic field of competitors rather than only the peptides a search engine
identified.

Traces are built in fixed logarithmic m/z bins (10 ppm wide), not by splitting
the pooled peak list on gaps: with 20 million peaks over 800 m/z the pooled
axis has no gaps at all, and gap splitting silently returns a handful of giant
traces.  Two interleaved bin sets, offset by half a bin, so a trace sitting on
a boundary is still caught in one of them.

No isotope grouping and no charge deconvolution: for linking, the unit is an
(m/z, RT) signal with a reproducible fragment signature.
"""
import sys, pickle
import numpy as np
from pyteomics import mzxml

RUN = sys.argv[1]
PPM = 10e-6
MINSCANS = 4
PCTL = 90.0           # per-scan intensity percentile kept
LOGSTEP = np.log1p(PPM)

rts, mzs, ints = [], [], []
for s in mzxml.read(RUN + '.mzXML'):
    if s['msLevel'] != 1:
        continue
    mz = np.asarray(s['m/z array'], np.float64)
    it = np.asarray(s['intensity array'], np.float32)
    if mz.size:
        thr = np.percentile(it, PCTL)
        k = it >= thr
        mz, it = mz[k], it[k]
    rts.append(float(s['retentionTime']) * 60.0)
    mzs.append(mz); ints.append(it)
rts = np.array(rts)
npk = np.array([len(m) for m in mzs])
print(f'{RUN}: {len(rts)} MS1 scans, kept {npk.sum()/1e6:.1f}M peaks '
      f'(top {100-PCTL:.0f}% per scan), median {np.median(npk):.0f}/scan')

allmz = np.concatenate(mzs); allint = np.concatenate(ints)
allscan = np.repeat(np.arange(len(rts)), npk)
base = np.log(allmz.min())
feats = []
traces = []
for off in (0.0, 0.5):
    b = np.floor((np.log(allmz) - base) / LOGSTEP + off).astype(np.int64)
    o = np.lexsort((allscan, b))
    bs, ss, is_, ms = b[o], allscan[o], allint[o], allmz[o]
    edges = np.concatenate([[0], np.where(np.diff(bs) != 0)[0] + 1, [len(bs)]])
    for i in range(len(edges) - 1):
        a, z = edges[i], edges[i + 1]
        if z - a < MINSCANS:
            continue
        sc, it, mz = ss[a:z], is_[a:z], ms[a:z]
        brk = np.where(np.diff(sc) > 2)[0]
        for seg in np.split(np.arange(len(sc)), brk + 1):
            if len(seg) < MINSCANS:
                continue
            sit = it[seg]
            k = seg[sit.argmax()]
            feats.append((float(np.average(mz[seg], weights=sit)),
                          float(rts[sc[k]]), float(rts[sc[seg[0]]]),
                          float(rts[sc[seg[-1]]]), float(sit.max()),
                          int(len(seg)), off))
            traces.append((sc[seg].astype(np.int32), sit.astype(np.float32)))
arr = np.array(feats)
print(f'{len(arr)} raw traces from both bin offsets')

# deduplicate: same m/z within 10 ppm and apex within 10 s, keep the strongest
o = np.lexsort((-arr[:, 4], arr[:, 0]))
arr = arr[o]
traces = [traces[i] for i in o]
keep = np.ones(len(arr), bool)
for i in range(len(arr)):
    if not keep[i]:
        continue
    j = i + 1
    while j < len(arr) and arr[j, 0] - arr[i, 0] <= arr[i, 0] * PPM * 2:
        if keep[j] and abs(arr[j, 1] - arr[i, 1]) <= 10.0:
            keep[j] = False
        j += 1
traces = [t for t, k in zip(traces, keep) if k]
arr = arr[keep][:, :6]
print(f'{len(arr)} MS1 features after deduplication')
print(f'  m/z {arr[:,0].min():.1f}-{arr[:,0].max():.1f}, '
      f'width median {np.median(arr[:,3]-arr[:,2]):.1f} s, '
      f'scans median {np.median(arr[:,5]):.0f}, '
      f'apex intensity median {np.median(arr[:,4]):.3g}')
pickle.dump(dict(feat=arr, traces=traces, rts=rts),
            open(f'ms1feat_{RUN}.pkl', 'wb'), 4)
print(f'saved ms1feat_{RUN}.pkl')
