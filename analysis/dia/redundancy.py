#!/usr/bin/env python3
"""
How many DISTINCT molecules are behind the confidently linked unidentified
features?

The feature detector used here has no isotope grouping and no charge
deconvolution, so one peptide appears several times: its isotope peaks, split
traces and shoulders are all separate features, and they all link across runs
because they are all the same molecule.  Counting features would overstate the
gain.  Features are clustered when their m/z differ by an isotope spacing at
some charge (or agree) and their apices are within 10 s.
"""
import pickle
import numpy as np

exec(open('yield.py').read().split('idx = np.where(have_sig)[0]')[0])

rng = np.random.default_rng(0)
idx = np.where(have_sig)[0]
if len(idx) > 30000:
    idx = rng.choice(idx, 30000, replace=False)
link = np.array([best_link(i) for i in idx])
null = np.array([best_link(i, exclude_near=300.0) for i in idx])
thr = np.percentile(null[null > 0], 95)
isid = ident_mask[idx]
conf = link >= thr

sel = idx[conf]
selid = isid[conf]
mz = f1[sel, 0]; rt = f1[sel, 1]
o = np.argsort(mz)
sel, selid, mz, rt = sel[o], selid[o], mz[o], rt[o]
print(f'{len(sel)} confidently linked features '
      f'({int(selid.sum())} identified, {int((~selid).sum())} not)')

parent = np.arange(len(sel))
def find(a):
    while parent[a] != a:
        parent[a] = parent[parent[a]]; a = parent[a]
    return a
def union(a, b):
    ra, rb = find(a), find(b)
    if ra != rb:
        parent[rb] = ra

C13 = 1.0033548
spacings = sorted({0.0} | {k * C13 / z for z in (1, 2, 3, 4) for k in (1, 2, 3)})
for i in range(len(sel)):
    tol = max(mz[i] * 10e-6 * 2, 0.01)
    j = i + 1
    while j < len(sel) and mz[j] - mz[i] <= max(spacings) + tol:
        if abs(rt[j] - rt[i]) <= 10.0:
            d = mz[j] - mz[i]
            if any(abs(d - s) <= tol for s in spacings):
                union(i, j)
        j += 1
roots = np.array([find(i) for i in range(len(sel))])
import collections
cl = collections.defaultdict(list)
for i, r in enumerate(roots):
    cl[r].append(i)
sizes = np.array([len(v) for v in cl.values()])
print(f'{len(cl)} clusters, median size {np.median(sizes):.0f}, '
      f'mean {sizes.mean():.2f}, max {sizes.max()}')

pure_unid = sum(1 for v in cl.values() if not selid[v].any())
mixed = sum(1 for v in cl.values() if selid[v].any())
print(f'\nclusters containing an identified feature : {mixed}')
print(f'clusters with no identified feature      : {pure_unid}')
print(f'ratio unidentified : identified          = {pure_unid/max(mixed,1):.1f}x')
print(f'\n(feature-level ratio was 8.1x; collapsing isotopes and split traces '
      f'brings it to {pure_unid/max(mixed,1):.1f}x)')
