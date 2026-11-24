"""For spectra of run 0 whose instrument precursor is among Quandenser's
candidates: is that precursor in a printed feature group (run 0 row with the
same charge, m/z within 30 ppm, apex within 0.5 min)?

    pasef_lost_precursor.py <feature_groups.tsv> <ms2_header.tsv> <candidates.tsv> <consensus_hdr.tsv> <clusters.tsv> <run>
"""
import sys

import numpy as np
import pandas as pd

fg_fn, hdr_fn, cand_fn, cons_fn, cl_fn, run = sys.argv[1:7]
fg = pd.read_csv(fg_fn, sep='\t', header=None, names=['file', 'mz', 'z', 'rt', 'inten', 'links'],
                 skip_blank_lines=True).dropna(subset=['mz'])
fg = fg[fg.file == 0].sort_values('mz').reset_index(drop=True)
hdr = pd.read_csv(hdr_fn, sep='\t').merge(pd.read_csv(cand_fn, sep='\t')[['scan', 'true_in']], on='scan')
cl = pd.read_csv(cl_fn, sep='\t', header=None, names=['file', 'scan', 'cluster']).dropna()
cl = cl[cl.file == run][['scan', 'cluster']].astype(int)
cons = pd.read_csv(cons_fn, sep='\t', header=None, names=['cscan', 'cmz', 'cz'])
cons['cluster'] = cons.cscan // 100
hdr = hdr.merge(cl, on='scan')
m = hdr[['scan', 'cluster', 'precmz', 'charge']].merge(cons, on='cluster')
m['true'] = (m.cz == m.charge) & (np.abs(m.cmz - m.precmz) <= 30e-6 * m.precmz)
hdr['true_cons'] = hdr.scan.isin(m[m.true].scan)

mz, z, rt = fg.mz.to_numpy(), fg.z.to_numpy(), fg.rt.to_numpy()


def in_group(s):
    a, b = np.searchsorted(mz, [s.precmz * (1 - 30e-6), s.precmz * (1 + 30e-6)])
    return bool(((z[a:b] == s.charge) & (np.abs(rt[a:b] - s.rt_min) <= 0.5)).any())


lost = hdr[(hdr.true_in == 1) & ~hdr.true_cons]
grp = np.array([in_group(s) for s in lost.itertuples(index=False)])
n = len(hdr)
print(f'spectra {n}; instrument precursor among candidates {100 * (hdr.true_in == 1).mean():.1f}%; '
      f'carried by a consensus spectrum {100 * hdr.true_cons.mean():.1f}%')
print(f'among candidates but not carried: {100 * len(lost) / n:.1f}% of spectra')
print(f'  its feature is in a printed feature group (lost the intensity contest): {100 * grp.sum() / n:.1f}%')
print(f'  in no printed feature group (fails --max_missing or unmatched): {100 * (~grp).sum() / n:.1f}%')
