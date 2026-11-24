"""Precursor candidates Quandenser gives each MS2 spectrum (src/SpectrumFiles.cpp:
features with m/z inside the isolation window and an RT span covering the
spectrum), against the instrument's own precursor.

    pasef_candidates.py <ms2_header.tsv> <features.tsv> <out.tsv>

"True" candidate: same charge as the instrument and m/z within 30 ppm of the
instrument's monoisotopic m/z.
"""
import sys

import numpy as np
import pandas as pd

hdr = pd.read_csv(sys.argv[1], sep='\t')
ft = pd.read_csv(sys.argv[2], sep='\t', usecols=['mz', 'charge', 'rtStart', 'rtEnd',
                                                  'intensityApex', 'ionMobility'])
ft = ft[ft.intensityApex > 0].sort_values('mz').reset_index(drop=True)
mz, z, rs, re_, im = (ft[c].to_numpy() for c in ['mz', 'charge', 'rtStart', 'rtEnd', 'ionMobility'])

rows = []
for s in hdr.itertuples(index=False):
    a, b = np.searchsorted(mz, [s.precmz - s.iso_lo, s.precmz + s.iso_hi])
    sl = slice(a, b)
    keep = (rs[sl] <= s.rt_min) & (re_[sl] >= s.rt_min)
    cz, cmz, cim = z[sl][keep], mz[sl][keep], im[sl][keep]
    dim = np.abs(cim - s.im)
    true = (cz == s.charge) & (np.abs(cmz - s.precmz) <= 30e-6 * s.precmz)
    rows.append((s.scan, s.charge, len(cz), int((cz == 1).sum()),
                 int((dim <= 0.05).sum()), int((dim <= 0.03).sum()),
                 int(true.any()), int((true & (dim <= 0.05)).any())))
out = pd.DataFrame(rows, columns=['scan', 'inst_charge', 'n_cand', 'n_cand_z1', 'n_cand_im05',
                                  'n_cand_im03', 'true_in', 'true_in_im05'])
out.to_csv(sys.argv[3], sep='\t', index=False)

n = len(out)
print(f'MS2 spectra: {n}')
print('instrument charge:', out.inst_charge.value_counts().sort_index().to_dict())
for col in ['n_cand', 'n_cand_im05', 'n_cand_im03']:
    q = out[col].quantile([.25, .5, .75, .95]).to_dict()
    print(f'{col}: zero {100 * (out[col] == 0).mean():.1f}%  mean {out[col].mean():.1f}  '
          f'quartiles {q[.25]:.0f}/{q[.5]:.0f}/{q[.75]:.0f}  p95 {q[.95]:.0f}')
print(f'candidates of charge 1: {100 * out.n_cand_z1.sum() / out.n_cand.sum():.1f}% of all candidates')
print(f'true precursor among candidates: {100 * out.true_in.mean():.1f}% of spectra')
print(f'true precursor among candidates within 1/K0 0.05: {100 * out.true_in_im05.mean():.1f}%')
