#!/usr/bin/env python3
"""
Channel 1 test: can a chromatographic profile separate co-eluting species?

Fragment-level co-variation cannot be tested on DDA data, because a DDA run
fragments each precursor once rather than across its elution.  But the binding
constraint on channel 1 is not fragment extraction, it is whether elution
profiles of distinct co-eluting species are distinguishable at all at DIA
sampling rates.  That is measurable on MS1, where the ground truth is free:

  positive pair  = monoisotopic XIC and M+1 XIC of the same Dinosaur feature.
                   Same molecule by construction, different m/z, so each
                   experiences its own interference.  This is the best
                   available stand-in for two fragments of one peptide.
  negative pair  = monoisotopic XICs of two DIFFERENT features that co-elute
                   and fall in the same simulated isolation window.  This is
                   exactly the interference a DIA deconvolution must reject.

Both are measured on the same scans, so nothing about the comparison depends
on a model of chromatography.
"""
import sys, pickle
import numpy as np
from pyteomics import mzml

RUN = sys.argv[1] if len(sys.argv) > 1 else '20150513_20_hTau_CSF_Frac9of12_rep1'
MZML = f'/home/user/csf/mzml/{RUN}.mzML'
FEAT = f'/home/user/csf/results/features/{RUN}.features.tsv'
RT_LO, RT_HI = 1800.0, 3600.0          # seconds; mid-gradient, highest complexity
PPM = 10e-6
C13 = 1.0033548


def load_features():
    rows = []
    with open(FEAT) as fh:
        hdr = fh.readline().rstrip('\n').split('\t')
        c = {k: i for i, k in enumerate(hdr)}
        for line in fh:
            f = line.rstrip('\n').split('\t')
            apex = float(f[c['rtApex']]) * 60.0
            if not (RT_LO <= apex <= RT_HI):
                continue
            z = int(f[c['charge']])
            if not (2 <= z <= 4):
                continue
            if int(f[c['nIsotopes']]) < 2 or float(f[c['averagineCorr']]) < 0.9:
                continue
            rows.append(dict(mz=float(f[c['mz']]), z=z, apex=apex,
                             lo=float(f[c['rtStart']]) * 60.0,
                             hi=float(f[c['rtEnd']]) * 60.0,
                             fwhm=float(f[c['fwhm']]) * 60.0,
                             inten=float(f[c['intensityApex']])))
    return rows


def main():
    feats = load_features()
    print(f'{RUN}: {len(feats)} well-formed features with apex in '
          f'[{RT_LO:.0f}, {RT_HI:.0f}] s')
    fw = np.array([f['fwhm'] for f in feats])
    print(f'  FWHM  median {np.median(fw):.1f} s   q10 {np.percentile(fw,10):.1f}'
          f'   q90 {np.percentile(fw,90):.1f}')

    # two XIC targets per feature: monoisotope and M+1
    targets = np.array([f['mz'] for f in feats] +
                       [f['mz'] + C13 / f['z'] for f in feats])
    order = np.argsort(targets)
    tsorted = targets[order]
    tol = tsorted * PPM

    cols, rts = [], []
    for s in mzml.read(MZML):
        if s.get('ms level') != 1:
            continue
        rt = float(s['scanList']['scan'][0]['scan start time']) * 60.0
        if rt < RT_LO - 30:
            continue
        if rt > RT_HI + 30:
            break
        smz = np.asarray(s['m/z array'], float)
        sit = np.asarray(s['intensity array'], float)
        if smz.size == 0:
            continue
        idx = np.searchsorted(smz, tsorted)
        lo = np.clip(idx - 1, 0, smz.size - 1)
        hi = np.clip(idx, 0, smz.size - 1)
        dl = np.abs(smz[lo] - tsorted)
        dh = np.abs(smz[hi] - tsorted)
        pick = np.where(dl < dh, lo, hi)
        d = np.minimum(dl, dh)
        cols.append(np.where(d <= tol, sit[pick], 0.0).astype(np.float32))
        rts.append(rt)

    X = np.empty((len(tsorted), len(cols)), np.float32)
    for j, c in enumerate(cols):
        X[:, j] = c
    inv = np.empty_like(order)
    inv[order] = np.arange(len(order))
    X = X[inv]                                   # back to feature order
    rts = np.array(rts)
    n = len(feats)
    print(f'  MS1 scans in window: {len(rts)}   cycle '
          f'{np.median(np.diff(rts)):.2f} s')
    print(f'  XIC matrix: {X.shape}  ({X.nbytes/1e6:.0f} MB)')
    with open(f'/home/user/csf/calib/xic_{RUN}.pkl', 'wb') as fh:
        pickle.dump(dict(feats=feats, rts=rts, mono=X[:n], iso1=X[n:]), fh, 4)
    print(f'  saved xic_{RUN}.pkl')


if __name__ == '__main__':
    main()
