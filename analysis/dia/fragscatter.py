#!/usr/bin/env python3
"""
Fragment m/z scatter between two MS2 scans of the same isolation window.

`ppmSigma` in MaRaCluster's ppm merge is defined as the scatter between two
spectra of one cluster, and its default of 2 ppm is an Orbitrap number.  On the
TripleTOF data used here it has to be measured, or the merge either splits one
fragment into several peaks or joins fragments it should not.

Consecutive scans of one window are the same mixture measured again a cycle
later, so the strongest peaks of one should reappear in the other.  Each strong
peak of a scan is matched to the nearest peak of the next scan of the same
window within MAX_PPM, and the signed difference is recorded.  A robust width
is reported rather than a standard deviation, because the tail is mismatches
rather than measurement error.

    python3 fragscatter.py <run> [pairs]
"""
import sys

import numpy as np
from pyteomics import mzxml

RUN = sys.argv[1] if len(sys.argv) > 1 else '18484_REP3_1ug_Ecoli_NewStock2_SWATH_1'
PAIRS = int(sys.argv[2]) if len(sys.argv) > 2 else 3000
TOPN = 40
MAX_PPM = 100.0


def match(prev_mz, mz, targets):
    """Nearest peak of the next scan to each target, in ppm."""
    j = np.searchsorted(mz, targets)
    a = np.clip(j - 1, 0, mz.size - 1)
    b = np.clip(j, 0, mz.size - 1)
    da = np.abs(mz[a] - targets)
    db = np.abs(mz[b] - targets)
    k = np.where(da < db, a, b)
    return (mz[k] - targets) / targets * 1e6


def main():
    prev = {}
    diffs, mzs = [], []
    pairs = 0
    for s in mzxml.read(RUN + '.mzXML'):
        if s['msLevel'] != 2:
            continue
        c = float(s['precursorMz'][0]['precursorMz'])
        mz = np.asarray(s['m/z array'], np.float64)
        it = np.asarray(s['intensity array'], np.float64)
        if mz.size < TOPN:
            prev[c] = (mz, it)
            continue
        if c in prev:
            pmz, pit = prev[c]
            sel = np.argpartition(pit, -TOPN)[-TOPN:]
            targets = np.sort(pmz[sel])
            d = match(pmz, mz, targets)
            keep = np.abs(d) <= MAX_PPM
            diffs.append(d[keep])
            mzs.append(targets[keep])
            pairs += 1
            if pairs >= PAIRS:
                break
        prev[c] = (mz, it)

    d = np.concatenate(diffs)
    m = np.concatenate(mzs)
    print(f'{pairs} scan pairs, {d.size} matched peaks within {MAX_PPM:g} ppm')
    mad = float(np.median(np.abs(d - np.median(d))))
    print('median %+.2f ppm, robust sigma %.2f ppm (1.4826 MAD), '
          'q05 %+.2f, q95 %+.2f'
          % (np.median(d), 1.4826 * mad, np.percentile(d, 5),
             np.percentile(d, 95)))
    print('\n%10s %8s %10s %10s' % ('m/z', 'n', 'median', 'sigma'))
    edges = np.percentile(m, [0, 20, 40, 60, 80, 100])
    for lo, hi in zip(edges[:-1], edges[1:]):
        k = (m >= lo) & (m < hi)
        if k.sum() < 100:
            continue
        s_mad = float(np.median(np.abs(d[k] - np.median(d[k]))))
        print('%4.0f-%-5.0f %8d %+10.2f %10.2f'
              % (lo, hi, k.sum(), np.median(d[k]), 1.4826 * s_mad))
    print('\nA sigma for --consensus-ppm-tol or diatest.py is the robust sigma '
          'above,\nnot the half width: the merge cuts at four of them.')


if __name__ == '__main__':
    main()
