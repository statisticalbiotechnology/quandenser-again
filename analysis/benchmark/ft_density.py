#!/usr/bin/env python3
"""
How often a feature in one run has a candidate partner in another, as
Quandenser's feature matcher sees it: same charge, m/z within PPM, RT apex
within RT minutes, and optionally 1/K0 within IM. Decoy candidates are the
same search at m/z + 5.0025, Quandenser's decoy offset.

When decoy candidates are nearly as common as real ones, Percolator has
nothing to learn from and match-between-runs finds no links.

    ft_density.py <query features.tsv> <target features.tsv>

Reads biosaur2's own output (ion mobility in `im`) or the Dinosaur-layout
files the pipeline writes (in `ionMobility`, when present); the ion-mobility
rows are skipped for files without it.
"""
import csv
import sys

import numpy as np

DECOY_OFFSET = 5.0 * 1.000508
TOLERANCES = [(20, 3.0, None), (20, 1.0, None), (10, 1.0, None), (10, 0.5, None),
              (20, 3.0, 0.03), (10, 1.0, 0.03), (10, 0.5, 0.02)]


def load(path):
    with open(path) as fh:
        rows = [(float(x['mz']), int(float(x['charge'])), float(x['rtApex']),
                 float(x.get('im') or x.get('ionMobility') or 'nan'))
                for x in csv.DictReader(fh, delimiter='\t')]
    a = np.array(rows)
    return a[np.argsort(a[:, 0])]


def has_candidate(query, target, offset, ppm, rt, im):
    mz = query[:, 0] + offset
    lo = np.searchsorted(target[:, 0], mz * (1 - ppm * 1e-6))
    hi = np.searchsorted(target[:, 0], mz * (1 + ppm * 1e-6))
    hit = np.zeros(len(query), dtype=bool)
    for i in range(len(query)):
        c = target[lo[i]:hi[i]]
        ok = (c[:, 1] == query[i, 1]) & (np.abs(c[:, 2] - query[i, 2]) <= rt)
        if im is not None:
            ok &= np.abs(c[:, 3] - query[i, 3]) <= im
        hit[i] = ok.any()
    return hit


def main():
    a, b = load(sys.argv[1]), load(sys.argv[2])
    rng = np.random.default_rng(1)
    query = a[rng.choice(len(a), min(40000, len(a)), replace=False)]
    print(f'{len(a)} vs {len(b)} features')
    for ppm, rt, im in TOLERANCES:
        if im is not None and np.all(np.isnan(b[:, 3])):
            continue
        t = has_candidate(query, b, 0.0, ppm, rt, im).mean()
        d = has_candidate(query, b, DECOY_OFFSET, ppm, rt, im).mean()
        print(f'  {ppm:>2} ppm, RT +-{rt} min, 1/K0 {im}: with a target {100 * t:.1f}%, '
              f'with a decoy {100 * d:.1f}%, excess {100 * (t - d):.1f}%')


if __name__ == '__main__':
    main()
