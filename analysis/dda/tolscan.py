#!/usr/bin/env python3
"""
Background pairs scale as n^2 * tol; true pairs do not scale with tol at all.
So the usable multiplexing depth should be set by MS2 mass accuracy.  Measured
here by rerunning the chimera experiment across tolerances.
"""
import runpy, sys, random, collections
import numpy as np
m = runpy.run_path('/home/user/csf/calib/chimera.py', run_name='_lib')
parse_mzid, read_mgf, true_pairs, count_11 = (
    m['parse_mzid'], m['read_mgf'], m['true_pairs'], m['count_11'])
PROTON = m['PROTON']
SHIFTS = m['SHIFTS']; NDEC = len(SHIFTS)
WANT = m['WANT']

TOLS = [0.002, 0.005, 0.01, 0.02]
KS = [1, 4, 8, 16]

rng = random.Random(1)
pep, psms = parse_mzid('/home/user/csf/calib/peaks.mzid')
targets, pools = [], {}
for run, byidx in psms.items():
    mgf = read_mgf(f'/home/user/csf/calib/{run}.mgf')
    pools[run] = [np.sort(x[1]) for x in mgf.values()]
    for idx, (pref, z, calcmz) in byidx.items():
        if idx in mgf:
            targets.append((run, z, calcmz * z - z * PROTON, np.sort(mgf[idx][1])))
print(f'{len(targets)} target spectra, {NDEC} decoys, '
      f'null beats-all {100/(NDEC+1):.2f}%\n')

# fix the mixtures once so tolerances are compared on identical spectra
mixes = {}
for k in KS:
    r = random.Random(1)
    mixes[k] = [np.sort(np.concatenate(
        [t[3]] + [pools[t[0]][r.randrange(len(pools[t[0]]))] for _ in range(k - 1)]))
        if k > 1 else t[3] for t in targets]

print(f'{"beats-all %":<14}' + ''.join(f'{f"+/-{t}Da":>12}' for t in TOLS))
for k in KS:
    row = []
    for tol in TOLS:
        top = 0
        for i, (run, z, M, _) in enumerate(targets):
            mix = mixes[k][i]
            s0 = M + 2 * PROTON
            t = count_11(mix, s0, tol)
            d = np.fromiter((count_11(mix, s0 + sh, tol) for sh in SHIFTS),
                            int, NDEC)
            if (d >= t).sum() == 0:
                top += 1
        row.append(100 * top / len(targets))
    print(f'{f"k={k:<3}":<14}' + ''.join(f'{v:>11.1f}%' for v in row))
