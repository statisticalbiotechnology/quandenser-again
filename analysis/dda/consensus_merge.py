#!/usr/bin/env python3
"""
Consensus spectra from the CSF DDA run under both merges, for comparison.

`range.py` and `comp_pairs.py` were run on the consensus spectra the pipeline
wrote, which used the Thomson merge. This builds the same kind of file with
either merge, so the two can be put through those scripts unchanged and the
caveat in `docs/dia-feasibility.md` can be checked rather than argued about.

The clusters here are defined by identity, not by MaRaCluster: spectra that
PEAKS assigned the same peptide and charge, across both replicates, are one
cluster. That isolates the merge from the clustering, which is what the caveat
is about, and it needs only the deposited MGF and mzIdentML rather than a run
of the pipeline. Cluster sizes are therefore an upper bound on what MaRaCluster
would produce, and no unidentified spectrum is in any cluster.

    python3 consensus_merge.py [--ppm-sigma 2.0] [--max-peaks 160] \\
                               [--min-size 2] [--out-dir .]

Writes one .ms2 per merge, and prints the cluster-size distribution.
"""
import argparse
import collections
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import consensusmerge                                          # noqa: E402

from identified import parse_mzid, read_mgf, PROTON            # noqa: E402

CALIB = '/home/user/csf/calib'


def clusters(min_size):
    """Spectra grouped by the peptide and charge PEAKS assigned them."""
    peptides, psms = parse_mzid(os.path.join(CALIB, 'peaks.mzid'))
    print(f'peptides in mzid: {len(peptides)}')
    bykey = collections.defaultdict(list)
    for run, byidx in psms.items():
        mgf = read_mgf(os.path.join(CALIB, f'{run}.mgf'))
        print(f'{run}: {len(mgf)} MGF spectra, {len(byidx)} rank-1 PSMs')
        for idx, (pref, z, calcmz, score) in byidx.items():
            if idx not in mgf:
                continue
            _, mz, it = mgf[idx]
            o = np.argsort(mz)
            bykey[(pref, z)].append((mz[o], it[o], calcmz))

    sizes = np.array([len(v) for v in bykey.values()])
    print(f'{len(bykey)} peptide/charge groups over both runs, '
          f'{sizes.sum()} spectra')
    print('  cluster size  1: %d, 2: %d, 3-5: %d, 6-10: %d, >10: %d'
          % ((sizes == 1).sum(), (sizes == 2).sum(),
             ((sizes >= 3) & (sizes <= 5)).sum(),
             ((sizes >= 6) & (sizes <= 10)).sum(), (sizes > 10).sum()))
    out = {k: v for k, v in bykey.items() if len(v) >= min_size}
    print(f'{len(out)} clusters of at least {min_size} spectra, '
          f'{sum(len(v) for v in out.values())} spectra in them\n')
    return out


def write_ms2(path, spectra):
    with open(path, 'w') as fh:
        fh.write('H\tCreator\tconsensus_merge.py\n')
        for scan, (z, mh, mz, it) in enumerate(spectra, 1):
            fh.write('S\t%d\t%d\t%.5f\n' % (scan, scan, (mh + (z - 1) * PROTON) / z))
            fh.write('Z\t%d\t%.5f\n' % (z, mh))
            for a, b in zip(mz, it):
                fh.write('%.5f %.4f\n' % (a, b))
    print(f'wrote {path}: {len(spectra)} consensus spectra')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--ppm-sigma', type=float, default=2.0)
    ap.add_argument('--num-sigmas', type=float, default=4.0)
    ap.add_argument('--max-peaks', type=int, default=160,
                    help='0 keeps every peak, as --consensus-max-peaks 0 does')
    ap.add_argument('--min-size', type=int, default=2)
    ap.add_argument('--out-dir', default='.')
    args = ap.parse_args()

    groups = clusters(args.min_size)
    thomson, ppm = [], []
    for (pref, z), members in sorted(groups.items()):
        cluster = [(mz, it) for mz, it, _ in members]
        calcmz = members[0][2]
        mh = calcmz * z - (z - 1) * PROTON        # M + proton, as the Z line
        mz_t, it_t = consensusmerge.merge_thomson(cluster)
        mz_p, it_p = consensusmerge.merge_ppm(cluster, args.ppm_sigma,
                                              args.num_sigmas, args.max_peaks)
        thomson.append((z, mh, mz_t, it_t))
        ppm.append((z, mh, mz_p, it_p))

    write_ms2(os.path.join(args.out_dir, 'consensus_thomson.ms2'), thomson)
    write_ms2(os.path.join(args.out_dir, 'consensus_ppm%g_max%d.ms2'
                           % (args.ppm_sigma, args.max_peaks)), ppm)

    n_t = np.array([len(s[2]) for s in thomson])
    n_p = np.array([len(s[2]) for s in ppm])
    print('\npeaks per consensus spectrum: thomson median %d, ppm median %d'
          % (np.median(n_t), np.median(n_p)))


if __name__ == '__main__':
    main()
