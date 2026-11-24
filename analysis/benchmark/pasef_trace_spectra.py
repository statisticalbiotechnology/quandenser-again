"""Follow the MS2 spectra of one run through Quandenser's clustering and
consensus spectra to the identifications.

    pasef_trace_spectra.py <clusters.tsv> <consensus_hdr.tsv> <psms.txt> <run.mzML name>
                          <ms2_header.tsv> <candidates.tsv> [raw_psms.txt]

Consensus scan = cluster number (third column of the cluster file) * 100 +
mass-charge state index. consensus_hdr.tsv (scan, precursor m/z, charge)
comes from the consensus .ms2 files:
    awk '/^S/{s=$2; m=$4} /^Z/{print s"\t"m"\t"$2}' consensus_spectra/*.ms2
"""
import sys

import numpy as np
import pandas as pd

cl_fn, cons_fn, psm_fn, run, hdr_fn, cand_fn = sys.argv[1:7]
raw_fn = sys.argv[7] if len(sys.argv) > 7 else None

cl = pd.read_csv(cl_fn, sep='\t', header=None, names=['file', 'scan', 'cluster'],
                 skip_blank_lines=True).dropna()
cl['cluster'] = cl.cluster.astype(int)
size = cl.groupby('cluster').size().rename('size')
cons = pd.read_csv(cons_fn, sep='\t', header=None, names=['cscan', 'cmz', 'cz'])
cons['cluster'] = cons.cscan // 100
psm = pd.read_csv(psm_fn, sep='\t', usecols=['scan', 'percolator q-value', 'sequence'])
psm = psm[psm['percolator q-value'] <= 0.01].drop_duplicates('scan')
cons['ided'] = cons.cscan.isin(psm.scan)

hdr = pd.read_csv(hdr_fn, sep='\t')
cand = pd.read_csv(cand_fn, sep='\t')
me = cl[cl.file == run][['scan', 'cluster']].merge(hdr, on='scan').merge(
    cand[['scan', 'n_cand', 'true_in']], on='scan').merge(size, left_on='cluster', right_index=True)
n = len(me)

per_cl = cons.groupby('cluster').agg(n_cons=('cscan', 'size'), any_ided=('ided', 'any'))
me = me.merge(per_cl, left_on='cluster', right_index=True, how='left').fillna(
    {'n_cons': 0, 'any_ided': False})
# consensus spectrum of the spectrum's cluster carrying the instrument precursor
m = me[['scan', 'cluster', 'precmz', 'charge']].merge(cons, on='cluster')
m['true'] = (m.cz == m.charge) & (np.abs(m.cmz - m.precmz) <= 30e-6 * m.precmz)
m['true_ided'] = m.true & m.ided
tr = m.groupby('scan').agg(true_cons=('true', 'any'), true_ided=('true_ided', 'any'))
me = me.merge(tr, left_on='scan', right_index=True, how='left').fillna(
    {'true_cons': False, 'true_ided': False})


def pct(x):
    return f'{100 * x.mean():.1f}%'


print(f'{run}: {n} MS2 spectra')
print(f'  in singleton clusters: {pct(me["size"] == 1)}')
print(f'  cluster has no consensus spectrum: {pct(me.n_cons == 0)}'
      f'  (spectrum has no candidate: {pct(me.n_cand == 0)})')
print(f'  consensus spectra per cluster, mean over spectra: {me.n_cons.mean():.1f}')
print(f'  instrument precursor among the spectrum\'s candidates: {pct(me.true_in == 1)}')
print(f'  a consensus spectrum of its cluster carries the instrument precursor: {pct(me.true_cons)}')
print(f'  cluster identified (any consensus spectrum, q<=0.01): {pct(me.any_ided)}')
print(f'  ... through the consensus spectrum with the instrument precursor: {pct(me.true_ided)}')
z = me.groupby('charge').agg(n=('scan', 'size'), ided=('any_ided', 'mean'))
print('  identified by instrument charge:',
      ', '.join(f'z={k} {100 * r.ided:.1f}% of {r.n}' for k, r in z.iterrows()))
print(f'consensus spectra: {len(cons)}, identified {cons.ided.sum()} ({pct(cons.ided)})')

if raw_fn:
    raw = pd.read_csv(raw_fn, sep='\t', usecols=['scan', 'percolator q-value', 'sequence'])
    raw = raw[raw['percolator q-value'] <= 0.01].drop_duplicates('scan')
    me['raw_ided'] = me.scan.isin(raw.scan)
    print(f'raw search with instrument precursors: {len(raw)} spectra at q<=0.01 ({pct(me.raw_ided)})')
    r = me[me.raw_ided]
    print(f'  of these: instrument precursor among candidates {pct(r.true_in == 1)}, '
          f'consensus with it exists {pct(r.true_cons)}, cluster identified {pct(r.any_ided)}, '
          f'via that consensus {pct(r.true_ided)}')
    print(f'  identified through Quandenser but not raw: {(me.any_ided & ~me.raw_ided).sum()}')
