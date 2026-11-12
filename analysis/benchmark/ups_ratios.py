#!/usr/bin/env python3
"""
Quantitative accuracy of a Quandenser run on PXD002370 (UPS1 in yeast).

UPS1 is spiked at 25, 10 and 5 fmol (conditions C, D, E; three runs each) into
a constant yeast background, so the expected log2 ratios are known:

    UPS1   C/E = log2(5) = 2.32   D/E = 1.00   C/D = 1.32
    yeast  all ratios 0

Identifications come from Percolator on the consensus spectra at 1% FDR. A
consensus spectrum's scan number is the spectrumId Quandenser printed in
Quandenser.feature_groups.tsv, which ties a peptide to a feature group and so
to one MS1 intensity per run. Per feature group, intensities are averaged per
condition in log2 space over the runs where the group has a feature, and a
ratio is reported when both conditions have at least two runs.

The median ratio measures accuracy; the robust sd (1.4826 MAD) of the yeast
ratios, whose truth is 0, measures precision.

    ups_ratios.py <results_dir> <percolator_dir>
"""
import collections
import csv
import math
import statistics as st
import sys

RES, PERC = sys.argv[1], sys.argv[2]
EXPECT = {('C', 'E'): math.log2(5), ('D', 'E'): 1.0, ('C', 'D'): math.log2(2.5)}

batch = [l.split('\t')[0] for l in open(f'{RES}/pipeline_info/batch.tsv') if l.strip()]
cond = {i: name[0] for i, name in enumerate(batch)}            # C/D/E from the file name

# consensus scan -> (peptide, UPS1 or yeast) at 1% FDR
ident = {}
for r in csv.DictReader(open(f'{PERC}/percolator.target.psms.txt'), delimiter='\t'):
    if float(r['percolator q-value']) > 0.01:
        continue
    prots = [p for p in r['protein id'].split(',') if p]
    ups = bool(prots) and all('ups' in p.lower() for p in prots)
    yeast = bool(prots) and not any('ups' in p.lower() for p in prots)
    if ups or yeast:
        ident[int(r['scan'])] = (r['sequence'], 'ups' if ups else 'yeast')

# feature groups: intensity per run, and the consensus scans they link to
groups, cur = [], None
for line in open(f'{RES}/quandenser/Quandenser.feature_groups.tsv'):
    line = line.rstrip('\n')
    if not line:
        cur = None
        continue
    if cur is None:
        cur = {'int': collections.defaultdict(float), 'scans': set()}
        groups.append(cur)
    f = line.split('\t')
    cur['int'][int(f[0])] += float(f[4])
    if len(f) > 5 and f[5]:
        for pair in f[5].split(','):
            spec, _ = pair.split(';')
            if spec != '0':
                cur['scans'].add(int(spec))

# peptide -> its best-covered feature group
best = {}
for g in groups:
    hits = {ident[s] for s in g['scans'] if s in ident}
    if len(hits) != 1:
        continue                                            # unidentified or ambiguous
    (pep, kind), = hits
    n = sum(1 for v in g['int'].values() if v > 0)
    if pep not in best or n > best[pep][2]:
        best[pep] = (kind, g['int'], n)

print(RES)
print(f'  identified feature groups: {len(best)} peptides '
      f'({sum(1 for v in best.values() if v[0] == "ups")} UPS1, '
      f'{sum(1 for v in best.values() if v[0] == "yeast")} yeast)')
for (a, b), exp in EXPECT.items():
    for kind, target in (('ups', exp), ('yeast', 0.0)):
        ratios = []
        for k, ints, _ in best.values():
            if k != kind:
                continue
            la = [math.log2(v) for i, v in ints.items() if cond[i] == a and v > 0]
            lb = [math.log2(v) for i, v in ints.items() if cond[i] == b and v > 0]
            if len(la) >= 2 and len(lb) >= 2:
                ratios.append(st.mean(la) - st.mean(lb))
        if not ratios:
            continue
        med = st.median(ratios)
        mad = 1.4826 * st.median(abs(x - med) for x in ratios)
        print(f'  {kind:5s} {a}/{b}: n={len(ratios):5d}  median log2 {med:+.2f} '
              f'(expected {target:+.2f}, error {med - target:+.2f})  robust sd {mad:.2f}')
