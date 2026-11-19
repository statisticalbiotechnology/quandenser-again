#!/usr/bin/env python3
"""
Quantitative accuracy of a Quandenser run on a human/yeast/E. coli (HYE)
benchmark such as PXD028735.

Condition A is 65% human, 30% yeast and 5% E. coli; condition B is 65%
human, 15% yeast and 20% E. coli. The expected log2 A/B ratios are therefore

    human 0     yeast +1     E. coli -2

Identifications come from Percolator on the consensus spectra at 1% FDR; a
peptide counts for one species only if every protein it maps to is from that
species. As in ups_ratios.py, a consensus spectrum's scan number is the
spectrumId in Quandenser.feature_groups.tsv, which ties the peptide to a
feature group and one MS1 intensity per run. Per feature group, intensities
are averaged per condition in log2 space over the runs where the group has a
feature, and a ratio is reported when both conditions have at least two runs.
The condition is the first letter of the file name.

The median ratio measures accuracy; the robust sd (1.4826 MAD) measures
precision.

    hye_ratios.py <results_dir> <percolator_dir>
"""
import collections
import csv
import math
import statistics as st
import sys

RES, PERC = sys.argv[1], sys.argv[2]
EXPECT = {'HUMAN': 0.0, 'YEAST': 1.0, 'ECOLI': -2.0}

batch = [l.split('\t')[0] for l in open(f'{RES}/pipeline_info/batch.tsv') if l.strip()]
cond = {i: name.split('/')[-1][0] for i, name in enumerate(batch)}


def species(protein):
    tail = protein.split()[0].rsplit('_', 1)
    return tail[1] if len(tail) == 2 and tail[1] in EXPECT else None


# consensus scan -> (peptide, species) at 1% FDR
ident = {}
for r in csv.DictReader(open(f'{PERC}/percolator.target.psms.txt'), delimiter='\t'):
    if float(r['percolator q-value']) > 0.01:
        continue
    kinds = {species(p) for p in r['protein id'].split(',') if p}
    if len(kinds) == 1 and None not in kinds:
        ident[int(r['scan'])] = (r['sequence'], kinds.pop())

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
counts = collections.Counter(v[0] for v in best.values())
print(f'  identified feature groups: {len(best)} peptides ('
      + ', '.join(f'{counts[k]} {k.lower()}' for k in EXPECT) + ')')
complete = sum(1 for v in best.values() if v[2] == len(batch))
print(f'  quantified in all {len(batch)} runs: {complete} peptides')
for kind, target in EXPECT.items():
    ratios = []
    for k, ints, _ in best.values():
        if k != kind:
            continue
        la = [math.log2(v) for i, v in ints.items() if cond[i] == 'A' and v > 0]
        lb = [math.log2(v) for i, v in ints.items() if cond[i] == 'B' and v > 0]
        if len(la) >= 2 and len(lb) >= 2:
            ratios.append(st.mean(la) - st.mean(lb))
    if not ratios:
        continue
    med = st.median(ratios)
    mad = 1.4826 * st.median(abs(x - med) for x in ratios)
    print(f'  {kind.lower():5s} A/B: n={len(ratios):5d}  median log2 {med:+.2f} '
          f'(expected {target:+.2f}, error {med - target:+.2f})  robust sd {mad:.2f}')
