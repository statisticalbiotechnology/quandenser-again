#!/usr/bin/env python3
"""
OpenSWATH's own sub-scores on a real DIA run, with target/decoy labels.

var_xcorr_shape / var_xcorr_coelution are channel 1: how well the fragment
traces co-elute with each other.  var_bseries_score / var_yseries_score and
var_library_* are fragment-identity evidence, the family channel 2 belongs to.
Decoys are shuffled-sequence transition groups, i.e. fragment sets that
correspond to no real peptide, which is the right null.
"""
import numpy as np

F = 'split_18484_REP3_1ug_Ecoli_NewStock2_SWATH_1_mprophet_all_peakgroups_DDA.xls'
rows = open(F).read().rstrip('\n').split('\n')
hdr = rows[0].split('\t')
ix = {k: i for i, k in enumerate(hdr)}
data = [r.split('\t') for r in rows[1:]]
print(f'{len(data)} peak groups')


def col(name, cast=float):
    out = []
    for r in data:
        try:
            out.append(cast(r[ix[name]]))
        except (ValueError, IndexError):
            out.append(np.nan)
    return np.array(out)


decoy = np.array([r[ix['decoy']].strip().upper() == 'TRUE' for r in data])
print(f'targets {int((~decoy).sum())}, decoys {int(decoy.sum())}\n')


def auc(s, lab):
    m = ~np.isnan(s)
    s, lab = s[m], lab[m]
    p, n = s[lab], s[~lab]
    if not len(p) or not len(n):
        return np.nan
    o = np.argsort(s, kind='mergesort'); sv = s[o]
    rr = np.empty(len(s)); i = 0
    while i < len(sv):
        j = i
        while j + 1 < len(sv) and sv[j + 1] == sv[i]:
            j += 1
        rr[o[i:j + 1]] = (i + j) / 2.0 + 1.0
        i = j + 1
    rp = rr[lab].sum()
    return (rp - len(p) * (len(p) + 1) / 2.0) / (len(p) * len(n))


feats = [c for c in hdr if c.startswith('var_')]
res = []
for c in feats:
    s = col(c)
    a = auc(s, ~decoy)
    if not np.isnan(a):
        res.append((max(a, 1 - a), a, c))
res.sort(reverse=True)
print(f'{"feature":<34}{"AUC":>8}   (target vs decoy, direction-corrected)')
for absa, a, c in res:
    fam = ('co-elution' if 'xcorr' in c or 'elution' in c else
           'fragment identity' if 'series' in c or 'library' in c or 'dotprod' in c else
           'other')
    print(f'  {c:<32}{absa:>8.3f}   {fam}')
