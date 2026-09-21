#!/usr/bin/env python3
"""
Cross-run pooling of channel 1 on real DIA data.

Peptides confident in both SWATH replicates have their fragment traces and
their MS1 precursor traces concatenated across runs before correlating.  Each
run contributes its own retention time from its own OpenSWATH peak group, so
the alignment is exact and the measurement isolates the pooling gain rather
than testing an aligner.
"""
import pickle
import numpy as np

R1 = 'traces_18484_REP3_1ug_Ecoli_NewStock2_SWATH_1.pkl'
R2 = 'traces_18486_REP3_1ug_Ecoli_NewStock2_SWATH_2.pkl'


def pearson(a, b):
    a = a - a.mean(); b = b - b.mean()
    da = np.sqrt((a * a).sum()); db = np.sqrt((b * b).sum())
    return np.nan if da <= 0 or db <= 0 else float((a * b).sum() / (da * db))


def auc(pos, neg):
    pos = np.asarray(pos); neg = np.asarray(neg)
    pos = pos[~np.isnan(pos)]; neg = neg[~np.isnan(neg)]
    if not len(pos) or not len(neg):
        return np.nan, np.nan, np.nan
    allv = np.concatenate([pos, neg])
    o = np.argsort(allv, kind='mergesort'); sv = allv[o]
    rr = np.empty(len(allv)); i = 0
    while i < len(sv):
        j = i
        while j + 1 < len(sv) and sv[j + 1] == sv[i]:
            j += 1
        rr[o[i:j + 1]] = (i + j) / 2.0 + 1.0
        i = j + 1
    a = (rr[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2.0) / (len(pos) * len(neg))
    return a, np.mean(pos >= np.percentile(neg, 95)), np.mean(pos >= np.percentile(neg, 99))


def window(p):
    ms1, T, D = p['ms1'], p['tr_true'], p['tr_decoy']
    if ms1.size < 5 or T.size == 0:
        return None
    n = min(len(ms1), len(T))
    if n < 5:
        return None
    w = p['hi'] - p['lo']; half = max(1.5 * w, 15.0)
    keep = (p['rts'][:n] >= p['rt'] - half) & (p['rts'][:n] <= p['rt'] + half)
    if keep.sum() < 5:
        return None
    Dk = D[:n][keep] if D.size and len(D) >= n else None
    return ms1[:n][keep].astype(float), T[:n][keep], Dk


a = {p['key']: p for p in pickle.load(open(R1, 'rb'))}
b = {p['key']: p for p in pickle.load(open(R2, 'rb'))}
shared = sorted(set(a) & set(b))
print(f'run1 {len(a)} peptides, run2 {len(b)}, confident in both: {len(shared)}\n')

res = {'one': ([], []), 'two': ([], [])}
used = 0
for k in shared:
    wa, wb = window(a[k]), window(b[k])
    if wa is None or wb is None:
        continue
    # true fragment lists are the same library entry, so columns correspond
    if wa[1].shape[1] != wb[1].shape[1]:
        continue
    used += 1
    for tag in ('one', 'two'):
        if tag == 'one':
            ref, T = wa[0], wa[1]
            D = wa[2]
        else:
            n = min(len(wa[0]), len(wb[0]))
            ref = np.concatenate([wa[0][:n], wb[0][:n]])
            T = np.concatenate([wa[1][:n], wb[1][:n]])
            D = (np.concatenate([wa[2][:n], wb[2][:n]])
                 if wa[2] is not None and wb[2] is not None
                 and wa[2].shape[1] == wb[2].shape[1] else None)
        for j in range(T.shape[1]):
            r = pearson(ref, T[:, j].astype(float))
            if not np.isnan(r):
                res[tag][0].append(r)
        if D is not None:
            for j in range(D.shape[1]):
                r = pearson(ref, D[:, j].astype(float))
                if not np.isnan(r):
                    res[tag][1].append(r)

print(f'peptides usable in both runs: {used}\n')
print(f'{"":<14}{"n true":>9}{"n decoy":>9}{"med r true":>12}{"med r decoy":>13}'
      f'{"AUC":>8}{"TPR@5%":>9}{"TPR@1%":>9}')
for tag, label in (('one', 'single run'), ('two', 'two runs pooled')):
    p, n = res[tag]
    a_, t5, t1 = auc(p, n)
    print(f'{label:<14}{len(p):>9}{len(n):>9}{np.median(p):>12.3f}'
          f'{np.median(n):>13.3f}{a_:>8.3f}{100*t5:>8.1f}%{100*t1:>8.1f}%')
