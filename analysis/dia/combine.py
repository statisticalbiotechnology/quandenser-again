#!/usr/bin/env python3
"""
Are the co-elution and fragment-identity channels independent enough that
combining them reaches a usable operating point?

Peak groups are reduced to the best-scoring one per transition group, the way
mProphet does, so that "target" means a real peptide's peak group rather than
any of the several candidate peak groups per peptide.  Cross-validation is
grouped by peptide so the same transition group never appears in both folds.
"""
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, roc_curve
from scipy.stats import spearmanr

F = 'split_18484_REP3_1ug_Ecoli_NewStock2_SWATH_1_mprophet_all_peakgroups_DDA.xls'
rows = open(F).read().rstrip('\n').split('\n')
hdr = rows[0].split('\t'); ix = {k: i for i, k in enumerate(hdr)}
data = [r.split('\t') for r in rows[1:]]

COEL = ['var_xcorr_shape', 'var_xcorr_shape_weighted', 'var_xcorr_coelution',
        'var_xcorr_coelution_weighted', 'var_elution_model_fit_score']
IDENT = ['var_library_corr', 'var_library_dotprod', 'var_library_manhattan',
         'var_library_rmsd', 'var_dotprod_score', 'var_manhatt_score',
         'var_bseries_score', 'var_yseries_score']
OTHER = ['var_intensity_score', 'var_log_sn_score', 'var_isotope_correlation_score',
         'var_isotope_overlap_score', 'var_massdev_score',
         'var_massdev_score_weighted', 'var_norm_rt_score']
ALLF = COEL + IDENT + OTHER


def num(r, c):
    try:
        return float(r[ix[c]])
    except (ValueError, IndexError):
        return np.nan


X = np.array([[num(r, c) for c in ALLF] for r in data])
y = np.array([r[ix['decoy']].strip().upper() == 'FALSE' for r in data])
grp = np.array([r[ix['transition_group_id']].rsplit('_run', 1)[0] for r in data])
main = np.array([num(r, 'main_var_xx_swath_prelim_score') for r in data])

good = ~np.isnan(X).any(1) & ~np.isnan(main)
X, y, grp, main = X[good], y[good], grp[good], main[good]

# best peak group per transition group
best = {}
for i, g in enumerate(grp):
    if g not in best or main[i] > main[best[g]]:
        best[g] = i
sel = np.array(sorted(best.values()))
X, y, grp = X[sel], y[sel], grp[sel]
print(f'{len(y)} transition groups after keeping the best peak group each: '
      f'{int(y.sum())} target, {int((~y).sum())} decoy\n')

col = {c: i for i, c in enumerate(ALLF)}


def cv_auc(cols):
    idx = [col[c] for c in cols]
    Z = X[:, idx]
    oof = np.zeros(len(y))
    for tr, te in GroupKFold(n_splits=5).split(Z, y, grp):
        sc = StandardScaler().fit(Z[tr])
        m = LogisticRegression(max_iter=2000, C=1.0).fit(sc.transform(Z[tr]), y[tr])
        oof[te] = m.decision_function(sc.transform(Z[te]))
    a = roc_auc_score(y, oof)
    fpr, tpr, _ = roc_curve(y, oof)
    return a, np.interp(0.01, fpr, tpr), np.interp(0.05, fpr, tpr)


print(f'{"feature set":<38}{"AUC":>8}{"TPR@1%FPR":>12}{"TPR@5%FPR":>12}')
for name, cols in [('co-elution only (channel 1)', COEL),
                   ('fragment identity only', IDENT),
                   ('co-elution + identity', COEL + IDENT),
                   ('everything OpenSWATH has', ALLF)]:
    a, t1, t5 = cv_auc(cols)
    print(f'  {name:<36}{a:>8.3f}{100*t1:>11.1f}%{100*t5:>11.1f}%')

print('\nrank correlation between the two channels (best feature of each):')
for lab, m in (('targets', y), ('decoys', ~y)):
    r = spearmanr(X[m, col['var_xcorr_shape']], X[m, col['var_library_manhattan']]).statistic
    r2 = spearmanr(X[m, col['var_xcorr_shape']], X[m, col['var_yseries_score']]).statistic
    print(f'  {lab:<9} xcorr_shape vs library_manhattan  rho={r:+.3f}'
          f'     vs yseries  rho={r2:+.3f}')
