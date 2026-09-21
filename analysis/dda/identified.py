#!/usr/bin/env python3
"""
Definitive version of the calibration, restricted to spectra that PEAKS
identified in the original study, so the true peptide sequence is known.

Answers three things:
  1. How many complementary b/y pairs are actually observable in real
     high-resolution HCD spectra (the ground truth the statistic is chasing).
  2. Whether the target-vs-decoy statistic detects them without knowing the
     sequence.
  3. How that depends on precursor charge, which is what decides whether the
     statistic is usable on DIA data.
"""
import sys, re, collections
import xml.etree.ElementTree as ET
import numpy as np

PROTON = 1.00727646
H2O = 18.0105646
AA = dict(G=57.02146, A=71.03711, S=87.03203, P=97.05276, V=99.06841,
          T=101.04768, C=103.00919, L=113.08406, I=113.08406, N=114.04293,
          D=115.02694, Q=128.05858, K=128.09496, E=129.04259, M=131.04049,
          H=137.05891, F=147.06841, R=156.10111, Y=163.06333, W=186.07931)
SHIFTS = np.array([-43,-41,-37,-31,-29,-23,-19,-13,-9,-6,
                     6,  9, 13, 19, 23, 29, 31, 37, 41, 43], float)
TOL = 0.01
NS = '{http://psidev.info/psi/pi/mzIdentML/1.1}'
WANT = {'SPECTRADATA_21': '20150513_20_hTau_CSF_Frac9of12_rep1',
        'SPECTRADATA_33': '20150513_48_hTau_CSF_Frac9of12_rep2'}


def parse_mzid(path):
    peptides, psms = {}, collections.defaultdict(dict)
    for ev, el in ET.iterparse(path, events=('end',)):
        if el.tag == NS + 'Peptide':
            s = el.findtext(NS + 'PeptideSequence')
            mods = [(int(m.get('location', 0)), float(m.get('monoisotopicMassDelta', 0)))
                    for m in el.findall(NS + 'Modification')]
            peptides[el.get('id')] = (s, mods)
            el.clear()
        elif el.tag == NS + 'SpectrumIdentificationResult':
            sd = el.get('spectraData_ref')
            if sd in WANT:
                sii = el.find(NS + 'SpectrumIdentificationItem')
                if sii is not None and sii.get('rank') == '1':
                    idx = int(el.get('spectrumID').split('=')[1])
                    sc = sii.find(NS + 'cvParam')
                    psms[WANT[sd]][idx] = (
                        sii.get('peptide_ref'), int(sii.get('chargeState')),
                        float(sii.get('calculatedMassToCharge')),
                        float(sc.get('value')) if sc is not None else 0.0)
            el.clear()
    return peptides, psms


def read_mgf(path):
    out, idx, z, mz, it, inblock = {}, None, None, [], [], False
    for line in open(path):
        line = line.strip()
        if line == 'BEGIN IONS':
            idx, z, mz, it, inblock = None, None, [], [], True
        elif line == 'END IONS':
            if idx is not None and mz:
                out[idx] = (z, np.array(mz), np.array(it))
            inblock = False
        elif not inblock:
            continue
        elif line.startswith('TITLE='):
            idx = int(line.split('index=')[1])
        elif line.startswith('CHARGE='):
            z = int(line[7:].rstrip('+'))
        elif line and line[0].isdigit():
            a = line.split()
            mz.append(float(a[0])); it.append(float(a[1]))
    return out


def by_ions(seq, mods):
    res = [AA.get(c, 0.0) for c in seq]
    for loc, d in mods:
        if 1 <= loc <= len(res):
            res[loc - 1] += d
        elif loc == 0:
            res[0] += d
        else:
            res[-1] += d
    cum = np.cumsum(res)
    tot = cum[-1]
    b = cum[:-1] + PROTON                      # b1..b(n-1)
    y = (tot - cum[:-1]) + H2O + PROTON        # y(n-1)..y1
    return b, y, tot + H2O                     # neutral peptide mass


def count_11(mz, target, tol):
    lo = np.searchsorted(mz, target - tol - mz, 'left')
    hi = np.searchsorted(mz, target + tol - mz, 'right')
    tot = int((hi - lo).sum())
    self_hits = int(np.sum(np.abs(2 * mz - target) <= tol))
    return (tot - self_hits) // 2


def count_12(mz, target, tol):
    twice = 2.0 * mz
    lo = np.searchsorted(mz, target - tol - twice, 'left')
    hi = np.searchsorted(mz, target + tol - twice, 'right')
    tot = int((hi - lo).sum())
    self_hits = int(np.sum(np.abs(3 * mz - target) <= tol))
    return tot - self_hits


def main():
    peptides, psms = parse_mzid('/home/user/csf/calib/peaks.mzid')
    print(f'peptides in mzid: {len(peptides)}')
    rows = []
    for run, byidx in psms.items():
        mgf = read_mgf(f'/home/user/csf/calib/{run}.mgf')
        print(f'{run}: {len(mgf)} MGF spectra, {len(byidx)} rank-1 PSMs')
        for idx, (pref, z, calcmz, score) in byidx.items():
            if idx not in mgf:
                continue
            seq, mods = peptides[pref]
            _, mzs, ints = mgf[idx]
            order = np.argsort(mzs)
            mzs, ints = mzs[order], ints[order]
            b, y, Mseq = by_ions(seq, mods)
            M = calcmz * z - z * PROTON
            # ground truth: which b/y are observed, and which pairs are complete
            obs_b = np.abs(mzs[np.searchsorted(mzs, b).clip(0, len(mzs) - 1)] - b) <= TOL
            j = np.searchsorted(mzs, b)
            obs_b = np.zeros(len(b), bool)
            for k, v in enumerate(b):
                p = np.searchsorted(mzs, v)
                obs_b[k] = ((p < len(mzs) and abs(mzs[p] - v) <= TOL) or
                            (p > 0 and abs(mzs[p - 1] - v) <= TOL))
            obs_y = np.zeros(len(y), bool)
            for k, v in enumerate(y):
                p = np.searchsorted(mzs, v)
                obs_y[k] = ((p < len(mzs) and abs(mzs[p] - v) <= TOL) or
                            (p > 0 and abs(mzs[p - 1] - v) <= TOL))
            true_pairs = int(np.sum(obs_b & obs_y))     # b_i and y_{n-i} both seen
            s0 = M + 2 * PROTON
            t11 = count_11(mzs, s0, TOL)
            d11 = [count_11(mzs, s0 + d, TOL) for d in SHIFTS]
            s1 = M + 3 * PROTON
            t12 = count_12(mzs, s1, TOL)
            d12 = [count_12(mzs, s1 + d, TOL) for d in SHIFTS]
            rows.append((z, len(seq), score, len(mzs), int(obs_b.sum()),
                         int(obs_y.sum()), true_pairs,
                         t11, int(np.sum(np.array(d11) >= t11)), float(np.mean(d11)),
                         t12, int(np.sum(np.array(d12) >= t12)), float(np.mean(d12))))
    a = np.array(rows, float)
    np.save('/home/user/csf/calib/identified.npy', a)
    print(f'\nidentified spectra analysed: {len(a)}')
    hdr = ('z','n','len','peaks','b obs','y obs','true by pairs','found (1,1)',
           'decoy mean','beats all','+ (1,2)','beats all')
    z, plen, score, npk = a[:,0], a[:,1], a[:,2], a[:,3]
    nb, ny, tp = a[:,4], a[:,5], a[:,6]
    t11, be11, dm11 = a[:,7], a[:,8], a[:,9]
    t12, be12, dm12 = a[:,10], a[:,11], a[:,12]
    top11 = be11 == 0
    comb_t = t11 + t12
    comb_beats = ((a[:,8] + a[:,11]) == 0)
    print()
    print(f'{"z":<4}{"n":>7}{"pep len":>9}{"peaks":>8}{"b seen":>8}{"y seen":>8}'
          f'{"true b/y pairs":>16}{"found(1,1)":>12}{"decoy":>8}{"beats-all":>11}')
    for c in (1,2,3,4,5,6):
        m = z == c
        if m.sum() < 20: continue
        print(f'{c:<4}{m.sum():>7}{plen[m].mean():>9.1f}{npk[m].mean():>8.0f}'
              f'{nb[m].mean():>8.1f}{ny[m].mean():>8.1f}{tp[m].mean():>16.2f}'
              f'{t11[m].mean():>12.2f}{dm11[m].mean():>8.2f}'
              f'{100*top11[m].mean():>10.1f}%')
    m = np.ones(len(a), bool)
    print(f'{"all":<4}{m.sum():>7}{plen.mean():>9.1f}{npk.mean():>8.0f}'
          f'{nb.mean():>8.1f}{ny.mean():>8.1f}{tp.mean():>16.2f}'
          f'{t11.mean():>12.2f}{dm11.mean():>8.2f}{100*top11.mean():>10.1f}%')
    print(f'\nnull beats-all rate = {100/21:.2f}%')
    print(f'adding the (1,2) relation, beats-all on the combined count: '
          f'{100*comb_beats.mean():.1f}%')
    print(f'\nfraction of identified spectra with at least one TRUE b/y '
          f'complementary pair: {100*np.mean(tp>0):.1f}%')
    for k in (1,2,3,5):
        print(f'  >= {k} true pairs: {100*np.mean(tp>=k):5.1f}%'
              f'   (beats-all among these: {100*top11[tp>=k].mean():5.1f}%)')


if __name__ == '__main__':
    main()
