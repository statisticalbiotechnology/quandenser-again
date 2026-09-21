#!/usr/bin/env python3
"""
Does the complementary-pair statistic survive multiplexing?

DDA isolates one precursor, DIA co-isolates everything in the window.  We
simulate that by pooling each identified spectrum's peaks with those of k-1
randomly drawn spectra from the same run, then asking whether the statistic
still finds the known peptide mass of the target among 110 decoy masses.

110 decoy shifts (all integers from 6 to 60 Da, both signs) give a minimum
permutation p-value of 1/111 = 0.90%, which leaves headroom for the multiple
testing that DIA imposes: in a real window you test every MS1 feature, not one.
"""
import sys, collections, random
import xml.etree.ElementTree as ET
import numpy as np

PROTON = 1.00727646
H2O = 18.0105646
AA = dict(G=57.02146, A=71.03711, S=87.03203, P=97.05276, V=99.06841,
          T=101.04768, C=103.00919, L=113.08406, I=113.08406, N=114.04293,
          D=115.02694, Q=128.05858, K=128.09496, E=129.04259, M=131.04049,
          H=137.05891, F=147.06841, R=156.10111, Y=163.06333, W=186.07931)
SHIFTS = np.array([s * d for d in range(6, 61) for s in (-1, 1)], float)
NDEC = len(SHIFTS)
TOL = 0.01
KS = [1, 2, 4, 8, 16, 32]
NS = '{http://psidev.info/psi/pi/mzIdentML/1.1}'
WANT = {'SPECTRADATA_21': '20150513_20_hTau_CSF_Frac9of12_rep1',
        'SPECTRADATA_33': '20150513_48_hTau_CSF_Frac9of12_rep2'}


def parse_mzid(path):
    pep, psms = {}, collections.defaultdict(dict)
    for _, el in ET.iterparse(path, events=('end',)):
        if el.tag == NS + 'Peptide':
            pep[el.get('id')] = (el.findtext(NS + 'PeptideSequence'),
                                 [(int(m.get('location', 0)),
                                   float(m.get('monoisotopicMassDelta', 0)))
                                  for m in el.findall(NS + 'Modification')])
            el.clear()
        elif el.tag == NS + 'SpectrumIdentificationResult':
            sd = el.get('spectraData_ref')
            if sd in WANT:
                sii = el.find(NS + 'SpectrumIdentificationItem')
                if sii is not None and sii.get('rank') == '1':
                    psms[WANT[sd]][int(el.get('spectrumID').split('=')[1])] = (
                        sii.get('peptide_ref'), int(sii.get('chargeState')),
                        float(sii.get('calculatedMassToCharge')))
            el.clear()
    return pep, psms


def read_mgf(path):
    out, idx, z, mz, inblock = {}, None, None, [], False
    for line in open(path):
        line = line.strip()
        if line == 'BEGIN IONS':
            idx, z, mz, inblock = None, None, [], True
        elif line == 'END IONS':
            if idx is not None and mz:
                out[idx] = (z, np.array(mz))
            inblock = False
        elif not inblock:
            continue
        elif line.startswith('TITLE='):
            idx = int(line.split('index=')[1])
        elif line.startswith('CHARGE='):
            z = int(line[7:].rstrip('+'))
        elif line and line[0].isdigit():
            mz.append(float(line.split(None, 1)[0]))
    return out


def true_pairs(seq, mods, mzs):
    res = [AA.get(c, 0.0) for c in seq]
    for loc, d in mods:
        res[min(max(loc - 1, 0), len(res) - 1)] += d
    cum = np.cumsum(res)
    b = cum[:-1] + PROTON
    y = (cum[-1] - cum[:-1]) + H2O + PROTON
    def seen(v):
        p = np.searchsorted(mzs, v)
        return ((p < len(mzs) and abs(mzs[p] - v) <= TOL) or
                (p > 0 and abs(mzs[p - 1] - v) <= TOL))
    return int(sum(seen(bi) and seen(yi) for bi, yi in zip(b, y)))


def count_11(mz, target, tol=TOL):
    lo = np.searchsorted(mz, target - tol - mz, 'left')
    hi = np.searchsorted(mz, target + tol - mz, 'right')
    tot = int((hi - lo).sum())
    self_hits = int(np.sum(np.abs(2 * mz - target) <= tol))
    return (tot - self_hits) // 2


def main():
    rng = random.Random(1)
    pep, psms = parse_mzid('/home/user/csf/calib/peaks.mzid')
    targets, pools = [], {}
    for run, byidx in psms.items():
        mgf = read_mgf(f'/home/user/csf/calib/{run}.mgf')
        pools[run] = [np.sort(m) for _, m in mgf.values()]
        for idx, (pref, z, calcmz) in byidx.items():
            if idx not in mgf:
                continue
            seq, mods = pep[pref]
            mzs = np.sort(mgf[idx][1])
            targets.append((run, z, calcmz * z - z * PROTON, mzs,
                            true_pairs(seq, mods, mzs), len(seq)))
    print(f'identified target spectra: {len(targets)}')
    print(f'decoy masses per spectrum: {NDEC}   null beats-all rate '
          f'{100/(NDEC+1):.2f}%\n')

    tp = np.array([t[4] for t in targets])
    zs = np.array([t[1] for t in targets])
    print(f'{"co-isolated":>12}{"peaks":>8}{"target":>9}{"decoy":>8}{"ratio":>8}'
          f'{"beats-all":>11}{"z=2":>9}{"tp>=2":>9}')
    for k in KS:
        tcnt = np.zeros(len(targets))
        beats = np.zeros(len(targets), int)
        npk = np.zeros(len(targets))
        dmean = np.zeros(len(targets))
        for i, (run, z, M, mzs, _, _) in enumerate(targets):
            if k > 1:
                extra = [mzs] + [pools[run][rng.randrange(len(pools[run]))]
                                 for _ in range(k - 1)]
                mix = np.sort(np.concatenate(extra))
            else:
                mix = mzs
            npk[i] = len(mix)
            s0 = M + 2 * PROTON
            t = count_11(mix, s0)
            d = np.array([count_11(mix, s0 + sh) for sh in SHIFTS])
            tcnt[i] = t
            dmean[i] = d.mean()
            beats[i] = int((d >= t).sum())
        top = beats == 0
        print(f'{k:>12}{npk.mean():>8.0f}{tcnt.mean():>9.2f}{dmean.mean():>8.2f}'
              f'{tcnt.mean()/max(dmean.mean(),1e-9):>8.1f}'
              f'{100*top.mean():>10.1f}%{100*top[zs==2].mean():>8.1f}%'
              f'{100*top[tp>=2].mean():>8.1f}%')
    print(f'\n(z=2 column: n={int((zs==2).sum())}.  '
          f'tp>=2 column: n={int((tp>=2).sum())} spectra that actually contain '
          f'two or more true b/y complementary pairs.)')


if __name__ == '__main__':
    main()
