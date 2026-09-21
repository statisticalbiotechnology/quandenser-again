#!/usr/bin/env python3
"""
Both channels measured on real DIA data, at the fragment level, with ground
truth from the study's own DDA spectral library and OpenSWATH peak groups.

Channel 1 (co-elution).  For each confidently identified peptide, correlate
each candidate fragment's MS2 trace with the peptide's MS1 precursor trace.
True fragments are the library's annotated b/y ions for that peptide; decoy
fragments are the library ions of a different peptide whose precursor falls in
the same isolation window but which elutes far away.  Using the MS1 trace as
the reference keeps this non-circular.

Channel 2 (complementary pairs).  On the summed MS2 spectrum of the same
isolation window over the peptide's elution, count peak pairs summing to
M + 2*proton, against 110 integer-Da shifted masses.  OpenSWATH's own decoys
are shuffled sequences, which preserve precursor mass, so they are useless as
a null here; the mass-shift null is used instead.

Independence.  The two scores are computed on the same peptides, so their
rank correlation answers directly whether combining them can help.
"""
import sys, re, random, pickle, collections, zlib
import numpy as np
from pyteomics import mzxml

RUN = sys.argv[1] if len(sys.argv) > 1 else '18484_REP3_1ug_Ecoli_NewStock2_SWATH_1'
MZXML = RUN + '.mzXML'
MPROPH = f'split_{RUN}_mprophet_all_peakgroups_DDA.xls'
OUTPKL = f'traces_{RUN}.pkl'
SPTXT = 'Ecoli_DDA_CombinedLib.sptxt'
PROTON = 1.00727646
QCUT = 0.01
FRAG_TOL = 0.05          # TripleTOF MS2, absolute Da
MAX_PEP = 6000
SHIFTS = np.array([s * d for d in range(6, 61) for s in (-1, 1)], float)


def load_library():
    lib = {}
    name = mz = librt = None
    peaks = []
    for line in open(SPTXT):
        if line.startswith('Name:'):
            if name and peaks:
                lib[name] = (mz, peaks, librt)
            name = line.split(':', 1)[1].strip()
            peaks = []
            librt = None
        elif line.startswith('PrecursorMZ:'):
            mz = float(line.split(':', 1)[1])
        elif line.startswith('Comment:') and 'RetentionTime=' in line:
            try:
                librt = float(line.split('RetentionTime=')[1].split()[0].split(',')[0])
            except (ValueError, IndexError):
                librt = None
        elif line and line[0].isdigit() and name:
            f = line.split('\t')
            if len(f) >= 3:
                ann = f[2].split(',')[0]
                if re.match(r'^[by]\d+(\^\d+)?(/|$)', ann):
                    peaks.append((float(f[0]), float(f[1]), ann.split('/')[0]))
    if name and peaks:
        lib[name] = (mz, peaks, librt)
    return lib


def load_peptides(lib):
    rows = open(MPROPH).read().rstrip('\n').split('\n')
    hdr = rows[0].split('\t'); ix = {k: i for i, k in enumerate(hdr)}
    out = {}
    for r in rows[1:]:
        f = r.split('\t')
        if f[ix['decoy']].strip().upper() != 'FALSE':
            continue
        try:
            q = float(f[ix['m_score']])
            if q > QCUT:
                continue
            key = f[ix['Sequence']] + '/' + f[ix['Charge']]
            if key not in lib:
                continue
            rec = dict(seq=f[ix['Sequence']], z=int(f[ix['Charge']]),
                       mz=float(f[ix['m.z']]), rt=float(f[ix['RT']]),
                       lo=float(f[ix['leftWidth']]), hi=float(f[ix['rightWidth']]),
                       q=q, xcorr=float(f[ix['var_xcorr_shape']]), key=key)
        except (ValueError, KeyError, IndexError):
            continue
        if key not in out or q < out[key]['q']:
            out[key] = rec
    return list(out.values())


def main():
    lib = load_library()
    print(f'library entries with annotated b/y ions: {len(lib)}')
    peps = load_peptides(lib)
    print(f'confident peak groups (m_score < {QCUT}) matched to the library: {len(peps)}')
    rng = random.Random(0)
    peps.sort(key=lambda p: p['key'])
    rng.shuffle(peps)
    peps = peps[:MAX_PEP]

    # window scheme from the first cycle
    windows, seen = [], set()
    n1 = 0
    for s in mzxml.read(MZXML):
        if s['msLevel'] == 1:
            n1 += 1
            if n1 > 1:
                break
            continue
        p = s['precursorMz'][0]
        c = float(p['precursorMz']); w = float(p.get('windowWideness', 25.0))
        if (c, w) not in seen:
            seen.add((c, w)); windows.append((c - w / 2, c + w / 2))
    windows.sort()
    print(f'{len(windows)} isolation windows, '
          f'{windows[0][0]:.1f}-{windows[-1][1]:.1f} m/z, '
          f'width {windows[0][1]-windows[0][0]:.1f}\n')

    def win_of(mz):
        for i, (a, b) in enumerate(windows):
            if a <= mz < b:
                return i
        return -1

    for p in peps:
        p['w'] = win_of(p['mz'])
    peps = [p for p in peps if p['w'] >= 0]
    bywin = collections.defaultdict(list)
    for i, p in enumerate(peps):
        bywin[p['w']].append(i)

    # Decoy fragments: library ions of a different peptide whose precursor falls
    # in the same isolation window but which elutes far away.  The partner is
    # picked from the LIBRARY by a deterministic rule, not from this run's
    # confident list, so a peptide gets the same decoy fragments in every run.
    # Drawing it per run instead makes the decoy columns incomparable and
    # silently destroys any cross-run pooling comparison.
    libwin = collections.defaultdict(list)
    for key, (lmz, lpk, lrt) in lib.items():
        if lrt is None:
            continue
        w = win_of(lmz)
        if w >= 0:
            libwin[w].append((key, lrt))
    for w in libwin:
        libwin[w].sort()

    for p in peps:
        p['true_frag'] = np.array(sorted(m for m, _, _ in lib[p['key']][1]))
        pool = [(k, r) for k, r in libwin.get(p['w'], [])
                if k != p['key'] and abs(r - lib[p['key']][2]) > 600]
        if pool:
            # zlib.crc32, not hash(): str hashing is randomised per
            # process, so the two runs would get different decoys.
            k = pool[zlib.crc32(p['key'].encode()) % len(pool)][0]
            p['decoy_src'] = k
            p['decoy_frag'] = np.array(sorted(m for m, _, _ in lib[k][1]))
        else:
            p['decoy_src'] = None
            p['decoy_frag'] = np.array([])
        p['ms1'] = []
        p['tr_true'] = []
        p['tr_decoy'] = []
        p['sum_mz'] = []
        p['sum_int'] = []
        p['rts'] = []
    print(f'{len(peps)} peptides assigned to a window, '
          f'{sum(1 for p in peps if p["decoy_src"] is not None)} with a decoy partner')

    pad = 45.0
    lo = np.array([p['lo'] - pad for p in peps])
    hi = np.array([p['hi'] + pad for p in peps])
    wof = np.array([p['w'] for p in peps])
    pmz = np.array([p['mz'] for p in peps])

    def pick(smz, sint, targets, tol):
        if smz.size == 0 or len(targets) == 0:
            return np.zeros(len(targets))
        j = np.searchsorted(smz, targets)
        a = np.clip(j - 1, 0, smz.size - 1); b = np.clip(j, 0, smz.size - 1)
        da = np.abs(smz[a] - targets); db = np.abs(smz[b] - targets)
        k = np.where(da < db, a, b); d = np.minimum(da, db)
        return np.where(d <= tol, sint[k], 0.0)

    nms1 = nms2 = 0
    for s in mzxml.read(MZXML):
        # pyteomics reports mzXML retentionTime in MINUTES; OpenSWATH's RT,
        # leftWidth and rightWidth are in seconds.  Mixing them matches nothing
        # at all, silently.
        rt = float(s['retentionTime']) * 60.0
        smz = np.asarray(s['m/z array'], float)
        sint = np.asarray(s['intensity array'], float)
        if s['msLevel'] == 1:
            nms1 += 1
            act = np.where((rt >= lo) & (rt <= hi))[0]
            if act.size:
                v = pick(smz, sint, pmz[act], 0.05)
                for k, i in enumerate(act):
                    peps[i]['ms1'].append(v[k])
                    peps[i]['rts'].append(rt)
        else:
            nms2 += 1
            c = float(s['precursorMz'][0]['precursorMz'])
            w = win_of(c)
            if w < 0:
                continue
            act = [i for i in bywin[w] if lo[i] <= rt <= hi[i]]
            for i in act:
                p = peps[i]
                p['tr_true'].append(pick(smz, sint, p['true_frag'], FRAG_TOL))
                if p['decoy_frag'].size:
                    p['tr_decoy'].append(pick(smz, sint, p['decoy_frag'], FRAG_TOL))
                if p['lo'] <= rt <= p['hi']:
                    p['sum_mz'].append(smz); p['sum_int'].append(sint)
    print(f'streamed {nms1} MS1 and {nms2} MS2 scans')
    for p in peps:
        for k in ('ms1', 'rts'):
            p[k] = np.array(p[k])
        p['tr_true'] = np.array(p['tr_true']) if p['tr_true'] else np.zeros((0, 0))
        p['tr_decoy'] = np.array(p['tr_decoy']) if p['tr_decoy'] else np.zeros((0, 0))
        # Merge the MS2 scans across the elution peak into one spectrum rather
        # than concatenating them: concatenated scans repeat every fragment once
        # per cycle, which multiplies the accidental pair count without adding
        # any information.  Peaks within FRAG_TOL are summed.
        if p['sum_mz']:
            mz = np.concatenate(p['sum_mz']); it = np.concatenate(p['sum_int'])
            o = np.argsort(mz); mz, it = mz[o], it[o]
            grp = np.concatenate([[0], np.cumsum(np.diff(mz) > FRAG_TOL)])
            n = grp[-1] + 1
            wsum = np.bincount(grp, weights=mz * it, minlength=n)
            isum = np.bincount(grp, weights=it, minlength=n)
            keep = isum > 0
            p['spec_mz'] = wsum[keep] / isum[keep]
            p['spec_int'] = isum[keep]
        else:
            p['spec_mz'] = np.array([]); p['spec_int'] = np.array([])
        p.pop('sum_mz', None); p.pop('sum_int', None)
    pickle.dump(peps, open(OUTPKL, 'wb'), 4)
    print(f'saved {OUTPKL} ({len(peps)} peptides)')


if __name__ == '__main__':
    main()
