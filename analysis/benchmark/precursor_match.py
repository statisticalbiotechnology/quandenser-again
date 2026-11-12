#!/usr/bin/env python3
"""
How well an MS1 feature list covers the run's own MS2 precursors.

A precursor counts as matched when a feature has the same charge, an m/z
within TOL_PPM of the precursor's monoisotopic m/z, and an RT range
[rtStart, rtEnd] that contains the precursor's RT. That is the condition under
which Quandenser can tie the MS2 spectrum to an MS1 feature at all, so it is a
fast proxy for what feature detection does to consensus spectra and
identifications.

Precursors not matched at M are also tested at M+1 and M-1 (one 13C spacing
up or down at the same charge). A precursor matched only at M+1 means the
feature detector placed its monoisotopic peak one isotope too high; that is
the error behind biosaur2's shifted-cosine branch above 1500 Da.

Precursor sources:
  *.d/ or analysis.tdf   Bruker DDA-PASEF: Precursors.MonoisotopicMz and
                         Charge, RT from the parent frame (minutes)
  *.mzML                 MS2 spectra: selected ion m/z, charge state, scan
                         start time (converted to minutes if in seconds)

Feature lists are read by column name, so both biosaur2's own output and the
Dinosaur-layout file the pipeline writes work.

    precursor_match.py <precursors: .d | analysis.tdf | .mzML> <features.tsv> [label]
"""
import csv
import os
import re
import sqlite3
import sys

import numpy as np

TOL_PPM = 10.0
C13 = 1.003355


def bruker_precursors(path):
    tdf = path if path.endswith('.tdf') else os.path.join(path, 'analysis.tdf')
    con = sqlite3.connect(tdf)
    rows = con.execute(
        'SELECT p.MonoisotopicMz, p.Charge, f.Time FROM Precursors p '
        'JOIN Frames f ON f.Id = p.Parent '
        'WHERE p.MonoisotopicMz IS NOT NULL AND p.Charge IS NOT NULL').fetchall()
    return np.array([(mz, z, t / 60.0) for mz, z, t in rows], dtype=float)


CV_NAME = re.compile(rb'name="([^"]+)"')
CV_VALUE = re.compile(rb'value="([^"]*)"')
CV_UNIT = re.compile(rb'unitName="([^"]+)"')
WANTED = {b'ms level': 'level', b'scan start time': 'rt',
          b'selected ion m/z': 'mz', b'charge state': 'z'}


def mzml_precursors(path):
    """MS2 precursors from an mzML. Attribute order differs between writers
    (ThermoRawFileParser puts value before name), so name and value are
    matched separately rather than as one pattern."""
    out = []
    cur = None
    with open(path, 'rb') as fh:
        for line in fh:
            if b'<spectrum ' in line:
                cur = {}
            elif cur is None:
                continue
            elif b'</spectrum>' in line:
                if cur.get('level') == 2 and all(k in cur for k in ('mz', 'z', 'rt')):
                    out.append((cur['mz'], cur['z'], cur['rt']))
                cur = None
            elif b'<cvParam' in line:
                m = CV_NAME.search(line)
                key = WANTED.get(m.group(1)) if m else None
                if key is None or key in cur:
                    continue
                v = CV_VALUE.search(line).group(1)
                if key == 'rt':
                    u = CV_UNIT.search(line)
                    cur[key] = float(v) / 60.0 if u and u.group(1) == b'second' else float(v)
                elif key in ('level', 'z'):
                    cur[key] = int(v)
                else:
                    cur[key] = float(v)
    return np.array(out, dtype=float).reshape(-1, 3)


def features(path):
    with open(path) as fh:
        r = csv.DictReader(fh, delimiter='\t')
        f = [(float(x['mz']), int(float(x['charge'])), float(x['rtStart']), float(x['rtEnd']))
             for x in r]
    a = np.array(f, dtype=float).reshape(-1, 4)
    return a[np.argsort(a[:, 0], kind='stable')]


def matched(feat, mz, z, rt):
    """Boolean per precursor: a feature at this m/z, charge and RT exists."""
    fmz = feat[:, 0]
    lo = np.searchsorted(fmz, mz * (1 - TOL_PPM * 1e-6), 'left')
    hi = np.searchsorted(fmz, mz * (1 + TOL_PPM * 1e-6), 'right')
    hit = np.zeros(len(mz), dtype=bool)
    for i in range(len(mz)):
        if hi[i] <= lo[i]:
            continue
        c = feat[lo[i]:hi[i]]
        hit[i] = np.any((c[:, 1] == z[i]) & (c[:, 2] <= rt[i]) & (c[:, 3] >= rt[i]))
    return hit


def main():
    src, ftsv = sys.argv[1], sys.argv[2]
    label = sys.argv[3] if len(sys.argv) > 3 else os.path.basename(ftsv)
    pre = mzml_precursors(src) if src.lower().endswith('.mzml') else bruker_precursors(src)
    feat = features(ftsv)
    mz, z, rt = pre[:, 0], pre[:, 1], pre[:, 2]
    at_m = matched(feat, mz, z, rt)
    at_up = matched(feat, mz + C13 / z, z, rt) & ~at_m
    at_dn = matched(feat, mz - C13 / z, z, rt) & ~at_m & ~at_up
    n = len(mz)
    heavy = mz * z >= 1500
    # why the rest missed: right m/z and charge but outside the RT range, right
    # m/z but another charge, or no feature at that m/z at all
    miss = ~(at_m | at_up | at_dn)
    fmz = feat[:, 0]
    lo = np.searchsorted(fmz, mz * (1 - TOL_PPM * 1e-6), 'left')
    hi = np.searchsorted(fmz, mz * (1 + TOL_PPM * 1e-6), 'right')
    why = {'rt': 0, 'charge': 0, 'none': 0}
    gaps = []
    for i in np.flatnonzero(miss):
        c = feat[lo[i]:hi[i]]
        same_z = c[c[:, 1] == z[i]]
        if len(same_z):
            why['rt'] += 1
            gaps.append(np.min(np.maximum(same_z[:, 2] - rt[i], rt[i] - same_z[:, 3])))
        elif len(c):
            why['charge'] += 1
        else:
            why['none'] += 1
    gap = f', median {60 * np.median(gaps):.0f} s outside' if gaps else ''
    print(f'  misses: RT outside range {100 * why["rt"] / n:.1f}%{gap} '
          f'| other charge {100 * why["charge"] / n:.1f}% | no feature at that m/z {100 * why["none"] / n:.1f}%')
    summary = (f'{label}: {len(feat)} features, {n} precursors | matched at M {100 * at_m.mean():.1f}% '
          f'| only at M+1 {100 * at_up.mean():.1f}% | only at M-1 {100 * at_dn.mean():.1f}% '
          f'| >=1500 Da: M {100 * at_m[heavy].mean():.1f}%, M+1 {100 * at_up[heavy].mean():.1f}%')
    print(summary)


if __name__ == '__main__':
    main()
