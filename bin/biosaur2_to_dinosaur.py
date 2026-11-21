#!/usr/bin/env python3
"""
Rewrite a biosaur2 features.tsv in the column layout Quandenser reads.

Quandenser parses feature files by column *position*, in Dinosaur's order
(src/DinosaurIO.cpp, parseDinosaurFeatureRow), and uses six fields: mz,
charge, rtStart, rtApex, rtEnd and intensitySum. biosaur2 writes the same
names in a different order, so its file would be read without error and with
every value in the wrong place. The remaining Dinosaur columns are filled with
placeholders, since nothing reads them.

biosaur2's ion mobility (`im`, 1/K0) is appended as a 15th column,
ionMobility, which Quandenser reads when matching features between runs with
--align-im-tol. Quandenser builds without that option stop reading after
intensitySum and ignore it. The other extras (FAIMS, hill lists, isotope
errors) are dropped.

With MIN_SCANS, features whose monoisotopic hill spans fewer MS1 scans are
dropped. On two-hour timsTOF runs biosaur2 reports about 1.8 million
features per run, and Quandenser's feature matching then outgrows 32 GB;
MIN_SCANS 6 keeps 54% of them and loses under one point of MS2 precursor
coverage (analysis/RESULTS-biosaur2.md).

    biosaur2_to_dinosaur.py biosaur2.features.tsv [MIN_SCANS] > <stem>.features.tsv
"""
import sys

DINOSAUR = ['mz', 'mostAbundantMz', 'charge', 'rtStart', 'rtApex', 'rtEnd',
            'fwhm', 'nIsotopes', 'nScans', 'averagineCorr', 'mass', 'massCalib',
            'intensityApex', 'intensitySum']
FROM = {'mostAbundantMz': 'mz', 'mass': 'massCalib'}
PLACEHOLDER = {'fwhm': '0', 'averagineCorr': '0'}


def main(path, min_scans=0):
    with open(path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        ix = {name: i for i, name in enumerate(header)}
        missing = [c for c in DINOSAUR if c not in PLACEHOLDER
                   and FROM.get(c, c) not in ix]
        if missing:
            sys.exit(f'{path}: biosaur2 output lacks {missing}')
        im = ix.get('im')
        scans = ix['nScans']
        out = sys.stdout
        out.write('\t'.join(DINOSAUR + (['ionMobility'] if im is not None else [])) + '\n')
        n = dropped = 0
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if int(f[scans]) < min_scans:
                dropped += 1
                continue
            cols = [PLACEHOLDER[c] if c in PLACEHOLDER else f[ix[FROM.get(c, c)]]
                    for c in DINOSAUR]
            if im is not None:
                cols.append(f[im])
            out.write('\t'.join(cols) + '\n')
            n += 1
    kept = f', {dropped} with fewer than {min_scans} scans dropped' if min_scans else ''
    print(f'{n} features rewritten in Dinosaur column order{kept}', file=sys.stderr)


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 0)
