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

    biosaur2_to_dinosaur.py biosaur2.features.tsv > <stem>.features.tsv
"""
import sys

DINOSAUR = ['mz', 'mostAbundantMz', 'charge', 'rtStart', 'rtApex', 'rtEnd',
            'fwhm', 'nIsotopes', 'nScans', 'averagineCorr', 'mass', 'massCalib',
            'intensityApex', 'intensitySum']
FROM = {'mostAbundantMz': 'mz', 'mass': 'massCalib'}
PLACEHOLDER = {'fwhm': '0', 'averagineCorr': '0'}


def main(path):
    with open(path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        ix = {name: i for i, name in enumerate(header)}
        missing = [c for c in DINOSAUR if c not in PLACEHOLDER
                   and FROM.get(c, c) not in ix]
        if missing:
            sys.exit(f'{path}: biosaur2 output lacks {missing}')
        im = ix.get('im')
        out = sys.stdout
        out.write('\t'.join(DINOSAUR + (['ionMobility'] if im is not None else [])) + '\n')
        n = 0
        for line in fh:
            f = line.rstrip('\n').split('\t')
            cols = [PLACEHOLDER[c] if c in PLACEHOLDER else f[ix[FROM.get(c, c)]]
                    for c in DINOSAUR]
            if im is not None:
                cols.append(f[im])
            out.write('\t'.join(cols) + '\n')
            n += 1
    print(f'{n} features rewritten in Dinosaur column order', file=sys.stderr)


if __name__ == '__main__':
    main(sys.argv[1])
