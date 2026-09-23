#!/usr/bin/env python3
"""
DIA-NN's own identifications, rewritten as the two files diatest.py reads.

The channel-2 test scores library peptides against decoys, and what counts as
a library peptide is decided by the OpenSWATH spectral library and its
mProphet peak groups. DIA-NN, run library-free on the same file, identifies
considerably more peptides than that library contains, so its output is a
second truth set worth reporting next to the first.

It is not a better one, and should not replace it. DIA-NN identifies by
fragment co-elution and retention-time prediction, which is channel 1, and the
fragments in its generated library were selected and calibrated on this very
run. Scoring channel 2 against them is therefore optimistic in a way the
OpenSWATH set, which comes from separate DDA runs, is not quite. Report both
and let the difference between them say how much the truth set is doing.

    python3 diann_as_library.py <diann_out_dir> <run> <out_dir>

Writes <out_dir>/Ecoli_DDA_CombinedLib.sptxt and
<out_dir>/split_<run>_mprophet_all_peakgroups_DDA.xls, the names diatest.py
expects, so no analysis script has to change.
"""
import collections
import csv
import os
import sys

DIANN, RUN, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
QCUT = 0.01


def library(path):
    """(sequence, charge) -> precursor m/z, RT, [(m/z, intensity, annotation)]"""
    entries = collections.defaultdict(list)
    meta = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            if r.get('decoy', '0') not in ('0', 'False', 'FALSE', ''):
                continue
            if r.get('FragmentLossType', 'noloss') not in ('noloss', ''):
                continue
            kind = r['FragmentType']
            if kind not in ('b', 'y'):
                continue
            key = (r['PeptideSequence'], int(r['PrecursorCharge']))
            z = int(r['FragmentCharge'])
            ann = f"{kind}{int(r['FragmentSeriesNumber'])}" + (f'^{z}' if z > 1 else '')
            entries[key].append((float(r['ProductMz']),
                                 float(r['LibraryIntensity']), ann))
            meta[key] = (float(r['PrecursorMz']), float(r['Tr_recalibrated']))
    return entries, meta


def write_sptxt(entries, meta, path):
    with open(path, 'w') as fh:
        for (seq, z), peaks in sorted(entries.items()):
            mz, rt = meta[(seq, z)]
            peaks = sorted(peaks)
            fh.write(f'Name: {seq}/{z}\n')
            fh.write(f'LibID: {seq}_{z}\n')
            fh.write(f'PrecursorMZ: {mz:.6f}\n')
            fh.write(f'Comment: RetentionTime={rt:.4f} Spec=DIA-NN\n')
            fh.write(f'NumPeaks: {len(peaks)}\n')
            for pmz, inten, ann in peaks:
                fh.write(f'{pmz:.6f}\t{inten:.4f}\t{ann}\n')
            fh.write('\n')
    return len(entries)


def write_peakgroups(report, entries, meta, path):
    """report.tsv -> the mProphet columns diatest.py reads.

    RT, RT.Start and RT.Stop are in minutes; leftWidth and rightWidth are read
    as seconds, which is the unit trap the analysis README warns about.
    """
    cols = ['transition_group_id', 'decoy', 'Sequence', 'Charge', 'm.z', 'RT',
            'leftWidth', 'rightWidth', 'm_score', 'var_xcorr_shape']
    best = {}
    with open(report) as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            q = float(r['Q.Value'])
            if q > QCUT:
                continue
            key = (r['Stripped.Sequence'], int(r['Precursor.Charge']))
            if key not in entries or key not in meta:
                continue
            if key in best and best[key][0] <= q:
                continue
            best[key] = (q, r)
    with open(path, 'w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for (seq, z), (q, r) in sorted(best.items()):
            fh.write('\t'.join([
                # report.tsv carries no precursor m/z; the library does
                f'{seq}_{z}', 'FALSE', seq, str(z), f'{meta[(seq, z)][0]:.6f}',
                f"{float(r['RT']) * 60.0:.3f}",
                f"{float(r['RT.Start']) * 60.0:.3f}",
                f"{float(r['RT.Stop']) * 60.0:.3f}",
                f'{q:.6g}',
                '0',            # var_xcorr_shape: OpenSWATH-only, unused here
            ]) + '\n')
    return len(best)


def main():
    os.makedirs(OUT, exist_ok=True)
    entries, meta = library(os.path.join(DIANN, 'lib.tsv'))
    n = write_sptxt(entries, meta, os.path.join(OUT, 'Ecoli_DDA_CombinedLib.sptxt'))
    m = write_peakgroups(os.path.join(DIANN, 'report.tsv'), entries, meta,
                         os.path.join(OUT, f'split_{RUN}_mprophet_all_peakgroups_DDA.xls'))
    print(f'{n} precursors written to the library, '
          f'{m} peak groups at q <= {QCUT}')
    print('note: var_xcorr_shape is written as 0; analyse.py uses it, '
          'ch2_sweep.py does not')


if __name__ == '__main__':
    main()
