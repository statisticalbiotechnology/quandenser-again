#!/usr/bin/env python3
"""
Score a MaRaCluster clustering with identifications of the raw MS2 spectra.

Spectra are identified by a separate search of the unclustered runs (Crux Tide +
Percolator, q <= Q). Over pairs of identified spectra:

  sensitivity  pairs with the same peptide and charge that share a cluster,
               out of all such pairs
  impurity     pairs that share a cluster but carry different peptides, out of
               all identified pairs that share a cluster

A spectrum can sit in several clusters (Quandenser splits mass-charge states);
a pair counts as clustered when it shares any of them. I and L are treated as
one residue.

    cluster_ids.py <tide.log> <percolator.target.psms.txt> <clusters.tsv> [<clusters.tsv> ...]

tide.log supplies the file order behind Percolator's file_idx.
"""
import collections
import os
import re
import sys

Q = 0.01


def file_order(tide_log):
    seen = []
    for line in open(tide_log):
        m = re.search(r'Converting (\S+) to spectrumrecords', line) or re.search(r'Reading (\S+\.mzML)', line)
        if m and os.path.basename(m.group(1)) not in seen:
            seen.append(os.path.basename(m.group(1)))
    return seen


def identifications(psms, files):
    ids = {}
    with open(psms) as fh:
        head = fh.readline().rstrip('\n').split('\t')
        c = {k: head.index(k) for k in ('file_idx', 'scan', 'charge', 'percolator q-value', 'sequence')}
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if float(p[c['percolator q-value']]) > Q:
                continue
            seq = p[c['sequence']].replace('I', 'L')
            ids[(files[int(p[c['file_idx']])], int(p[c['scan']]))] = (seq, int(p[c['charge']]))
    return ids


def score(clusters, ids):
    member, sizes = collections.defaultdict(list), collections.Counter()
    spectra = set()
    cid = 0
    with open(clusters) as fh:
        for line in fh:
            p = line.split('\t')
            if len(p) < 2:
                cid += 1
                continue
            key = (os.path.basename(p[0]), int(p[1]))
            spectra.add(key)
            sizes[cid] += 1
            if key in ids:
                member[cid].append(key)
    index = {k: i for i, k in enumerate(ids)}
    same, diff = set(), set()
    for keys in member.values():
        keys = sorted(set(keys), key=index.get)
        for a in range(len(keys)):
            for b in range(a + 1, len(keys)):
                pair = (index[keys[a]], index[keys[b]])
                sa, sb = ids[keys[a]], ids[keys[b]]
                if sa == sb:
                    same.add(pair)
                elif sa[0] != sb[0]:
                    diff.add(pair)
    by_pep = collections.Counter(ids[k] for k in ids if k in spectra)
    total = sum(n * (n - 1) // 2 for n in by_pep.values())
    singletons = sum(1 for s in sizes.values() if s == 1)
    print(f'{os.path.basename(clusters)}: {len(spectra)} spectra, {len(sizes)} clusters, '
          f'{100 * singletons / len(spectra):.1f}% singletons | identified {sum(by_pep.values())} '
          f'| sensitivity {100 * len(same) / max(total, 1):.1f}% of {total} same-peptide pairs '
          f'| impurity {100 * len(diff) / max(len(same) + len(diff), 1):.2f}%')


def main():
    files = file_order(sys.argv[1])
    ids = identifications(sys.argv[2], files)
    for c in sys.argv[3:]:
        score(c, ids)


if __name__ == '__main__':
    main()
