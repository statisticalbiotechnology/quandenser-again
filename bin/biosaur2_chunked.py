#!/usr/bin/env python3
"""
Run biosaur2 on a long run in overlapping retention-time chunks.

biosaur2 holds a whole run's hills and candidate isotope clusters in memory:
about 10 GB per GB of raw ion-mobility MS1 mzML on PXD039875, so a two-hour
timsTOF run would need some 70 GB. It reads only the frames inside
-rtmin/-rtmax, though, so a chunk needs memory in proportion to its length.

Chunk i covers [start - OVERLAP, end + OVERLAP] and keeps the features whose
apex lies in its core [start, end). Every feature is therefore reported once,
by a chunk that saw all of it, as long as it reaches no further than OVERLAP
from its apex.

    biosaur2_chunked.py <ms1.mzML> <out.tsv> <chunk minutes> [biosaur2 args ...]
"""
import os
import re
import subprocess
import sys
import tempfile

OVERLAP = 2.0  # minutes
# attribute order differs between writers, so the cvParam is found first
RT = re.compile(rb'<cvParam[^>]*name="scan start time"[^>]*>')
VALUE = re.compile(rb'value="([^"]+)"')
UNIT = re.compile(rb'unitName="(\w+)"')


def last_rt(path):
    """Retention time of the last spectrum, in minutes, read from the file's end."""
    size = os.path.getsize(path)
    block = 1 << 24
    with open(path, 'rb') as fh:
        while True:
            fh.seek(max(0, size - block))
            hits = RT.findall(fh.read())
            if hits or block >= size:
                break
            block *= 4
    if not hits:
        sys.exit(f'{path}: no scan start time found')
    value = float(VALUE.search(hits[-1]).group(1))
    unit = UNIT.search(hits[-1])
    return value / 60.0 if unit and unit.group(1) == b'second' else value


def main():
    mzml, out, chunk = sys.argv[1], sys.argv[2], float(sys.argv[3])
    args = sys.argv[4:]
    end_rt = last_rt(mzml)
    header, n = None, 0
    with open(out, 'w') as dst, tempfile.TemporaryDirectory(dir='.') as tmp:
        start = 0.0
        while start <= end_rt:
            end = start + chunk
            part = os.path.join(tmp, 'chunk.tsv')
            cmd = ['biosaur2', mzml, '-o', part, '-rtmin', str(max(0.0, start - OVERLAP))]
            if end < end_rt:
                cmd += ['-rtmax', str(end + OVERLAP)]
            print(f'chunk {start:.1f}-{min(end, end_rt):.1f} min: {" ".join(cmd + args)}',
                  file=sys.stderr, flush=True)
            subprocess.run(cmd + args, check=True)
            if os.path.exists(part):
                with open(part) as src:
                    first = src.readline()
                    if first:
                        if header is None:
                            header = first
                            dst.write(header)
                            apex = header.rstrip('\n').split('\t').index('rtApex')
                        for line in src:
                            rt = float(line.split('\t', apex + 1)[apex])
                            if start <= rt < end or (end >= end_rt and rt >= start):
                                dst.write(line)
                                n += 1
                os.remove(part)
            start = end
    print(f'{n} features from {mzml} in chunks of {chunk:g} min', file=sys.stderr)


if __name__ == '__main__':
    main()
