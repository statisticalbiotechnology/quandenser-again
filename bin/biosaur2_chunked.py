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

Chunks are independent, so -j runs that many at a time, each a biosaur2
process of its own; memory grows with the number running. The output does
not depend on -j: chunks are merged in retention-time order, and each
chunk's log is printed when it finishes.

    biosaur2_chunked.py [-j WORKERS] <ms1.mzML> <out.tsv> <chunk minutes> [biosaur2 args ...]
"""
import os
import re
import subprocess
import sys
import tempfile
from concurrent.futures import FIRST_EXCEPTION, ThreadPoolExecutor, wait

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


def run_chunk(mzml, part, start, end, end_rt, args):
    """biosaur2 on [start - OVERLAP, end + OVERLAP]; its log is printed when done."""
    cmd = ['biosaur2', mzml, '-o', part, '-rtmin', str(max(0.0, start - OVERLAP))]
    if end < end_rt:
        cmd += ['-rtmax', str(end + OVERLAP)]
    with open(part + '.log', 'w+') as log:
        done = subprocess.run(cmd + args, stdout=log, stderr=subprocess.STDOUT)
        log.seek(0)
        sys.stderr.write(f'chunk {start:.1f}-{min(end, end_rt):.1f} min: {" ".join(cmd + args)}\n'
                         + log.read())
        sys.stderr.flush()
    if done.returncode:
        raise subprocess.CalledProcessError(done.returncode, cmd + args)


def main():
    argv = sys.argv[1:]
    workers = 1
    if argv[0] == '-j':
        workers, argv = int(argv[1]), argv[2:]
    mzml, out, chunk = argv[0], argv[1], float(argv[2])
    args = argv[3:]
    end_rt = last_rt(mzml)
    cores, start = [], 0.0
    while start <= end_rt:
        cores.append((start, start + chunk))
        start += chunk
    header, n = None, 0
    with open(out, 'w') as dst, tempfile.TemporaryDirectory(dir='.') as tmp:
        parts = [os.path.join(tmp, f'chunk{i}.tsv') for i in range(len(cores))]
        with ThreadPoolExecutor(max_workers=workers) as pool:
            jobs = [pool.submit(run_chunk, mzml, part, s, e, end_rt, args)
                    for part, (s, e) in zip(parts, cores)]
            finished, _ = wait(jobs, return_when=FIRST_EXCEPTION)
            failed = [j for j in finished if j.exception()]
            if failed:
                pool.shutdown(cancel_futures=True)
                raise failed[0].exception()
        for part, (start, end) in zip(parts, cores):
            if not os.path.exists(part):
                continue
            with open(part) as src:
                first = src.readline()
                if not first:
                    continue
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
    print(f'{n} features from {mzml} in chunks of {chunk:g} min', file=sys.stderr)


if __name__ == '__main__':
    main()
