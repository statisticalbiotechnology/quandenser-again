"""Write the MS2 spectra of one tdf2mzml mzML as .ms2 chunks with the
instrument's precursor m/z and charge, plus a header table.

    pasef_extract_ms2.py <mzML> <out_dir> [spectra_per_chunk]

Header table columns: scan rt_min im precmz charge iso_lo iso_hi npeaks.
Scan is N from the native id "index=N". Spectra without an instrument charge
get Z lines for charges 2 and 3.
"""
import os
import sys

from pyteomics import mzml

PROTON = 1.00727646688

src, out = sys.argv[1], sys.argv[2]
per_chunk = int(sys.argv[3]) if len(sys.argv) > 3 else 60000
os.makedirs(out, exist_ok=True)
hdr = open(os.path.join(out, 'ms2_header.tsv'), 'w')
hdr.write('scan\trt_min\tim\tprecmz\tcharge\tiso_lo\tiso_hi\tnpeaks\n')
chunk, n, fh = 0, 0, None
for s in mzml.read(src, decode_binary=True):
    if s.get('ms level') != 2:
        continue
    scan = int(s['id'].split('=')[1])
    sc = s['scanList']['scan'][0]
    rt = float(sc['scan start time'])
    im = float(sc.get('inverse reduced ion mobility', 'nan'))
    p = s['precursorList']['precursor'][0]
    ion = p['selectedIonList']['selectedIon'][0]
    mz = float(ion['selected ion m/z'])
    z = int(ion['charge state']) if 'charge state' in ion else 0
    iw = p['isolationWindow']
    lo = float(iw.get('isolation window lower offset', 0))
    hi = float(iw.get('isolation window upper offset', 0))
    mzs, ints = s['m/z array'], s['intensity array']
    hdr.write(f'{scan}\t{rt:.5f}\t{im:.4f}\t{mz:.5f}\t{z}\t{lo}\t{hi}\t{len(mzs)}\n')
    if n % per_chunk == 0:
        if fh:
            fh.close()
        chunk += 1
        fh = open(os.path.join(out, f'raw.part{chunk}.ms2'), 'w')
    n += 1
    fh.write(f'S\t{scan}\t{scan}\t{mz:.5f}\n')
    fh.write(f'I\tRTime\t{rt:.5f}\n')
    for c in ([z] if z > 0 else [2, 3]):
        fh.write(f'Z\t{c}\t{(mz - PROTON) * c + PROTON:.5f}\n')
    fh.write(''.join(f'{a:.5f} {b:.1f}\n' for a, b in zip(mzs, ints)))
fh.close()
hdr.close()
print(f'{n} MS2 spectra in {chunk} chunks')
