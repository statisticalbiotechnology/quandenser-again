# The ppm consensus merge, measured

Question: does MaRaCluster's ppm consensus merge
(`--consensus-method ppm`, statisticalbiotechnology/maracluster#35) change
either of the two measurements in `docs/dia-feasibility.md` that depend on how
peaks are merged?

The scripts were written before the merge existed and had never been run on
real data; what follows is the first run of them, on a machine with network
access to PRIDE.

Nothing here edits `docs/dia-feasibility.md`. One bug was found and fixed in
`dia/fragscatter.py`; it is described under Step 3.

Machine: 4 cores (i5-7500), 39 GB. Everything ran in containers; the analysis
image is Ubuntu 24.04 with numpy, scipy, pyteomics, psims, lxml and
scikit-learn, plus g++ and libboost-dev for the port check. The scripts expect
`/home/user/csf` and `/home/user/dia`, which were provided as bind mounts
rather than by editing the scripts.

---

## Step 0: are the Python ports of the two merges faithful

    cd analysis && python3 consensusmerge_check.py 50

    unit test PASSED
    50 trials
      ppm      largest disagreement 2.73e-09 ppm in m/z, 4.64e-16 relative in intensity
      thomson  largest disagreement 0 ppm in m/z, 0 relative in intensity
    identical

It compiles `ext/maracluster/src/{PpmConsensusMerge,MSClusterMerge}.cpp`
unmodified and compares against them, so the numbers below rest on the same
code the pipeline runs.

---

## Step 3: DIA, channel 2 on merged window spectra

Data: PXD001587, `18484_REP3_1ug_Ecoli_NewStock2_SWATH_1`, 1.85 GB mzXML,
1,995 MS1 and 67,830 MS2 scans, 34 fixed isolation windows of 25 Th. Library
and mProphet table taken out of `Ecoli.zip` by HTTP range request, as the
README describes.

### What ppmSigma has to be on this instrument

    cd analysis/dia && python3 fragscatter.py 18484_REP3_1ug_Ecoli_NewStock2_SWATH_1

    3000 scan pairs, 69317 matched peaks within 100 ppm
    median -0.17 ppm, robust sigma 13.39 ppm (1.4826 MAD), q05 -37.20, q95 +35.51

           m/z        n     median      sigma
     112-267      13863      -0.10      12.17
     267-346      13864      -0.18      10.86
     346-386      13863      -0.44      12.82
     386-448      13863      +0.00      12.35
     448-1512     13863      -0.06      21.81

**13.39 ppm, against the 2 ppm default, which is an Orbitrap number.** At the
merge's 4-sigma cut that is a +/-54 ppm window, 0.027 Th at m/z 500, the same
order as the fixed 0.05 Th grouping the DIA scripts use anyway. On this
instrument a ppm-scale merge is therefore not obviously the right shape: the
tolerance it needs is wide enough that "constant in ppm" and "constant in Th"
are nearly the same window over the range where most fragments are. The scatter
is flat in ppm to about m/z 450 and then roughly doubles, which is the one
feature a fixed-width grouping cannot follow.

Independent corroboration: DIA-NN, run on the same file (below), chose
12.42 ppm as its recommended MS1 mass accuracy without being told anything
about this.

`fragscatter.py` crashed on the first real file it saw:

    ValueError: kth(=-4) out of bounds (36)

Its size guard tests the current scan (`if mz.size < TOPN`), but the peaks come
from the *previous* scan of that window, which the same guard may have stored
with fewer than `TOPN` peaks. Fixed by taking `min(TOPN, pit.size)`, which
preserves the intent. This is the only change made to the analysis scripts.

### Channel 2 under the two merges

    python3 diatest.py 18484_REP3_1ug_Ecoli_NewStock2_SWATH_1
    python3 diatest.py 18484_REP3_1ug_Ecoli_NewStock2_SWATH_1 ppm 13.39
    python3 ch2_sweep.py traces_18484_REP3_1ug_Ecoli_NewStock2_SWATH_1.pkl
    python3 ch2_sweep.py traces_18484_REP3_1ug_Ecoli_NewStock2_SWATH_1_ppm13.39.pkl

4,710 library peptides assigned to a window, 4,707 with a decoy partner. Null
beats-all rate 0.90%.

Peaks per merged window spectrum: centroid median 3,139 (q10 2,540, q90 4,105);
ppm median 14,444 (q10 9,617, q90 19,248).

| top-N | tol | centroid target/decoy | beats-all | enrich | ppm target/decoy | beats-all | enrich |
|---|---|---|---|---|---|---|---|
| 50 | 0.02 | 0.50 / 0.19 | 4.0% | 4.5x | 0.77 / 0.27 | **6.1%** | 6.7x |
| 50 | 0.05 | 0.89 / 0.46 | 4.1% | 4.6x | 1.21 / 0.60 | 4.6% | 5.1x |
| 100 | 0.02 | 1.33 / 0.73 | 3.8% | 4.2x | 2.14 / 1.07 | **6.8%** | 7.6x |
| 100 | 0.05 | 2.62 / 1.72 | 3.7% | 4.1x | 3.81 / 2.44 | 4.0% | 4.5x |
| 200 | 0.02 | 3.78 / 2.81 | 2.5% | 2.7x | 6.24 / 4.14 | **6.4%** | 7.1x |
| 200 | 0.05 | 8.33 / 6.66 | 2.9% | 3.3x | 12.49 / 9.46 | 4.6% | 5.2x |
| 500 | 0.02 | 17.33 / 16.07 | 0.8% | 0.8x | 28.29 / 23.62 | 3.5% | 3.9x |
| 500 | 0.05 | 42.41 / 39.85 | 0.8% | 0.8x | 62.29 / 54.69 | 3.3% | 3.7x |
| all | 0.02 | 172.25 / 172.33 | 0.3% | 0.3x | 5080.56 / 5049.14 | 0.5% | 0.6x |
| all | 0.05 | 437.53 / 436.74 | 0.2% | 0.2x | 12664.87 / 12585.59 | 0.6% | 0.7x |

The ppm merge roughly doubles the beats-all rate at matched denoising, 3.8% to
6.8% at top-100 and 0.02 Da, and it keeps signal much further out: at top-500
the centroid merge has fallen below the null while ppm still shows 3.5%. That
follows from what the merge does — it does not collapse peaks that the
instrument resolved, so the strongest 100 peaks of a merged window spectrum are
more often real fragments rather than centroids of several.

It does not change the conclusion. 6.8% of peptides beating all 110 decoys,
against a 0.90% null, is a real but weak signal, and it is measured on library
peptides with a known answer. The merge makes channel 2 less bad; it does not
make it a grouping mechanism.

### Benchmark: DIA-NN on the same file

    # mzXML -> mzML, isolation windows repaired (see below), then
    diann --f swath1_diann.mzML --fasta ecoli_sp.fasta --fasta-search --predictor \
          --gen-spec-lib --out report.tsv --out-lib lib.tsv --threads 4 --qvalue 0.01 \
          --met-excision --cut 'K*,R*' --missed-cleavages 1 --min-pep-len 7 \
          --max-pep-len 30 --var-mods 1 --unimod4 --reanalyse --smart-profiling

DIA-NN 1.8.1, library-free against reviewed E. coli (4,531 entries), 40 minutes
on 4 cores:

| | at 1% FDR |
|---|---:|
| precursors | 9,560 |
| peptides (stripped sequences) | 7,671 |
| protein groups | 1,523 |
| proteins | 1,333 |

For scale, the OpenSWATH library the channel-2 test measures against supplies
4,710 peptides in these windows: a current DIA engine identifies more peptides
from this one run without any library than the library contains. DIA-NN's own
diagnostics on the run: median MS1 mass accuracy 4.13 ppm (1.41 corrected),
median MS2 6.03 ppm (2.44 corrected), FWHM 3.02 scans, and a recommended MS1
setting of 12.42 ppm, which is the closest thing it reports to the 13.39 ppm
scan-to-scan scatter measured above. The two are not the same quantity --
DIA-NN's figure is the accuracy of matched fragments after its own calibration,
the 13.39 ppm is the reproducibility of a peak between two scans of one window,
which is what a merge tolerance has to cover -- but they agree on the order.

Read against this, channel 2 is not a competitor. Its best operating point,
6.8% of library peptides beating all 110 decoys, is a weak feature on peptides
whose identity is already known, while DIA-NN reconstructs 7,671 peptides from
the same file with no library at all. The measurement stands as a statement
about the evidence channel, not as an identification method.

Two things had to be fixed before DIA-NN would run at all, both worth knowing
for anyone repeating this:

- The deposited mzXML carries `<precursorMz>` without a `windowWideness`
  attribute, so converting it to mzML yields isolation window *targets* with no
  width. DIA-NN then has no window to search in and reports **0 precursors**,
  exiting 0 as if it had succeeded.
- The scheme was reconstructed from the targets themselves (34 windows, 25.0 Th
  spacing), `isolation window lower/upper offset` of 12.5 Th written into all
  67,830 isolation windows, and the file re-indexed with msconvert.

---

## Appendix: DDA measurements on PXD004863

Recorded for completeness. These were run before it was decided that PXD004863
is not a fair basis for this question: it is an endopeptidome study, the file
used here is fraction 9 of 12, and identification of its spectra needs a
no-enzyme search. Read them as a comparison of merges on one fractionated
peptidomics run, not as evidence about either merge in general.

    cd analysis/dda
    python3 range.py                                  # the pipeline's own consensus file
    python3 consensus_merge.py --max-peaks 160
    python3 consensus_merge.py --max-peaks 0
    python3 range.py consensus_thomson.ms2
    python3 range.py consensus_ppm2_max160.ms2
    python3 range.py consensus_ppm2_max0.ms2
    python3 comp_pairs.py consensus_thomson.ms2
    python3 comp_pairs.py consensus_ppm2_max0.ms2

`consensus_merge.py` clusters by identity: 1,329 peptide/charge groups over
both replicates, 564 clusters of at least 2 spectra, 1,510 spectra in them.
These are not MaRaCluster's clusters. Peaks per consensus spectrum: Thomson
median 96, ppm median 121.

Complementary pairs, target beats all 20 decoys, null 4.76%:

| consensus | spectra | 0.005 Da beats-all | ratio | 0.01 Da beats-all | ratio |
|---|---:|---:|---:|---:|---:|
| identity clusters, Thomson | 564 | 59.57% | 26.24 | 59.57% | 17.08 |
| identity clusters, ppm cap 160 | 564 | 62.23% | 23.20 | 62.23% | 13.52 |
| identity clusters, ppm no cap | 564 | 62.41% | 20.56 | 62.23% | 13.52 |
| pipeline clusters, Thomson | 14,515 | 9.32% | 5.28 | — | — |
| pipeline clusters, ppm no cap | 14,515 | 10.07% | 3.17 | — | — |

The ppm merge raises the beats-all rate by about 3 points on identity clusters
and lowers the target/decoy ratio, because its extra peaks give the decoys more
to match as well.

`range.py`, fraction of spectra whose highest peak exceeds S*/2, by charge:

| consensus | z=2 | z=3 | z=4 |
|---|---:|---:|---:|
| pipeline clusters, Thomson | 1.00 | 0.88 | 0.40 |
| pipeline clusters, ppm no cap | 1.00 | 0.88 | 0.42 |
| identity clusters, Thomson | 1.00 | 1.00 | 1.00 |
| identity clusters, ppm cap 160 | 1.00 | 1.00 | 1.00 |
| identity clusters, ppm no cap | 1.00 | 1.00 | 1.00 |

The charge-4 shortfall the caveat rests on is unchanged by the merge (0.40 to
0.42) and absent entirely from identity clusters. Whatever causes it, it is not
the merge: the two populations differ (median S* at charge 4 is 3,469 on
MaRaCluster's clusters against 2,597 on identity ones), and only identified
spectra enter the identity clusters.

The pipeline-faithful ppm consensus file used above was produced by rerunning
Quandenser's consensus step over the two fraction-9 mzML files with
`--consensus-method ppm --consensus-max-peaks 0`, reusing the existing
clusters, so both rows rest on one clustering.
