# Calibration of the complementary-pair statistic (PXD004863, fraction 9)

Question: can peak pairs satisfying m(b_i) + m(y_{n-i}) = M + 2*proton identify
which fragments in a chimeric MS2 spectrum belong to a precursor of mass M,
without a spectral library or a search engine?

Data: the two fraction-9 hTau CSF runs, high-resolution MGF as deposited, with
PEAKS rank-1 identifications from the study (2,275 identified spectra of 36,592).
Null: shift M by an integer number of Da, which preserves the nominal-mass phase
of the pair-sum distribution. 110 decoy masses, so the null "beats all" rate is
1/111 = 0.90%.

Scripts: comp_pairs.py (consensus, superseded), comp_pairs_mzml.py (all MS2),
identified.py (known sequences), chimera.py (simulated co-isolation),
tolscan.py (mass-accuracy dependence), compare_mgf_mzml.py (control).

## Conclusion

The relation is real and cleanly detectable on isolated spectra, and it decays
faster than the co-isolation depth of real DIA. It is a usable Percolator
feature up to roughly 4 to 8 co-isolated peptides, not a grouping mechanism.
Tighter MS2 mass accuracy does not rescue it.

---

# Channel 1: separation of co-eluting species by profile correlation

Fragment co-variation cannot be measured on DDA data, since a DDA run fragments
each precursor once rather than across its elution. But the binding constraint
on channel 1 is whether elution profiles of distinct co-eluting species are
separable at all at DIA sampling rates, and that is measurable on MS1, where
ground truth is free.

  positive = monoisotopic vs M+1 XIC of one Dinosaur feature (same molecule)
  negative = monoisotopic XICs of two different features that co-elute and fall
             in the same simulated 25 m/z isolation window

Both over the candidate feature's own elution window. Same-molecule pairs
(same peptide at two charge states, isotope confusions) are excluded from the
negative set by neutral mass; leaving them in inflates the interference tail
and destroys the measurement.

The positive control is optimistic: mono and M+1 lie 1/z apart and share their
local interference, while two fragments of one peptide are far apart in m/z and
do not. Every number here is an upper bound on real fragment pairs.

Scripts: xic.py, corr3.py (single run), crossrun2.py (two runs).

## Result

Discrimination is dominated by trace intensity, not by sampling rate. Going
from 89 to 13 points per peak costs about 4 AUC points; going from the top
intensity quartile to the bottom costs about 20.

Pooling two aligned runs gives a consistent gain, and at a 7 s cycle a single
run has fewer than 5 points across a median peak, so the correlation cannot be
computed at all from one run.

The retention-time warp between the two replicates, acquired back to back, is
about 290 s and varies with retention time from -324 to -238 s. Without
correcting it only 7.9% of mass-matched features pair within 45 s. That is the
concrete cost of not having an alignment, and the concrete argument for doing
this inside Quandenser.

---

# Real DIA data: PXD001587 (MassIVE-hosted PXD026600 was unreachable)

SWATH gold-standard E. coli, TripleTOF, 34 windows of 25 m/z, matched DDA runs.
Ground truth from the study's own artefacts: a SpectraST library built from the
DDA runs (6,409 entries with annotated b/y ions, median 14 per peptide) and
OpenSWATH/mProphet peak groups giving each peptide's retention time in the DIA
run itself, so no alignment step is needed. 4,711 peak groups at m_score < 0.01.

Scripts in /home/user/dia: scores.py, combine.py, diatest.py, analyse.py,
ch2_sweep.py, indep.py.

## Channel 1

Fragment MS2 trace against the peptide's own MS1 precursor trace, true library
ions versus the library ions of a different peptide in the same window eluting
more than 600 s away: AUC 0.741, TPR 37.4% at 5% FPR, 19.5% at 1% FPR.

The MS1-isotope proxy measured on the CSF data predicted AUC 0.739. The proxy
was well calibrated.

OpenSWATH's own five co-elution sub-scores combined give AUC 0.826, which also
matches the CSF two-run estimate of 0.829.

## Channel 2

Dead on raw merged window spectra (~3200 peaks): ratio 1.00, below null. With
denoising to the top 50-100 peaks it recovers 4x enrichment over null, which is
what the CSF chimera simulation predicted for this co-isolation depth (3.7x at
k=8). Still only a 3.6% hit rate, so a weak feature, not a mechanism.

## Independence: no

Channel 2's hits concentrate where channel 1 is already strong.

    split by OpenSWATH xcorr_shape   weak 2.0x   middle 2.9x   strong 7.1x
    split by own co-elution score    weak 0.7x   middle 1.5x   strong 9.9x

Rank correlation +0.17 to +0.35. Both channels depend on the peptide being
abundant and well fragmented, so they fail on the same spectra. The same holds
for OpenSWATH's library-based identity features: co-elution alone AUC 0.826,
identity alone 0.946, both together 0.948, i.e. the combination adds 0.002.

## Cross-run pooling on real DIA: much smaller than the CSF proxy suggested

Both SWATH replicates, 4,384 peptides confident in both, decoy partners assigned
from the library by a stable hash so the same peptide gets the same decoy
fragments in each run (drawing them per run instead drops 94% of decoys on
pooling and makes the two rows incomparable).

    single run        AUC 0.740   TPR@5% 37.5%   TPR@1% 19.4%
    two runs pooled   AUC 0.753   TPR@5% 41.0%   TPR@1% 21.9%

+0.013 AUC, against the +0.09 the CSF MS1 proxy predicted. The proxy was right
about the single-run number and wrong about the pooling gain, because MS1
isotope pairs are limited by random noise, which averaging reduces, while real
fragment traces are limited by interference from co-eluting peptides, which is a
property of the sample and therefore reproduces in every replicate. Averaging
replicates does not remove it.

---

# Quant first: linking unidentified DIA features across runs

The identification framing failed because a library-free pseudo-spectrum has to
be correct. Linking only requires it to be reproducible, and the thing that
killed cross-run pooling above, that interference is a property of the sample
and so repeats in every run, is exactly what makes a signature reproducible.

Method, entirely library-free. Candidate fragments are the 150 most intense
peaks of each (isolation window, 60 s block) merged spectrum, so every feature
in a block is scored against the same pool and all specificity has to come from
the filter. A fragment joins a feature's signature if its MS2 trace correlates
above 0.5 with that feature's own MS1 trace; at most 20 fragments are kept.
MS1 features come from peak picking on the MS1 survey scans, 97k per run.
Identity is used only to label the right answer.

Scripts in /home/user/dia: pseudospec.py, linktest.py, linkstrat.py, ms1feat.py,
blockspec.py, linkfull.py, linkreal.py, linkvalue.py, yield.py, redundancy.py.

## Signatures are reproducible and specific

Same peptide across runs, median similarity 0.540. A different feature
co-eluting in the same 25 m/z window, 0.000. The correlation filter makes a
signature specific to its own feature even though the candidate pool is shared.

## Link validation

True cross-run pair against a same-m/z feature more than 120 s away, which is
the false link that m/z alone accepts:

    AUC 0.925   TPR 82.8% at 5% FPR   76.5% at 1% FPR
    weak quartile 73.0% above the 5% FPR threshold, strong 91.2%

Compare the identification framing on the same data: per-fragment AUC 0.740.

## Candidate selection needs alignment error to be worth it

On these two replicates, acquired back to back with a 2.5 s offset, m/z plus
nearest retention time is already 95.3% correct and the signature does not beat
it. Adding a retention-time perturbation to stand for alignment residual:

    residual   0s  10s  20s  40s  60s
    RT-nearest 91.2 88.5 79.6 65.1 50.5 %
    signature  85.9 85.8 86.1 86.5 86.6 %

Crossover near 20 s. The CSF replicates in this session, also run back to back,
differ by 290 s varying with retention time, so a residual in that range is not
exotic, but it has to be demonstrated on runs that actually have it.

## Yield

At a false-link rate calibrated on the null, from a 30,000-feature sample:

    identified     789 / 1552   linked (50.8%)
    unidentified  6407 / 28448  linked (22.5%)

The gap is entirely the intensity distribution. Within a decile the two are the
same: decile 8, 32.4% against 33.8%; decile 9, 43.7% against 44.2%; decile 10,
61.8% against 66.6%. Being unidentified is not a handicap for linking.

Collapsing isotope peaks and split traces, which this crude feature detector
counts separately, 5,124 distinct unidentified signals against 773 identified,
a ratio of 6.6x rather than the 8.1x at feature level.

## Caveats

Two technical replicates from one day; the alignment-residual row is a
simulation on top of them. Every query with known truth is a peptide that was
identified, so the unidentified yield rests on a decoy-calibrated null rather
than on truth. The feature detector has no isotope grouping or charge
deconvolution; Dinosaur would change the feature set.
