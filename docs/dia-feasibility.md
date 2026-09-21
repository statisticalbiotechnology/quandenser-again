# Quandenser and DIA: what was tested and what came of it

This is a record of a feasibility study, not a design document. It exists
because the conclusions are mostly negative and the negative results cost real
work to obtain. Anyone proposing to extend Quandenser to data-independent
acquisition should read it before starting, and should in particular not
re-derive the parts that were measured and found wanting.

All numbers below come from public data. The scripts are in `analysis/`.

## The problem

Quandenser quantifies with match-between-runs: MaRaCluster clusters MS2 spectra
to establish that two features in different runs are the same peptide, the
runs are aligned in retention time, and features are linked with a posterior
error probability. DIA breaks exactly one link in that chain. DIA runs still
acquire MS1 survey scans, so feature detection, alignment and feature linking
work unchanged. What fails is the assumption behind MaRaCluster, that a
spectrum is one peptide. In a DIA isolation window an MS2 spectrum is a
superposition of everything co-eluting in that window, so two runs of the same
peptide share only the subset of peaks belonging to it, diluted by a different
mixture each time.

The question was whether the clustering step can be replaced by something that
works on chimeric spectra, without a spectral library.

## What was tried

Two sources of evidence, referred to below as channels.

**Channel 1, co-elution.** Fragments of one peptide co-vary across the elution
peak. Well-trodden ground: DIA-Umpire, CorrDec, a PARAFAC formulation on the
(m/z, retention time, sample) tensor, and Group-DIA for the cross-run case.

**Channel 2, complementary pairs.** For a peptide of neutral monoisotopic mass
M, singly charged complementary ions satisfy

    m(b_i) + m(y_{n-i}) = M + 2 * proton

independent of i. Counting peak pairs summing to that value tests a candidate
precursor mass against a spectrum, with no library, no search engine and no
predicted spectra. The null comes free: shift M by an integer number of Da and
recount. Integer shifts matter, because the pair-sum distribution is periodic
at 1 Da (fragment masses follow the averagine line) and a non-integer shift
compares against a different phase, manufacturing a signal.

The plan was to combine the two as features in a Percolator model, on the
reasoning that they fail on different spectra. That reasoning was wrong, and
establishing that was most of the work.

## Channel 2

### Calibration on DDA, where the truth is known

PXD004863, fraction 9 of the hTau CSF series, high-resolution MGF as deposited,
with the study's own PEAKS rank-1 identifications. 2,275 identified spectra of
36,592. 110 decoy masses, so the null "beats all" rate is 0.90%.

The relation is real and cleanly detectable on isolated spectra:

| z | n | pep len | b seen | y seen | true b/y pairs | found | decoy |
|---|---|---|---|---|---|---|---|
| 2 | 1097 | 12.3 | 3.3 | 4.3 | 1.75 | 1.87 | 0.07 |
| 3 | 797 | 18.4 | 3.0 | 5.5 | 1.27 | 1.32 | 0.03 |
| 4 | 310 | 24.9 | 2.8 | 5.9 | 0.96 | 0.97 | 0.03 |
| all | 2275 | 16.8 | 3.2 | 5.0 | 1.43 | 1.51 | 0.05 |

77.4% of identified spectra contain at least one genuine complementary pair.
The blind count tracks the annotated count almost exactly, 1.51 against 1.43
with a decoy background of 0.05, so the statistic measures what it is supposed
to and nothing else. Enrichment over decoy is 22-fold.

Best tolerance is 0.005 to 0.01 Da, which incidentally confirms the MS2 data is
high resolution without having to trust the metadata.

### It does not survive multiplexing

Pooling each spectrum with k-1 randomly drawn spectra from the same run:

| co-isolated | peaks | target | decoy | excess | Poisson SNR | beats-all | vs null |
|---|---|---|---|---|---|---|---|
| 1 | 71 | 1.51 | 0.07 | 1.44 | 5.4 | 32.0% | 36x |
| 2 | 138 | 1.66 | 0.18 | 1.48 | 3.5 | 18.7% | 21x |
| 4 | 270 | 2.16 | 0.55 | 1.61 | 2.2 | 8.0% | 8.9x |
| 8 | 535 | 3.71 | 1.87 | 1.84 | 1.3 | 3.3% | 3.7x |
| 16 | 1062 | 8.92 | 6.92 | 2.00 | 0.8 | 1.8% | 2.0x |
| 32 | 2121 | 29.51 | 26.63 | 2.88 | 0.6 | 1.5% | 1.7x |

The excess is essentially constant, as it must be: the peptide's own pairs do
not disappear when clutter is added. The background grows as the square of the
peak count. Signal-to-noise crosses 1 between 8 and 16 co-isolated peptides,
which is the regime of a 25 m/z window on a real digest.

Tightening mass accuracy does not rescue it. The optimum sits at 0.005 to
0.01 Da, and below that true pairs are lost faster than background is removed,
because the error on a sum is the quadrature sum of two peak errors:

| beats-all | ±0.002 Da | ±0.005 Da | ±0.01 Da | ±0.02 Da |
|---|---|---|---|---|
| k=1 | 25.5% | 34.4% | 32.0% | 25.9% |
| k=8 | 2.4% | 3.3% | 3.4% | 3.5% |
| k=16 | 1.5% | 1.9% | 2.2% | 2.3% |

A better mass spectrometer does not buy this.

### Confirmed on real DIA

On raw merged window spectra from PXD001587, about 3,200 peaks, the statistic
is dead: ratio 1.00, hit rate below null. That is not a fair test, since the
DDA calibration used 65-peak denoised spectra. Denoised:

| top-N peaks | target | decoy | ratio | beats-all | enrichment |
|---|---|---|---|---|---|
| 50 | 0.87 | 0.46 | 1.89 | 3.9% | 4.3x |
| 100 | 1.29 | 0.74 | 1.76 | 3.6% | 4.0x |
| 200 | 3.78 | 2.84 | 1.33 | 2.4% | 2.7x |
| 500 | 17.5 | 16.3 | 1.08 | 1.0% | 1.1x |
| all | 442.9 | 442.0 | 1.00 | 0.2% | 0.3x |

The simulation predicted 3.7x at eight co-isolated peptides; the real data gives
4.0x. A weak feature, not a mechanism.

## Channel 1

### Measured on DIA fragment traces

PXD001587, SWATH gold-standard E. coli, TripleTOF, 34 windows of 25 m/z.
Ground truth from the study's own artefacts: a SpectraST library built from its
DDA runs (6,409 entries with annotated b/y ions, median 14 per peptide) and
OpenSWATH/mProphet peak groups giving each peptide's retention time in the DIA
run itself, so no alignment step enters the measurement.

Each library fragment's MS2 trace against the peptide's own MS1 precursor
trace. Decoys are the library ions of a different peptide in the same window
eluting more than 600 s away.

    AUC 0.740   TPR 37.5% at 5% FPR   19.4% at 1% FPR

OpenSWATH's own five co-elution sub-scores combined give AUC 0.826.

A proxy measured on MS1 isotope pairs in the CSF DDA data predicted AUC 0.739,
and was well calibrated including the direction of its bias. Two things that
proxy got right and one it did not are recorded below, because the one it got
wrong is the interesting one.

What it got right: the single-run number, and that discrimination is set by
trace intensity rather than by sampling rate. Going from 89 points per peak to
6 costs about 4 AUC points; going from the top intensity quartile to the bottom
costs about 20. Cycle time in the 2 to 7 s range is not the binding constraint.

### Cross-run pooling is worth less than the noise argument suggests

Both SWATH replicates, 4,384 peptides confident in both, decoy partners
assigned from the library by a stable hash so the same peptide gets the same
decoy fragments in each run. Drawing them per run instead drops 94% of decoys
on pooling and makes the two rows incomparable, which is a trap worth naming.

| | AUC | TPR@5% | TPR@1% |
|---|---|---|---|
| single run | 0.740 | 37.5% | 19.4% |
| two runs pooled | 0.753 | 41.0% | 21.9% |

+0.013, against the +0.09 the MS1 proxy predicted. MS1 isotope pairs are
limited by random noise, which averaging reduces. Real fragment traces are
limited by interference from co-eluting peptides, which is a property of the
sample and therefore reproduces in every replicate. Averaging replicates does
not remove it, so the sqrt(N) argument does not apply to the dominant error
term.

## The channels are not independent

This was the load-bearing assumption of the combination plan, and it is false.

Channel 2's hits concentrate where channel 1 is already strong:

| split by | weak | middle | strong |
|---|---|---|---|
| OpenSWATH `xcorr_shape` | 2.0x | 2.9x | 7.1x |
| own co-elution score | 0.7x | 1.5x | 9.9x |

Rank correlation +0.17 to +0.35. In the weakest co-elution tercile channel 2 is
at or below null. Both channels need the peptide to be abundant and well
fragmented, so they fail on the same spectra.

OpenSWATH's own sub-scores say the same, trained with peptide-grouped
cross-validation on its target and decoy peak groups:

| feature set | AUC | TPR@1%FPR |
|---|---|---|
| co-elution only | 0.826 | 29.4% |
| fragment identity only | 0.946 | 68.1% |
| both | 0.948 | 68.0% |
| everything | 0.953 | 71.3% |

Adding co-elution to identity buys 0.002 AUC. The two are correlated at
rho = -0.81 among targets and -0.60 among decoys.

**Conclusion for the identification framing.** A library-free decomposition
built on these two channels will not reach a usable operating point. The
information in DIA lives in the fragment-identity channel, and at AUC 0.946 the
library-based route already extracts it. That is what DIA-NN and Spectronaut
do, increasingly with predicted rather than measured libraries.

## The reframing that works

Quandenser's purpose is quantification, and the features that matter most are
the ones no search engine identifies. For those, the spectrum does not have to
be correct. It has to be reproducible. And the finding that killed cross-run
pooling, that interference is a property of the sample and repeats in every
run, is exactly what makes a library-free signature reproducible.

### Method

Entirely library-free. Candidate fragments are the 150 most intense peaks of
each (isolation window, 60 s block) merged spectrum, so every feature in a
block is scored against the same pool and all specificity has to come from the
filter. A fragment joins a feature's signature if its MS2 trace correlates
above 0.5 with that feature's own MS1 trace; at most 20 are kept. MS1 features
come from peak picking on the survey scans, 97k per run. Identity is used only
to label the right answer.

### Signatures are specific

Same peptide across runs, median similarity 0.540. A different feature
co-eluting in the same 25 m/z window, 0.000. The correlation filter separates
two features sharing a window and a retention time even though they draw
candidates from one pool.

### Link validation

True cross-run pair against a same-m/z feature more than 120 s away, which is
the false link that m/z alone accepts:

    AUC 0.925   TPR 82.8% at 5% FPR   76.5% at 1% FPR
    weak quartile 73.0% above the 5% FPR threshold, strong 91.2%

The same data gave AUC 0.740 per fragment under the identification framing.
This is the linkPEP analogue, obtained without a targeted search and without an
identification.

### Candidate selection needs real alignment error

On two replicates acquired back to back with a 2.5 s offset, m/z plus nearest
retention time is already 95.3% correct and the signature does not beat it.
Adding a perturbation to stand for alignment residual:

| residual | 0 s | 10 s | 20 s | 40 s | 60 s |
|---|---|---|---|---|---|
| RT-nearest | 91.2% | 88.5% | 79.6% | 65.1% | 50.5% |
| signature | 85.9% | 85.8% | 86.1% | 86.5% | 86.6% |

Crossover near 20 s. The signature is flat because it does not use retention
time. Two CSF replicates run back to back on the same instrument differ by
about 290 s varying with retention time, so a residual in that range is not
exotic, but it has to be shown on runs that actually have it.

### Yield

At a false-link rate calibrated on the null, from a 30,000-feature sample:

| | linked | of | rate |
|---|---|---|---|
| identified | 789 | 1,552 | 50.8% |
| unidentified | 6,407 | 28,448 | 22.5% |

The gap is entirely the intensity distribution. Within a decile the two are
indistinguishable: decile 8, 32.4% against 33.8%; decile 9, 43.7% against
44.2%; decile 10, 61.8% against 66.6%. Being unidentified is not a handicap for
linking. Identification failure and linking failure are different failures,
which is the premise of quantification-first and is here measured rather than
assumed.

Collapsing isotope peaks and split traces, which this crude feature detector
counts separately, 5,124 distinct unidentified signals against 773 identified,
a ratio of 6.6x rather than the 8.1x at feature level.

## Proposed architecture

Most of it already exists.

1. MS1 feature detection per run. Dinosaur, as in the DDA path.
2. Per (isolation window, retention-time block), a merged spectrum and traces
   for its most intense peaks. New, and cheap: two streaming passes, about a
   minute per run on a 1.8 GB file.
3. Per feature, a signature from correlating those traces against the feature's
   own MS1 trace. New, and cheap.
4. Retention-time alignment across runs. Already in Quandenser.
5. Feature linking, with the signature similarity supplying the linkPEP.
   Replaces the MaRaCluster step, which cannot run on chimeric spectra.
6. Triqler downstream, unchanged.

Identifications from DIA-NN, Spectronaut, OpenSWATH or anything else attach
afterwards to whichever clusters have them. The pipeline does not need them to
produce a quantification matrix, which is the point.

## Caveats

- Two technical replicates from one day. The alignment-residual row is a
  simulation on top of them, not a measurement of them.
- Every query with known truth is a peptide that was identified. The
  unidentified yield rests on a decoy-calibrated null rather than on truth.
- The MS1 feature detection used here has no isotope grouping and no charge
  deconvolution. Dinosaur would change the feature set; the obvious redundancy
  was collapsed after the fact rather than avoided.
- Consensus spectra from MaRaCluster are binned and cannot support a 0.01 Da
  test. Any fragment-level work must use the original spectra. 60% of charge-4
  consensus spectra have no peak above half the precursor-pair sum, so no
  complementary pair is constructible in them at all.

## What to test next

Runs with genuine retention-time drift, which is the one condition under which
the signature was not shown to beat nearest-retention-time. A multi-day DIA
series, or several batches of the same cohort. That single experiment decides
whether step 5 above is worth building.

## Data

| Accession | Use |
|---|---|
| PXD004863 | hTau CSF, DDA, Orbitrap Fusion. Channel 2 calibration and the multiplexing simulation. High-resolution MGF and PEAKS identifications as deposited. |
| PXD001587 | SWATH gold standard, E. coli, TripleTOF. Matched DDA and DIA as mzXML, with a SpectraST library and OpenSWATH peak groups inside `Ecoli.zip`. |

PXD026600 (Orbitrap DIA with narrow and wide window schemes) would be a better
fit for the window-width question but is hosted on MassIVE, which was
unreachable from this environment.
