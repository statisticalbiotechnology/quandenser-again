# biosaur2 as Quandenser's feature detector, measured

Question: can biosaur2 replace Dinosaur in the pipeline, on Orbitrap without
losses and on timsTOF where Dinosaur cannot go?

Scripts are in `benchmark/`. Data and search outputs live outside the
repository, under `/home/lukask/data/PXD002370` and
`/home/lukask/data/PXD039875`.

---

## Orbitrap: PXD002370 (UPS1 in yeast, 9 runs)

UPS1 spiked at 25, 10 and 5 fmol (conditions C, D, E; three runs each) into a
yeast background. True log2 ratios are C/E 2.32, D/E 1.00 and yeast 0. Each
variant ran the whole pipeline. Its consensus spectra were searched with Crux
Tide and Percolator (`benchmark/search.sh`; UPS1 + reviewed yeast, 20 ppm,
0.02 Th bins). Quantities come from `benchmark/ups_ratios.py`: the mean log2
intensity per condition of each identified feature group (q ≤ 0.01).

| feature detector | features / run | consensus spectra | PSMs at 1% | id. feature groups | yeast C/E sd | yeast D/E sd | UPS C/E error | UPS D/E error |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Dinosaur | 104,054 | 101,540 | 25,761 | 11,954 | 0.32 | 0.22 | +0.56 | +0.13 |
| biosaur2 2e7c9a1 (released) | 114,761 | 91,864 | 24,511 | 11,325 | 0.43 | 0.34 | +0.30 | −0.01 |
| fix/known-issues caa800e | 112,486 | 91,657 | 24,564 | 11,323 | 0.40 | 0.31 | +0.32 | +0.00 |
| caa800e, `-iuse -1` | 112,486 | 91,387 | 24,571 | 11,307 | 0.36 | 0.27 | +0.58 | +0.14 |
| feat/hill-gaps 03ed8d5, `-hill_max_missing 1` | 72,430 | 93,285 | 24,507 | 11,394 | 0.30 | 0.21 | +0.30 | +0.02 |

"sd" is the robust standard deviation of yeast log2 ratios, so lower is more
precise. "error" is the median UPS1 log2 ratio minus the true ratio, so closer
to 0 is more accurate.

* **Gap bridging makes biosaur2's quantification as good as Dinosaur's or
  better.** Released biosaur2 ends a hill at the first scan in which the peak
  is missing. Its features were a median 3 s wide, against Dinosaur's 14 s.
  Letting a hill bridge one missing scan widens them to 17 s, and the yeast
  spread drops to Dinosaur's level (0.21 vs 0.22). UPS1 ratios keep
  biosaur2's smaller compression (+0.30 vs +0.56).
* **Summing isotopes (`-iuse -1`) is the wrong lever.** It recovers part of
  the precision but gives up the accuracy advantage.
* **Identifications stay about 5% below Dinosaur for every biosaur2
  variant.** Wider features raised the share of MS2 precursors inside a
  feature (`benchmark/precursor_match.py`) from 47.9% to 50.1%, against
  52.9% for Dinosaur. PSMs did not move. The gap is in the consensus
  spectra (93k vs 102k), and the next place to look is how many precursor
  candidates each MS2 spectrum gets: 1.3–1.5 mass-charge states per spectrum
  with biosaur2 against 1.8 with Dinosaur.
* **The monoisotopic rescue is not worth it.** Moving features that biosaur2
  places one isotope high shifted 615 of 119,227 features and gained
  0.1 points of precursor coverage (branch `exp/mono-rescue`, experiment
  only).

---

## timsTOF: PXD039875 (plasma, dda-PASEF)

Six plasma runs on the same gradient (4324, 4342, 4353, 4403, 4448, 4471),
about 59,000 PASEF precursors each. Features came from biosaur2 caa800e;
Quandenser ran with `--max_missing 2 --target_search_threshold 1.0`. The
pipeline finished, but two stages that work on Orbitrap fail here.

### Match-between-runs finds almost nothing

| | features / run | links at q < 0.02, first run pairs | feature-group rows |
|---|---:|---:|---:|
| PXD002370, Dinosaur | 104k | 10,000–78,000 | 1,175,909 |
| PXD039875, biosaur2 | 600k–750k | 0–1,201 | 7,878 |

biosaur2 reports six times more features per run on PASEF data than either
detector does on Orbitrap. In run 4324, 42% of them span two MS1 frames,
and 94% are narrower than 5 s. A peptide eluting over 10–20 s is cut into several short
features. Quandenser's feature-to-feature matching then cannot separate true
partners from decoys. RT alignment itself works: every run pair has
15,600–25,900 alignment pairs and a combined RMSE of 0.17–0.24 min.

### MaRaCluster clusters PASEF spectra much less readily

Quality is scored with identifications of the raw spectra
(`benchmark/cluster_ids.py`). Over pairs of spectra identified as the same
peptide and charge, it reports the share that MaRaCluster puts in one
cluster ("sensitivity"). Over identified pairs that share a cluster, it
reports the share carrying different peptides ("impurity").

PASEF identifications come from Tide + Percolator against reviewed human
proteins (20 ppm, 0.02 Th bins). The search was killed by the 16 GB container
limit in the fifth file, so the scores use the four runs searched in full:
14,159 identified spectra. The search and the threshold sweep are
`/home/lukask/data/PXD039875/pasef_search.sh`.

| clustering | threshold | Orbitrap: singletons | sensitivity | impurity | PASEF: singletons | sensitivity | impurity |
|---|---|---:|---:|---:|---:|---:|---:|
| inside Quandenser | 1e-10 | 16.8% | 54.0% | 7.0% | 77.2% | 25.2% | 0.45% |
| MaRaCluster alone | 1e-15 | | | | 87.1% | 21.8% | 0.02% |
| MaRaCluster alone | 1e-10 | 10.7% | 64.7% | 0.22% | 79.9% | 30.9% | 0.03% |
| MaRaCluster alone | 1e-7 | | | | 71.9% | 38.7% | 0.04% |
| MaRaCluster alone | 1e-5 | 7.4% | 79.3% | 0.65% | 63.0% | 46.4% | 0.04% |
| MaRaCluster alone | 1e-3 | 5.3% | 88.8% | 1.81% | 47.5% | 57.2% | 0.04% |

The Orbitrap identifications are the same kind of search on PXD002370's nine
raw runs (149,760 spectra at q ≤ 0.01, a 28% identification rate).

Inside Quandenser a spectrum enters one clustering per precursor candidate
(mass-charge state) it gets from the features. A co-isolated second peptide
therefore puts it in a cluster with that peptide's spectra. This is why
impurity reads higher there (7.0% on Orbitrap) than for MaRaCluster alone,
and it is not necessarily an error.

PASEF MS2 spectra have a median of 134 peaks, against 304 on the Orbitrap.
Under identical MaRaCluster settings, same-peptide PASEF spectra cluster at
half the Orbitrap rate (31% vs 65% at 1e-10). Standalone MaRaCluster does no
better than Quandenser on PASEF, so the feature assignment is not the cause:
the spectra rarely reach significance. Loosening the threshold to 1e-3
brings PASEF to roughly what the Orbitrap gets at the default 1e-10.
Impurity among identified PASEF spectra stays below 0.1%. That is a floor:
unidentified spectra are noisier and are not scored. On the Orbitrap,
impurity rises from 0.2% to 1.8% over the same range. The pipeline
already exposes the threshold as `--maracluster_pvalue_cutoff`.

---

## What would make timsTOF work

1. Run with `--maracluster_pvalue_cutoff -5`, or `-3`. It is a single flag;
   the numbers above say it brings same-peptide clustering from 31% to
   46–57%.
2. Use `-hill_max_missing 1` on PASEF data, plus a minimum feature length.
   The aim is a feature list closer to Orbitrap density (about 100k per
   run), in which feature-to-feature matching can tell true partners from
   decoys.
3. Check charge handling in the PASEF MS2 mzML. Tide had to infer charge for
   some spectra, and its 6% identification rate is far below what plasma
   dda-PASEF usually gives. This only affects how many spectra the scores
   above rest on, not the comparison itself.
