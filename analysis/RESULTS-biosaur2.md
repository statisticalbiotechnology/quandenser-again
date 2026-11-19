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
about 59,000 PASEF precursors each. The method is a short one
(`RAH_short_method`): each run lasts 11.5 min, and 90% of features have
their apex between 3.2 and 11.2 min. The Orbitrap runs above span 12–131 min. Features came from biosaur2 caa800e;
Quandenser ran with `--max_missing 2 --target_search_threshold 1.0`. The
pipeline finished, but two stages that work on Orbitrap fail here.

### Match-between-runs finds almost nothing

| | features / run | links at q < 0.02, first run pairs | feature-group rows |
|---|---:|---:|---:|
| PXD002370, Dinosaur | 104k | 10,000–78,000 | 1,175,909 |
| PXD039875, biosaur2 | 600k–750k | 0–1,201 | 7,878 |

biosaur2 reports six times more features per run on PASEF data than either
detector does on Orbitrap, packed into an eight-minute window. In run 4324,
42% of them span two MS1 frames, and 94% are narrower than 5 s.

RT alignment itself works: every run pair has 15,600–25,900 alignment pairs
and a combined RMSE of 0.17–0.24 min. Quandenser searches for a feature's
partner within `--align_rtime_tol` (10) times that RMSE, 2–4 min, and within
`--align_mz_tol` (20 ppm). On an Orbitrap run that window is a small slice
of a two-hour gradient. Here it is a third of the run. Decoy candidates (the
same window at m/z + 5.0025) are then almost as common as real ones.

The table shows the share of features in one run that have a candidate
partner in the next, from `ft_density.py`. It was run on 4324 → 4342
(biosaur2 03ed8d5) and on PXD002370 C-R1 → C-R2 (Dinosaur).

| tolerance | Orbitrap: target | decoy | PASEF: target | decoy |
|---|---:|---:|---:|---:|
| 20 ppm, ±3 min | 78% | 24% | 95% | 93% |
| 10 ppm, ±1 min | 70% | 5% | 76% | 72% |
| 10 ppm, ±1 min, 1/K0 ±0.03 | | | 47% | 39% |

Tighter m/z and RT tolerances do not fix it. Neither does a feature filter.
Keeping the 50,000 most intense features per run, or those spanning eight
or more frames, still leaves 55–66% decoy candidates at the default
tolerances. Only ion mobility combined with tight tolerances separates true
from random partners, and Quandenser's matching does not read ion mobility.

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

### Rerun with gap bridging and p ≤ 1e-5

The same six runs with biosaur2 03ed8d5, `-hill_max_missing 1` and
`--maracluster_pvalue_cutoff -5` (`/home/lukask/data/PXD039875/chain7.sh`):

| | first run | rerun |
|---|---:|---:|
| features / run | 597k–754k | 526k–675k |
| precursors of run 4324 inside a feature | 52.1% | 54.7% |
| spectra in singleton clusters | 77.2% | 55.6% |
| same-peptide pairs clustered (impurity) | 25.2% (0.45%) | 35.5% (0.75%) |
| consensus spectra | 21,554 | 32,742 |
| links at q < 0.02, first run pairs | 0–1,201 | 0–613 |
| feature-group rows | 7,878 | 12,599 |

Clustering and consensus spectra improve as the threshold sweep predicted.
Inside Quandenser the gain is smaller than for MaRaCluster alone (35.5% vs
46.4%). Match-between-runs stays broken, for the reasons above.

### Ion mobility in feature matching

`--align_im_tol` (commit e10e222) makes Quandenser require a feature's
partner to lie within that 1/K0 tolerance and gives Percolator the 1/K0
difference as a feature. Two runs (4324 → 4342, biosaur2 03ed8d5,
p ≤ 1e-5), first matching direction:

| settings | candidate pairs, target : decoy | links at q < 0.02 |
|---|---:|---:|
| 20 ppm, 10 RMSE, no ion mobility (old and new build) | 669,716 : 596,004 | 0 |
| 20 ppm, 10 RMSE, 1/K0 ± 0.03 | 534,630 : 459,294 | 0 |
| 10 ppm, 3 RMSE, 1/K0 ± 0.03 | 312,957 : 250,509 | 0 |

With the option unset, the new build's feature groups, link PSMs and
consensus spectra are byte-identical to the old build's. On these
11.5-min runs, ion mobility does not rescue matching: random partners
remain almost as common as real ones, as the density table predicted.

### Next

1. Test on a timsTOF dataset with a gradient of 30 min or more, with and
   without `--align_im_tol`. PXD028735 is running: an HYE benchmark with
   six runs, a 120-min gradient, and 370k PASEF spectra per run. Its runs
   need `--biosaur2_chunk_minutes` (commit 415ecb3) to fit in memory.
2. Check charge handling in the PASEF MS2 mzML. Tide had to infer charge
   for some spectra, and its 6% identification rate is far below what
   plasma dda-PASEF usually gives. This only affects how many spectra the
   clustering scores rest on, not the comparison itself.
