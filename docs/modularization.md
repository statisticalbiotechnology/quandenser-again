# Modularizing Quandenser

Notes from an analysis of the code base, and the reasoning behind the pipeline
in this repository. Line references are to `1aec399` (v0.03.2).

## 1. Why the tool was hard to run

Three problems, in order of severity.

### The build breaks itself over time

`ext/maracluster/admin/builders/install_proteowizard.sh` fetches ProteoWizard
from a floating pointer:

```
https://teamcity.labkey.org/guestAuth/repository/download/bt81/.lastSuccessful/VERSION
```

Whatever ProteoWizard's CI built most recently is what you get. At the time of
writing that is **3.0.26258**; the `ChangeLog` records v0.03.2 as released
against **3.0.22231**. Four years of drift, unrequested.

Line 60 of the same script is the other half:

```sh
rsync -ap --include "*/" --include '*.ipp' --exclude '*' libraries/boost_1_76_0/boost/ ../include/boost
```

The Boost directory name is hard-coded. The tarball now ships `boost_1_86_0`,
so the source directory does not exist, rsync exits 23, and because the script
has no `set -e` the failure is discarded. The build then continues for another
minute and dies with

```
fatal error: boost/date_time/gregorian_calendar.ipp: No such file or directory
fatal error: boost/unordered/unordered_map.hpp: No such file or directory
```

which has no visible connection to its cause. This was reproduced here, not
inferred: a from-source build of this repository fails today, on a clean
machine, for this reason.

A floating dependency pinned against a fixed version constant is a structural
guarantee that a project breaks itself with no commits. That is how software
with real users becomes software with none.

### The release channel is dead

`.github/workflows/build_and_release.yml` builds six platform targets and is
non-functional on several independent axes: `actions/upload-artifact@v1` and
`download-artifact@v3` have been shut down by GitHub; `centos:centos7` and
`fedora:35` no longer have package mirrors; the Maven URLs used by the Windows
and macOS builders are 404; `macos-latest` is arm64 while the script hard-codes
`x86_64`. So there are no fresh binaries either, and the last seven commits in
the repository are a single day spent losing a fight with this workflow.

### Even a successful build leaves sharp edges

- **Dinosaur asks for a 24 GB Java heap by default** (`src/DinosaurIO.cpp:26`).
  No environment variable can lower it: the `-Xmx` is built into the command
  line Quandenser constructs, so `JAVA_TOOL_OPTIONS` is ignored. Only
  `--dinosaur-memory` works. Under any scheduler that enforces a memory
  request, the default is an immediate OOM kill.
- **The Dinosaur jar path is an absolute compile-time constant**
  (`CommonCMake.txt:83` → `src/Globals.h.cmake:26` → `src/DinosaurIO.cpp:22-24`),
  with no environment or CLI override. A staged install (`make install
  DESTDIR=...`) silently produces a binary that cannot find its own jar.
- **The repository directory must literally be named `quandenser`**, because
  `quickbuild.sh` calls the builder without `-s` and the builder appends
  `/quandenser` to the parent directory. This checkout, named
  `quandenser-again`, cannot build with its own script.
- **Quandenser POSTs to Google Analytics on every run**
  (`src/Quandenser.cpp:463-488`), over plaintext port 80, with no opt-out. It
  is wrapped in `catch(...)`, so it fails silently, but on a sandboxed executor
  it can add a connect timeout to every invocation.

A container fixes the first three problems outright and makes the build
reproducible. It does not fix the heap default, the jar path, or the analytics
call; the pipeline works around the first and the second happens to be correct
inside an image because the image's `/usr` *is* the compile-time `/usr`.

## 2. Where the code can be cut

`Quandenser::run()` (`src/Quandenser.cpp:490-768`) is already commented as ten
steps. Their shapes:

| Step | What | Parallelism | Persisted between steps? |
|---|---|---|---|
| 0 | Read batch file | trivial | — |
| 1 | Dinosaur MS1 feature detection | **per input file** | **yes** — `dinosaur/<stem>.features.tsv` |
| 2 | MaRaCluster MS2 clustering, round 1 | global, OpenMP inside | yes — `.dat`, cluster TSV, precursor map |
| 3 | RT alignment models + spanning tree | global, **O(n²)**, single-threaded | **no** — in memory only |
| 4 | Match features between runs | per file pair, O(n) pairs | partly — Percolator results per pair |
| 5 | Single-linkage feature grouping | global, sequential | no |
| 6 | MaRaCluster round 2 | global, OpenMP inside | yes |
| 7 | Filter consensus features | global, sequential | no |
| 8 | Write feature groups | global, sequential | output only |
| 9 | Write consensus spectra | global, memory-bound | output only |

Two observations that matter more than the table.

**The O(n²) cost is not where you would guess.** Pairwise feature *matching*
(step 4) walks a minimum spanning tree and is O(n): exactly `2(n-1)` pairs.
The quadratic blowup is in step 3, `AlignRetention::getAlignModels`
(`src/AlignRetention.cpp:123-166`), which fits two penalised splines for every
file pair with enough shared retention-time points — single-threaded, in
memory, never checkpointed. That is what breaks a 200-file experiment.

**MaRaCluster runs twice over the whole dataset.** Step 6 repeats step 2 into a
different output folder, so none of round 1's p-value work is reused. It is the
single largest lever in the pipeline.

### Cut points, ranked by value per unit of effort

1. **Dinosaur feature detection (step 1) — zero patch.** The output path is
   deterministic and Quandenser skips its own Dinosaur call when the file
   already exists (`src/Quandenser.cpp:558`). Compute the features in separate
   tasks, stage them in, and the monolith walks straight past its slowest
   stage. This is what the pipeline does.

2. **Fix `--target-search-threshold` — one line, no pipeline needed.**
   `src/Quandenser.cpp:293-295` assigns it to `maxFeatureCandidates_` instead of
   `linkPEPMbrSearchThreshold_`, which is never settable from the CLI at all.
   The documented behaviour — "setting this to 1.0 will cause the targeted
   search to be skipped" — therefore does not work, and the second most
   expensive stage in the tool cannot be turned off. Better value than any
   amount of workflow engineering.

3. **MaRaCluster p-value fan-out (steps 2 and 6) — roughly 150 lines.** The
   `.dat` files are already format-compatible with the stock `maracluster
   pvalue`/`overlap` subcommands, and MaRaCluster ships its own Nextflow
   decomposition at `ext/maracluster/admin/nextflow/run_maracluster.nf` as a
   template. The blocker is that the spectrum-to-precursor map is serialized
   only after the whole MaRaCluster run rather than after indexing, so there is
   no clean stop point. Moving that one `serialize()` call and splitting the
   guard at `src/Quandenser.cpp:342` would open the largest fan-out in the tool,
   applied twice.

4. **Per-pair feature matching (step 4) — real refactoring.** Blocked because
   the fitted spline models (`AlignRetention::alignments_`) are never written
   to disk. Also, the "round" field is an explicit parallelism hint but only
   half correct: rounds that align children into a shared parent mutate that
   parent's feature list, so siblings conflict. The reverse-direction rounds
   are already conflict-free.

5. **Step 8 → 9 — real refactoring, little gain.** `spectrumClusterToConsensus
   Features_` passes between them purely in memory.

## 3. What this pipeline does, and what it deliberately does not

It implements cut point 1 and wraps everything else as a single process.

```
samplesheet ─┬─> DINOSAUR (one task per run, parallel) ─┐
             └─> batch.tsv ───────────────────────────► QUANDENSER ─► feature groups
                                                                   └► consensus spectra
                                                                          │
                              (optional)  CRUX_TIDE_INDEX/SEARCH ◄────────┘
                                                  └─► CRUX_PERCOLATOR ─► TRIQLER_CONVERT ─► TRIQLER
```

### The constraint that shapes everything: do not fork the binary

There was already a Nextflow port,
[quandenser-pipeline](https://github.com/statisticalbiotechnology/quandenser-pipeline).
It had a twelve-process decomposition, a GUI, Docker and Singularity images. It
is unrunnable today, and the reason is worth stating plainly: every one of its
parallel processes invoked `quandenser --partial-1-dinosaur`,
`--partial-2-maracluster`, `--partial-3-match-round`, `--partial-4-consensus`.
Those flags do not exist in this repository. They lived on a branch
(`quandenser-pipeline`) that was never merged. When that branch stopped being
maintained, the pipeline stopped working — and no user could tell why, because
the binary they installed simply rejected the flags.

That is the failure mode to avoid, and it is why the decomposition here stops
at cut point 1. Cut point 1 requires no new flags: it uses behaviour the
released binary already has. Cut point 3 is more valuable, and should be done —
but **as a merged change to `src/`, released in a versioned binary**, not as a
flag that exists only for the pipeline. A pipeline may only ever depend on a
released interface.

Two further reasons the old pipeline did not survive, both worth not repeating:
it was written in Nextflow DSL1, which no current Nextflow can execute at all;
and it wrote outside the Nextflow work directory and keyed state on a global
hash, which defeated resume and made its own bespoke resume machinery
necessary.

### Correctness details the pipeline has to handle

- **File order is file identity.** `SpectrumFileList::initFromFile` assigns each
  run its `fileIdx` by position in the batch file, *and silently skips lines
  naming files that do not exist*. Column 0 of
  `Quandenser.feature_groups.tsv` is that integer, not a name. One dropped line
  renumbers every run with no error anywhere. The pipeline therefore sorts the
  batch file deterministically and refuses to start if a listed file is not
  staged.
- **One file serves two tools.** MaRaCluster reads only the first tab-separated
  field of each batch line; Triqler reads the first two. So a single
  `path <TAB> condition` file is both the Quandenser batch file and the Triqler
  file list, and there is no second file to drift out of sync.
- **Base names must be unique.** Quandenser keys intermediates on the file stem
  and has a standing TODO about collisions (`src/Quandenser.cpp:555`). The
  pipeline checks this up front.
- **`OMP_NUM_THREADS` must be set explicitly.** MaRaCluster has no thread
  option of its own and `--num-threads` never reaches it; without the
  environment variable it takes every core on the node regardless of the
  scheduler's allocation.
- **`--max-missing` defaults to `numFiles/4`**, so leaving it unset makes a
  result depend on how many runs happened to be in the batch.
- **Percolator throws above 128 threads** even though Quandenser's validation
  accepts up to 1000.

### Known-broken options, not exposed as parameters

- `-c/--ft-link-candidates` assigns to `linkPEPThreshold_` instead of
  `maxFeatureCandidates_` (`src/Quandenser.cpp:289-291`), so setting it to 2
  sets the PEP cut-off to 2.0 and disables it.
- `-B/--target-search-threshold` assigns to `maxFeatureCandidates_`
  (`:293-295`).
- `-t/--percolator-test-fdr` writes into the *train* FDR option
  (`:261-267`), so setting only `-t` passes an empty string to Percolator. The
  pipeline always passes both FDRs together, which avoids this.

These are upstream bugs, not pipeline limitations. Exposing them as parameters
would surface flags that do the wrong thing.

### The downstream half

`bin/prepare_input.py` is not used. It unpacks `parsers.parseFileList()` into
three values, but current Triqler returns a single list, so it raises
`ValueError` on any modern install. It and `bin/normalize_intensities.py` and
`bin/percolator.py` are stale forks of code that has since moved into the
`triqler` package. The pipeline calls `python -m triqler.convert.quandenser`
instead. The `bin/` scripts are also never installed by CMake, so a user
following the ReadMe would not have them anyway.

## 4. Suggested next steps for the code itself

In the order I would do them:

1. Pin ProteoWizard to a fixed version and make the Boost header copy resolve
   the directory rather than hard-code it. Without this nothing else matters,
   because nothing builds. `containers/quandenser/Dockerfile` contains a
   working repair that could be moved upstream into MaRaCluster's script.
2. Fix the three option-parsing bugs listed above. All are one-liners.
3. Lower the Dinosaur heap default from 24 GB to something a laptop has, and
   let an environment variable override the jar path.
4. Repair or replace the release workflow, and have it publish a container
   image. A Dockerfile that CI never builds will rot exactly as the current
   workflow has.
5. Then, and only then, add `--stop-after`/`--start-from` for the MaRaCluster
   fan-out — merged to master and released, so the pipeline can depend on it.
