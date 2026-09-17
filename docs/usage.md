# Usage

## Quick start

```bash
nextflow run statisticalbiotechnology/quandenser-again \
    -profile docker \
    --input samplesheet.csv \
    --outdir results
```

Requirements: Nextflow ≥ 23.10 and one of Docker, Singularity, Apptainer or
Podman. Nothing else. In particular you do not need a compiler, Java, Maven,
Boost or ProteoWizard on the host.

To check the plumbing before committing real data to it:

```bash
nextflow run statisticalbiotechnology/quandenser-again -profile test,docker --outdir results
```

That downloads three small public mzML files and runs the real tools on them,
through to feature groups and consensus spectra. Two settings in
`conf/test.config` accommodate the size of that dataset rather than
representing advice for a real analysis; both are explained where they are
set.

## Samplesheet

A CSV with a header. `sample` and `mzml` are required; `condition` is needed
only if you intend to run Triqler.

```csv
sample,mzml,condition
control_1,/data/control_1.mzML,control
control_2,/data/control_2.mzML,control
treated_1,/data/treated_1.mzML,treated
treated_2,/data/treated_2.mzML,treated
```

| Column | Required | Meaning |
|---|---|---|
| `sample` | yes | Run name. Must be unique. |
| `mzml` | yes | Path or URL to a peak-picked mzML file. |
| `condition` | for Triqler | Experimental group. Defaults to `A`. |

Two constraints that come from Quandenser itself:

- **At least two runs.** Quandenser aligns runs against each other and refuses
  to work with fewer than two.
- **Base names must be unique.** `/a/run.mzML` and `/b/run.mzML` in the same
  samplesheet will not work, because Quandenser keys its intermediate files on
  the file stem. The pipeline stops with an explicit error rather than letting
  one silently overwrite the other.

The mzML files must be peak-picked (centroided) on **both** MS1 and MS2 level.
Quandenser does not do this for you, and neither does this pipeline. With
ProteoWizard's `msconvert` that is `--filter "peakPicking true 1-"`.

### A note on retention-time units

Either unit is valid in mzML and the file records which it means in
`unitName`. Quandenser reads the number and ignores the unit, while Dinosaur
always reports minutes, so a file declaring seconds matches no features at
all: every MS2 scan is dropped for having no precursor and the run fails with
`could not find any ms2 spectra in the input files`, which is not what went
wrong.

The pipeline detects the unit and converts the feature file to match, so
either kind of file works. This only applies on the default path, where
feature detection is a separate step; with `--parallel_feature_detection
false` Quandenser runs Dinosaur itself and the problem is untouched. See
`docs/modularization.md`.

## What you get

```
results/
├── features/                          one <run>.features.tsv per input
├── quandenser/
│   ├── Quandenser.feature_groups.tsv  the quantification matrix
│   ├── consensus_spectra/             .ms2 files to search
│   ├── maracluster/                   MS2 clusters
│   └── dinosaur/                      features as Quandenser saw them
├── triqler/                           only with --run_identification
└── pipeline_info/
    ├── batch.tsv                      run order, i.e. the fileIdx mapping
    └── ...                            timeline, report, trace, DAG, versions
```

`Quandenser.feature_groups.tsv` has **no header**. It is a sequence of
blank-line-separated blocks, one per feature group, each line being

```
fileIdx  precMz  charge  rTime  intensity  spectrumId;linkPEP[,spectrumId;linkPEP...]
```

`fileIdx` is an integer index into the batch file, not a file name. The
pipeline sorts the batch file by run name, so index 0 is the alphabetically
first mzML, index 1 the second, and so on. The batch file itself is published
as `pipeline_info/batch.tsv`, which is the only record of that mapping — keep
it with the results.

## Identification and Triqler

Off by default, because it needs a protein database:

```bash
nextflow run statisticalbiotechnology/quandenser-again \
    -profile docker \
    --input samplesheet.csv \
    --fasta uniprot.fasta \
    --run_identification \
    --outdir results
```

This searches the consensus spectra with Crux Tide, rescores with Percolator,
maps the identifications back onto feature groups with Triqler's own converter,
and runs Triqler. Any other search engine works equally well — Quandenser is
agnostic about it — but you would need to add a module and keep the scan
numbering intact (see `docs/modularization.md`).

## Parameters

### Input/output

| Parameter | Default | Description |
|---|---|---|
| `--input` | — | Samplesheet CSV. Required. |
| `--outdir` | `results` | Output directory. |

### Quandenser

| Parameter | Default | Description |
|---|---|---|
| `--parallel_feature_detection` | `true` | Run Dinosaur as one task per file. Setting this false reverts to Quandenser's internal serial loop. |
| `--prefix` | `Quandenser` | Output file prefix. |
| `--seed` | `1` | Random seed for Dinosaur and Percolator. |
| `--max_missing` | `null` | Maximum missing values per quantification row. Unset means Quandenser's own default of `numFiles/4`, which makes the result depend on the batch size; set it for reproducibility across experiments. |
| `--maracluster_pvalue_cutoff` | `-10.0` | log(p) threshold for MS2 clustering. |
| `--align_mz_tol` | `20.0` | Mass tolerance for matching features between runs, ppm. |
| `--align_rtime_tol` | `10.0` | Retention-time tolerance, in standard deviations of the alignment. |
| `--intensity_score_cut` | `0.5` | Intensity score ratio cut-off per cluster. |
| `--ft_link_cut` | `0.25` | PEP cut-off for matching features between runs. |
| `--percolator_train_fdr` / `--percolator_test_fdr` | `0.02` | FDR thresholds for the internal Percolator. Always passed together, because passing only one triggers an upstream option-parsing bug. |
| `--ft_link_candidates` | `null` | Candidates considered when matching features between runs. |
| `--target_search_threshold` | `null` | Minimum PEP for re-searching a feature with a targeted search. `1.0` skips the targeted search, the second most expensive stage in the tool. |
| `--quandenser_args` | `null` | Anything else, appended verbatim. |

The last two need a Quandenser built from this repository. In released 0.03.2
they assign to each other's variables (`src/Quandenser.cpp:289-295`), so
`--target-search-threshold` silently sets a candidate count and
`--ft-link-candidates` silently overwrites the feature-link PEP cut-off. Both
are accepted by an older binary rather than rejected, so the pipeline still
runs against one; it just will not skip the targeted search. See
`docs/modularization.md`.

### Identification and Triqler

| Parameter | Default | Description |
|---|---|---|
| `--run_identification` | `false` | Search the consensus spectra. |
| `--fasta` | `null` | Protein database. Required with the above. |
| `--enzyme` | `trypsin/p` | Digestion enzyme. |
| `--allowed_missed_cleavages` | `2` | |
| `--precursor_mass_tolerance` | `20.0` | With `--precursor_mass_tolerance_unit`, default `ppm`. |
| `--skip_triqler` | `false` | Stop after producing the Triqler input file. |
| `--fold_change_eval` | `0.8` | Triqler fold-change evaluation threshold, log2. |
| `--retain_unidentified` | `true` | Keep feature groups with no identification. |
| `--skip_normalization` | `false` | Skip retention-time normalisation. Necessary on small datasets: the normalisation averages over a fixed 2000-point window and fails on inputs yielding fewer feature groups than that. |

### Resources

| Parameter | Default | Description |
|---|---|---|
| `--max_cpus` | `16` | |
| `--max_memory` | `128.GB` | |
| `--max_time` | `240.h` | |

On a laptop, lower them:

```bash
nextflow run . -profile docker --input s.csv --max_cpus 4 --max_memory 12.GB
```

Requests are clamped to these ceilings, so a task asking for 48 GB on a 12 GB
machine gets 12 GB rather than never being scheduled.

Memory matters more than usual here. Dinosaur defaults to asking the JVM for a
24 GB heap, and no environment variable can change it — the `-Xmx` is compiled
into the command line Quandenser builds. The pipeline derives the heap from
what each task was actually granted, which is the only thing that works. If you
run the container by hand, pass `--dinosaur-memory` yourself.

## Resuming

```bash
nextflow run ... -resume
```

Note that Quandenser has its own informal resume mechanism: it skips work whose
output file already exists. That is convenient by hand and wrong under a
workflow engine, because it is keyed on existence rather than on content or
parameters — change a parameter and it will happily reuse intermediates
computed with the old one. Every task here gets a clean directory so that
Nextflow's content-addressed caching is the only caching in play. Do not point
the pipeline at a pre-existing Quandenser output folder.

## Containers

There is no Quandenser package on bioconda or biocontainers, so the image is
built from this repository:

```bash
git submodule update --init --recursive
docker build -t quandenser:0.03.2 -f containers/quandenser/Dockerfile .
nextflow run . -profile docker --input s.csv --quandenser_container quandenser:0.03.2
```

Expect the build to take the better part of an hour, most of it compiling Boost
and ProteoWizard. The other tools (Crux, Triqler) come from existing
biocontainers.
