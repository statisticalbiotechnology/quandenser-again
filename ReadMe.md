# Quandenser: QUANtification by Distillation for ENhanced Signals with Error Regulation

Quandenser condenses quantification data from label-free mass spectrometry experiments.

## Running it as a pipeline (recommended)

The simplest way to use Quandenser is the Nextflow pipeline in this
repository. It needs Nextflow and a container runtime, and nothing else — no
compiler, no Java, no Boost, no ProteoWizard:

```bash
nextflow run statisticalbiotechnology/quandenser-again \
    -profile docker \
    --input samplesheet.csv \
    --outdir results
```

To check it works before using your own data, `-profile test,docker` runs the
whole thing on three small public files.

The pipeline also runs MS1 feature detection as one task per input file rather
than in a serial loop, which is where most of the wall time goes, and it
optionally continues through database search and Triqler to protein-level
differential abundance.

See [docs/usage.md](docs/usage.md) for the samplesheet format and parameters,
and [docs/modularization.md](docs/modularization.md) for how the pipeline is
put together and why.

## Installation

An installer for all major platforms (Windows, OS X, Ubuntu, etc.) can be found on the Release page. Java 8 or later has to be installed.

If you prefer to compile from source, or are running on a different operating system, [click here](#installation-from-source).

> **Note on building from source.** The build currently fails on a clean
> machine. `ext/maracluster/admin/builders/install_proteowizard.sh` downloads
> whatever ProteoWizard built most recently, but copies Boost headers from a
> hard-coded `boost_1_76_0` directory that recent ProteoWizard releases no
> longer contain. The copy fails, the script does not check for it, and the
> compile dies later with a missing-header error that does not point at the
> cause. `containers/quandenser/Dockerfile` works around it; see
> [docs/modularization.md](docs/modularization.md).

An older Docker and Singularity wrapper exists in the
[Quandenser-pipeline project](https://github.com/statisticalbiotechnology/quandenser-pipeline),
but it is written in Nextflow DSL1, which current Nextflow releases cannot run,
and it depends on Quandenser flags that were never merged into this
repository. Prefer the pipeline above.

## Interface

Usage:

```
  quandenser -b <batch_file_list> -f <output_directory>
```

`<batch_file_list>` is a flat text file with the absolute path to each of the mzML files (one per line). The mzML files should be peak picked on both MS1 and MS2 level.

If everything executed correctly, you will find a file with the quantified feature groups at `<output_directory>/Quandenser.feature_groups.tsv` and one or several files with consensus spectra in `<output_directory>/consensus_spectra/MaRaCluster.consensus.part<x>.ms2`. Different output formats for the consensus spectra files can be specified using the `-o` option and the consensus spectra can subsequently be searched by any search engine.

## Example

An example run including downstream analysis with [Triqler](https://github.com/statisticalbiotechnology/triqler) can be found here: https://app.box.com/s/kp4219dc22l3gq27014nms8oco594c2i

The folder contains a ReadMe file with instructions on how to run the example.

Note that this same pipeline is also implemented in the [Quandenser-pipeline project](https://github.com/statisticalbiotechnology/quandenser-pipeline), which allows users to run the entire pipeline from RAW files to Triqler results using [Nextflow](https://www.nextflow.io/) with pre-built [Singularity](https://sylabs.io/docs/) or [Docker](https://www.docker.com/) images.

## Installation from source

First, clone the repository using `git clone --recursive https://github.com/statisticalbiotechnology/quandenser.git`. Note the `--recursive` flag which is needed to pull in the submodules.

Quandenser depends on the Proteowizard and Boost libraries, which are automatically installed by the builder scripts. 

To install Quandenser, you can use the provided installation script `./quickbuild.sh` (Unix) or `./quickbuild.bat`/`./quickbuild64.bat` (Windows 32-bit and 64-bit respectively), which calls the appropriate build script for your platform located at `admin/builders/<platform>_build.<ext>`. By default, it will install the executables in the `/usr/bin` folder (needs superuser rights). If you do not have superuser rights or want to install the executable somewhere else, modify the script accordingly by setting the `-DCMAKE_INSTALL_PREFIX` flag to the desired location.
