# DIA feasibility analysis

Scripts behind `docs/dia-feasibility.md`. Read that first; this file only says
how to rerun things.

These are research scripts, not pipeline code. They hard-code the data
directories given below and are meant to be run by hand in order. They are kept
because the study's conclusions are largely negative and the numbers should be
checkable rather than taken on trust.

`RESULTS-raw.md` is the working log written as the numbers came in, including
the attempts that were wrong. `docs/dia-feasibility.md` is the version to read.

## Requirements

    pip install numpy scipy pyteomics psims lxml scikit-learn

## Data layout

The scripts expect two directories. Nothing here is committed; all of it is
public and re-downloadable.

### `/home/user/csf` (PXD004863, DDA)

    raw/    20150513_20_hTau_CSF_Frac9of12_rep1.raw
            20150513_48_hTau_CSF_Frac9of12_rep2.raw
    mzml/   the same two, converted
    calib/  20150513_20_hTau_CSF_Frac9of12_rep1.mgf
            20150513_48_hTau_CSF_Frac9of12_rep2.mgf
            peaks.mzid
    results/  output of the Nextflow pipeline on the two mzML files

Base URL for all four raw/mgf/mzid files:

    https://ftp.pride.ebi.ac.uk/pride/data/archive/2017/01/PXD004863/

The identifications are `PEAKS_SwissProt_Human_Data_peptides_1_1_0.mzid`, saved
as `calib/peaks.mzid`. Conversion of the RAW files:

    docker run --rm -v "$PWD":/data \
        quay.io/biocontainers/thermorawfileparser:1.4.5--ha8f3691_0 \
        ThermoRawFileParser.sh -i /data/raw/<f>.raw -b /data/mzml/<f>.mzML -f 2

Peak picking is on by default; do not pass `-p`, which disables it.

### `/home/user/dia` (PXD001587, DIA)

    18484_REP3_1ug_Ecoli_NewStock2_SWATH_1.mzXML
    18486_REP3_1ug_Ecoli_NewStock2_SWATH_2.mzXML
    Ecoli_DDA_CombinedLib.sptxt
    split_18484_REP3_1ug_Ecoli_NewStock2_SWATH_1_mprophet_all_peakgroups_DDA.xls
    split_18486_REP3_1ug_Ecoli_NewStock2_SWATH_2_mprophet_all_peakgroups_DDA.xls

Base URL:

    https://ftp.pride.ebi.ac.uk/pride/data/archive/2015/01/PXD001587/

The two mzXML files download directly. The library and the two mProphet tables
live inside `Ecoli.zip` (1.1 GB) at

    OpenSWTH/Ecoli_DDA_CombinedLib.sptxt
    OpenSWTH/lib-Ecoli_DDA_CombinedLib_mins1_trans6_decoymode-normal/<run>/split_<run>_mprophet_all_peakgroups_DDA.xls

Downloading the whole zip for 50 MB of content is avoidable: PRIDE's server
honours range requests, so the central directory can be read from the tail and
individual members fetched and inflated. About 30 seconds instead of 20
minutes.

## Order of execution

### DDA, channel 2

    comp_pairs.py          consensus spectra. Superseded: MaRaCluster bins
                           peaks, so consensus spectra cannot support a 0.01 Da
                           test. Kept because that is worth knowing.
    range.py               shows why: 60% of charge-4 consensus spectra have no
                           peak above half the precursor-pair sum.
    comp_pairs_mzml.py     all MS2 from the mzML, both pair relations.
    compare_mgf_mzml.py    control. Confirms the deposited MGF is a denoised
                           subset of the mzML peaks and not deisotoped.
    identified.py          the real calibration, on PEAKS-identified spectra.
    chimera.py             multiplexing simulation. The decisive one.
    tolscan.py             mass-accuracy dependence.

### DDA, channel 1 proxy

    xic.py <run>           MS1 XICs for mid-gradient features.
    corr3.py               discrimination by profile correlation, stratified.
    crossrun2.py           cross-run pooling, with the retention-time warp
                           estimated rather than assumed.

### DIA

    scores.py              OpenSWATH sub-scores against its own decoys.
    combine.py             are co-elution and identity independent. No.
    diatest.py <run>       library fragment traces and merged window spectra.
    analyse.py             both channels at the fragment level.
    ch2_sweep.py           channel 2 with denoising.
    indep.py               independence with a channel-2 score that has signal.
    crossdia.py            cross-run pooling on DIA.

### DIA, quantification-first linking

    pseudospec.py <run>    per-feature candidate traces (superseded by the
                           block form below, kept for linktest.py).
    linktest.py            are library-free signatures reproducible and
                           specific.
    linkstrat.py           the same as a retrieval task, by intensity.
    ms1feat.py <run>       MS1 feature detection, about 97k features per run.
    blockspec.py <run>     candidate traces per (window, 60 s block).
    linkfull.py            retrieval against the full detected field. Also
                           defines the helpers the next four import by exec.
    linkreal.py            retrieval as a pipeline actually poses it, with
                           candidates constrained by m/z.
    linkvalue.py           alignment-residual crossover, and link validation.
    yield.py               how many unidentified features link confidently.
    redundancy.py          collapses isotope and split-trace duplicates.

`linkreal.py`, `linkvalue.py`, `yield.py` and `redundancy.py` read helper
definitions out of `linkfull.py` and `yield.py` with `exec`, so those files
must stay next to them.

## Traps met along the way

Recorded because each cost a debugging cycle and each would recur.

- pyteomics reports mzXML `retentionTime` in minutes. OpenSWATH's `RT`,
  `leftWidth` and `rightWidth` are in seconds. Mixing them matches nothing at
  all, silently, and the scripts still run to completion.
- Decoy sets contaminate easily. Two features co-eluting in one isolation
  window are very often the same peptide at two charge states, or one feature
  split in two. Counting them as interference destroys the measurement.
- `hash()` on strings is randomised per process. Cross-run work that picks a
  partner by hash needs `zlib.crc32` or the two runs disagree.
- Decoy masses must be shifted by whole Da. The pair-sum distribution is
  periodic at 1 Da and a fractional shift compares against a different phase.
- A correlation window that is mostly baseline inflates the correlation of any
  two traces, true or decoy, because they agree on being zero.
