/*
 * MS1 feature detection with biosaur2, one task per run: the alternative to
 * DINOSAUR that also handles timsTOF (PASEF) data.
 *
 * The output takes Dinosaur's place exactly: dinosaur/<stem>.features.tsv,
 * named after the mzML that Quandenser itself reads, so the monolith finds it
 * and skips its own Dinosaur call (src/Quandenser.cpp:558).
 *
 * Two things differ from the Dinosaur path.
 *
 * Column order. Quandenser parses feature files by position in Dinosaur's
 * layout, and biosaur2 writes the same column names in another order, which
 * would be read without error and with every value misplaced. The output is
 * therefore rewritten by bin/biosaur2_to_dinosaur.py.
 *
 * Retention time. biosaur2 takes scan start time from the mzML as it stands,
 * and so does MaRaCluster, so the minutes-to-seconds conversion the Dinosaur
 * module needs does not apply here.
 *
 * For timsTOF, biosaur2 wants an MS1-only mzML with ion-mobility spectra
 * combined and *not* scan-summed, while MaRaCluster wants the PASEF MS2
 * spectra. The samplesheet's optional ms1_mzml column names the first; the
 * mzml column names the second, and is used for both when ms1_mzml is absent,
 * which is the Orbitrap case.
 *
 * biosaur2 runs single-process: its exact fast path (markmipt/biosaur2#22)
 * applies only with -nprocs 1, and is about fifty times faster on PASEF runs
 * than the multi-process original.
 */

process BIOSAUR2 {
    tag "$meta.id"

    container "${params.biosaur2_container}"

    input:
    tuple val(meta), path(mzml), path(ms1_mzml, stageAs: 'ms1/*')

    output:
    tuple val(meta), path("dinosaur/${mzml.baseName}.features.tsv"), emit: features
    path "versions.yml"                                            , emit: versions

    script:
    def args = task.ext.args ?: ''
    // Long runs in retention-time chunks, which bounds memory; see the script.
    def run = params.biosaur2_chunk_minutes ?
        "biosaur2_chunked.py ${ms1_mzml} biosaur2.features.tsv ${params.biosaur2_chunk_minutes}" :
        "biosaur2 ${ms1_mzml} -o biosaur2.features.tsv"
    def min_scans = params.biosaur2_min_scans ?: ''
    """
    mkdir -p dinosaur

    ${run} -nprocs 1 ${args}

    # Quandenser tests for existence, not for content (src/Quandenser.cpp:558).
    if [ ! -s biosaur2.features.tsv ] || [ \$(wc -l < biosaur2.features.tsv) -lt 2 ]; then
        echo "biosaur2 produced no features for ${ms1_mzml}" >&2
        exit 1
    fi

    biosaur2_to_dinosaur.py biosaur2.features.tsv ${min_scans} > dinosaur/${mzml.baseName}.features.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biosaur2: \$(pip show biosaur2 2>/dev/null | sed -n 's/^Version: //p')
        biosaur2_fast_path: \$(python -c "import biosaur2.fastpath as f; print(f.FAST_OK)")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p dinosaur
    printf 'mz\\tmostAbundantMz\\tcharge\\trtStart\\trtApex\\trtEnd\\tfwhm\\tnIsotopes\\tnScans\\taveragineCorr\\tmass\\tmassCalib\\tintensityApex\\tintensitySum\\n' \\
        > dinosaur/${mzml.baseName}.features.tsv
    printf '500.25\\t500.25\\t2\\t10.0\\t10.5\\t11.0\\t0\\t3\\t12\\t0\\t998.5\\t998.5\\t1000\\t5000\\n' \\
        >> dinosaur/${mzml.baseName}.features.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biosaur2: stub
    END_VERSIONS
    """
}
