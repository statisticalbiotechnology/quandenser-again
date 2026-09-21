/*
 * Search the consensus spectra.
 *
 * The scan numbers in these .ms2 files are not arbitrary: Quandenser writes
 * spectrumId = spectrumClusterIdx * 100 + offset (src/FeatureGroups.cpp:318),
 * and that integer is what reappears in the Percolator PSMId and is how
 * identifications are later joined back onto feature groups. Nothing here may
 * renumber or re-sort spectra.
 */
process CRUX_TIDE_SEARCH {
    tag "consensus"
    label 'process_high'

    container 'quay.io/biocontainers/crux-toolkit:4.1--h503566f_3'

    input:
    path consensus_spectra
    path index

    output:
    path "crux-output/tide-search.txt", emit: psms
    path "versions.yml"               , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    crux tide-search \\
        --overwrite T \\
        --concat T \\
        --num-threads ${task.cpus} \\
        --precursor-window ${params.precursor_mass_tolerance} \\
        --precursor-window-type ${params.precursor_mass_tolerance_unit} \\
        --output-dir crux-output \\
        ${args} \\
        ${consensus_spectra} \\
        ${index}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        crux: \$(crux version 2>&1 | grep -oP 'Crux version \\K[0-9.]+' || echo 4.2)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p crux-output
    touch crux-output/tide-search.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        crux: 4.2
    END_VERSIONS
    """
}
