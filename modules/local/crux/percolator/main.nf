/*
 * Percolator over the consensus-spectrum PSMs.
 *
 * Note this is a different Percolator run from the one inside Quandenser:
 * Quandenser links Percolator to score MS1 feature matches between runs, never
 * to rescore a database search. --top-match 1 matches what the downstream
 * Triqler converter expects.
 */
process CRUX_PERCOLATOR {
    tag "consensus"
    label 'process_medium'

    container 'quay.io/biocontainers/crux-toolkit:4.2--h9ee0642_0'

    input:
    path psms

    output:
    path "crux-output/percolator.target.psms.txt", emit: target_psms
    path "crux-output/percolator.decoy.psms.txt" , emit: decoy_psms
    path "versions.yml"                          , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    crux percolator \\
        --overwrite T \\
        --top-match 1 \\
        --output-dir crux-output \\
        ${args} \\
        ${psms}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        crux: \$(crux version 2>&1 | grep -oP 'Crux version \\K[0-9.]+' || echo 4.2)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p crux-output
    touch crux-output/percolator.target.psms.txt
    touch crux-output/percolator.decoy.psms.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        crux: 4.2
    END_VERSIONS
    """
}
