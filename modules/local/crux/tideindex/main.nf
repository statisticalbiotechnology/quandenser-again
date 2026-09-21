process CRUX_TIDE_INDEX {
    tag "${fasta.baseName}"
    label 'process_medium'

    container 'quay.io/biocontainers/crux-toolkit:4.1--h503566f_3'

    input:
    path fasta

    output:
    path "tide_index" , emit: index
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    crux tide-index \\
        --overwrite T \\
        --decoy-format protein-reverse \\
        --missed-cleavages ${params.allowed_missed_cleavages} \\
        --enzyme ${params.enzyme} \\
        ${args} \\
        ${fasta} \\
        tide_index

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        crux: \$(crux version 2>&1 | grep -oP 'Crux version \\K[0-9.]+' || echo 4.2)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p tide_index
    touch tide_index/.stub
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        crux: 4.2
    END_VERSIONS
    """
}
