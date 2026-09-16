process TRIQLER {
    tag "fc_${params.fold_change_eval}"
    label 'process_high'

    container 'quay.io/biocontainers/triqler:0.9.2--pyhdfd78af_0'

    input:
    path triqler_input

    output:
    path "proteins.*.tsv", emit: proteins
    path "versions.yml"  , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    python -m triqler \\
        --fold_change_eval ${params.fold_change_eval} \\
        --num_threads ${task.cpus} \\
        --out_file proteins.tsv \\
        ${args} \\
        ${triqler_input}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        triqler: \$(python -c "import triqler; print(triqler.__version__)" 2>/dev/null || echo 0.9.2)
    END_VERSIONS
    """

    stub:
    """
    touch proteins.1vs2.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        triqler: 0.9.2
    END_VERSIONS
    """
}
