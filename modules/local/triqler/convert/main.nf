/*
 * Join identifications onto Quandenser feature groups and emit Triqler input.
 *
 * This deliberately calls triqler's own converter rather than the repository's
 * bin/prepare_input.py. That script unpacks parsers.parseFileList() into three
 * values, but current Triqler returns a single list, so it raises ValueError
 * on any modern install. It is a fork of code that has since moved into the
 * triqler package; the packaged version is the maintained one.
 */
process TRIQLER_CONVERT {
    tag "triqler_input"
    label 'process_medium'

    container 'quay.io/biocontainers/triqler:0.9.2--pyhdfd78af_0'

    input:
    path feature_groups
    path batch_file
    path target_psms
    path decoy_psms

    output:
    path "triqler_input.tsv", emit: triqler_input
    path "versions.yml"     , emit: versions

    script:
    def args = task.ext.args ?: ''
    def retain = params.retain_unidentified ? '--retain_unidentified' : ''
    // Quandenser's own RT normalisation lives in triqler; skipping it is the
    // safe choice on small inputs, where the hard-coded 2000-point running
    // mean window has nothing to average over.
    def skip_norm = params.skip_normalization ? '--skip_normalization' : ''
    """
    python -m triqler.convert.quandenser \\
        --file_list_file ${batch_file} \\
        --psm_files ${target_psms},${decoy_psms} \\
        --out_file triqler_input.tsv \\
        ${retain} \\
        ${skip_norm} \\
        ${args} \\
        ${feature_groups}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        triqler: \$(python -c "import triqler; print(triqler.__version__)" 2>/dev/null || echo 0.9.2)
    END_VERSIONS
    """

    stub:
    """
    touch triqler_input.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        triqler: 0.9.2
    END_VERSIONS
    """
}
