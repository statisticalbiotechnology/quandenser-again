/*
 * Quandenser: MS2 clustering, retention-time alignment, match-between-runs,
 * feature grouping and consensus-spectrum generation.
 *
 * This runs the released binary unmodified. That is a deliberate constraint,
 * not a limitation we failed to overcome: the previous Nextflow port of this
 * tool (statisticalbiotechnology/quandenser-pipeline) drove Quandenser through
 * --partial-1-dinosaur, --partial-2-maracluster and friends, which exist only
 * on a branch that was never merged. When that branch was abandoned the
 * pipeline became unrunnable, and it is unrunnable today. Every flag used
 * here exists in src/Quandenser.cpp on master.
 *
 * Feature detection is normally done by the DINOSAUR module and staged in
 * below, so the loop at src/Quandenser.cpp:547-572 finds its outputs already
 * present and skips straight to clustering.
 *
 * Note what is *not* passed: --verbatim. Quandenser forwards it to both
 * MaRaCluster and Percolator (src/Quandenser.cpp:216-220), but Percolator has
 * no such option, so setting it aborts the run at the feature-matching stage
 * with "ERROR: the option --verbatim is invalid" — after Dinosaur and two
 * MaRaCluster passes have already done their work. The default verbosity is
 * what we want anyway.
 */

process QUANDENSER {
    tag "${params.prefix}"
    label 'process_high'

    container "${params.quandenser_container}"

    input:
    path mzmls
    path batch_file
    path dinosaur_features, stageAs: 'precomputed_features/*'

    output:
    path "${params.prefix}.feature_groups.tsv", emit: feature_groups
    path "consensus_spectra/*"                , emit: consensus_spectra
    path "maracluster/*.clusters_p*.tsv"      , emit: clusters, optional: true
    path "dinosaur/*.features.tsv"            , emit: features, optional: true
    path "versions.yml"                       , emit: versions

    script:
    def args     = task.ext.args ?: ''
    def prefix   = params.prefix
    // Percolator throws above 128 threads (ext/percolator/src/Caller.cpp:632)
    // even though Quandenser's own validation accepts up to 1000.
    def nthreads = Math.min(task.cpus as int, 128)
    def heap_mb  = Math.max(1024L, (task.memory.toMega() * 0.5) as long)
    // Quandenser defaults this to numFiles/4, so leaving it unset makes the
    // result depend on how many runs are in the batch. Pin it when given.
    def max_missing = params.max_missing != null ? "--max-missing ${params.max_missing}" : ''
    """
    # MaRaCluster is parallelised purely with OpenMP and has no thread option
    # of its own; --num-threads never reaches it. Without this it takes every
    # core on the node regardless of what the scheduler granted.
    export OMP_NUM_THREADS=${nthreads}

    # The batch file must name files that exist. SpectrumFileList::initFromFile
    # silently drops missing paths and then assigns each surviving run its
    # fileIdx by position, so one dropped line renumbers every run and every
    # downstream artefact refers to the wrong file, with no error anywhere.
    while IFS=\$'\\t' read -r f _rest; do
        [ -z "\$f" ] && continue
        if [ ! -e "\$f" ]; then
            echo "ERROR: batch file names a file that is not staged: \$f" >&2
            exit 1
        fi
    done < ${batch_file}

    # Hand Quandenser the features the DINOSAUR tasks already computed.
    mkdir -p dinosaur
    if [ -d precomputed_features ]; then
        for f in precomputed_features/*.features.tsv; do
            [ -e "\$f" ] || continue
            cp -L "\$f" dinosaur/
        done
    fi

    quandenser \\
        --batch ${batch_file} \\
        --output-folder . \\
        --prefix ${prefix} \\
        --seed ${params.seed} \\
        --num-threads ${nthreads} \\
        --dinosaur-threads ${nthreads} \\
        --dinosaur-memory ${heap_mb}M \\
        --maracluster-pval-threshold ${params.maracluster_pvalue_cutoff} \\
        --align-mz-tol ${params.align_mz_tol} \\
        --align-rtime-tol ${params.align_rtime_tol} \\
        --intensity-score-cut ${params.intensity_score_cut} \\
        --ft-link-cut ${params.ft_link_cut} \\
        --percolator-train-fdr ${params.percolator_train_fdr} \\
        --percolator-test-fdr ${params.percolator_test_fdr} \\
        ${max_missing} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quandenser: ${params.quandenser_version}
    END_VERSIONS
    """

    stub:
    def prefix = params.prefix
    """
    mkdir -p consensus_spectra maracluster dinosaur
    touch ${prefix}.feature_groups.tsv
    touch consensus_spectra/${prefix}.consensus.part1.ms2
    touch maracluster/${prefix}.clusters_p10.tsv
    for f in precomputed_features/*.features.tsv; do
        [ -e "\$f" ] || continue
        cp -L "\$f" dinosaur/
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quandenser: ${params.quandenser_version}
    END_VERSIONS
    """
}
