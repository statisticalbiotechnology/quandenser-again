/*
 * MS1 feature detection with Dinosaur, one task per run.
 *
 * This is the whole reason the pipeline is worth having. Quandenser detects
 * features by shelling out to a JVM once per input file, in a plain serial
 * loop (src/Quandenser.cpp:547-572), so on a 20-run experiment the slowest
 * stage runs 20 times back to back on one core's worth of useful work.
 *
 * Splitting it out needs no change to Quandenser. The binary skips its own
 * Dinosaur call whenever the expected output file already exists
 * (src/Quandenser.cpp:558), so if we compute the features here and stage them
 * into <output-folder>/dinosaur/ the monolith picks them up and moves on.
 *
 * The command below reproduces src/DinosaurIO.cpp:51-56 and :44 flag for
 * flag. Diverging from it would silently produce different features than a
 * monolithic run, so treat those two source lines as the contract.
 */

process DINOSAUR {
    tag "$meta.id"
    label 'process_medium'

    container "${params.quandenser_container}"

    input:
    tuple val(meta), path(mzml)

    output:
    tuple val(meta), path("dinosaur/${mzml.baseName}.features.tsv"), emit: features
    path "versions.yml"                                            , emit: versions

    script:
    def args     = task.ext.args ?: ''
    def jar_path = params.quandenser_jar_path
    // Dinosaur's own default is a 24 GB heap (src/DinosaurIO.cpp:26), which is
    // fatal under any scheduler that enforces the memory request. Derive the
    // heap from what the task was actually granted and leave the JVM room for
    // its non-heap overhead.
    def heap_mb  = Math.max(1024L, (task.memory.toMega() * 0.75) as long)
    """
    mkdir -p dinosaur

    java -Xmx${heap_mb}M -jar ${jar_path}/Dinosaur-1.2.1.free.jar \\
        --force \\
        --profiling=true \\
        --nReport=0 \\
        --concurrency=${task.cpus} \\
        --seed=${params.seed} \\
        --outDir=dinosaur \\
        --advParams=${jar_path}/advParams_dinosaur.txt \\
        ${args} \\
        ${mzml}

    # Quandenser tests for existence, not for content (src/Quandenser.cpp:558),
    # so a truncated file here would be silently accepted on the far side.
    if [ ! -s "dinosaur/${mzml.baseName}.features.tsv" ]; then
        echo "Dinosaur produced no features for ${mzml}" >&2
        exit 1
    fi

    # Put the feature retention times in the same unit as the mzML.
    #
    # Dinosaur always reports retention time in minutes. Quandenser compares
    # those values directly against the mzML's scan start time
    # (src/SpectrumFiles.cpp:38, it->rtStart <= rTime && it->rtEnd >= rTime),
    # and MaRaCluster's SpectrumHandler::getRetentionTime returns that cvParam
    # verbatim, with no unit conversion. An mzML that declares seconds, which
    # is valid and common, therefore matches no features at all: every MS2 scan
    # is discarded for having no precursor, and the run dies with "could not
    # find any ms2 spectra in the input files" — a message that blames the
    # spectra rather than the units. Converting here is the narrowest fix, and
    # it is only possible because feature detection is a separate step.
    rt_unit=\$(grep -m1 -o 'accession="MS:1000016"[^>]*' "${mzml}" \\
               | grep -o 'unitName="[a-zA-Z]*"' | cut -d'"' -f2 || true)
    echo "mzML scan start time unit: \${rt_unit:-unspecified}"
    if [ "\${rt_unit}" = "second" ]; then
        echo "Converting Dinosaur retention times from minutes to seconds"
        awk -F'\\t' -v OFS='\\t' 'NR==1 {print; next} {\$4*=60; \$5*=60; \$6*=60; \$7*=60; print}' \\
            "dinosaur/${mzml.baseName}.features.tsv" > "dinosaur/.rt.tmp"
        mv "dinosaur/.rt.tmp" "dinosaur/${mzml.baseName}.features.tsv"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dinosaur: 1.2.1
        java: \$(java -version 2>&1 | head -n1 | sed 's/.*"\\(.*\\)".*/\\1/')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p dinosaur
    printf 'mz\\tmostAbundantMz\\tcharge\\trtStart\\trtApex\\trtEnd\\tfwhm\\tnIsotopes\\tnScans\\taveragineCorr\\tmass\\tmassCalib\\tintensityApex\\tintensitySum\\n' \\
        > dinosaur/${mzml.baseName}.features.tsv
    printf '500.25\\t500.25\\t2\\t10.0\\t10.5\\t11.0\\t0.2\\t3\\t12\\t0.99\\t998.5\\t998.5\\t1000\\t5000\\n' \\
        >> dinosaur/${mzml.baseName}.features.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dinosaur: 1.2.1
        java: stub
    END_VERSIONS
    """
}
