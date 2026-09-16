/*
 * Main analysis workflow.
 */

include { DINOSAUR          } from '../modules/local/dinosaur/main'
include { QUANDENSER        } from '../modules/local/quandenser/main'
include { CRUX_TIDE_INDEX   } from '../modules/local/crux/tideindex/main'
include { CRUX_TIDE_SEARCH  } from '../modules/local/crux/tidesearch/main'
include { CRUX_PERCOLATOR   } from '../modules/local/crux/percolator/main'
include { TRIQLER_CONVERT   } from '../modules/local/triqler/convert/main'
include { TRIQLER           } from '../modules/local/triqler/main'

workflow QUANDENSER_PIPELINE {

    take:
    ch_samplesheet   // channel: [ val(meta), path(mzml) ]

    main:
    ch_versions = channel.empty()

    //
    // Build the batch file Quandenser reads with --batch.
    //
    // Two things make this more delicate than it looks.
    //
    // First, the line order of this file *is* the file index. MaRaCluster's
    // SpectrumFileList assigns each run its fileIdx by position, and column 0
    // of Quandenser.feature_groups.tsv is that integer rather than a file
    // name. Sorting the lines makes the mapping reproducible across runs
    // instead of depending on the order channel items happen to arrive in.
    //
    // Second, MaRaCluster reads only the first tab-separated field of each
    // line while Triqler reads the first two, so a single `path <TAB>
    // condition` file serves as both the Quandenser batch file and the Triqler
    // file list. That is why there is one file here and not two that could
    // drift apart.
    //
    // Published, because it is the only record of which integer in column 0 of
    // the feature-groups file corresponds to which run.
    ch_batch_file = ch_samplesheet
        .map { meta, mzml -> "${mzml.name}\t${meta.condition}" }
        .collectFile(
            name: 'batch.tsv',
            newLine: true,
            sort: true,
            storeDir: "${params.outdir}/pipeline_info"
        )

    ch_mzmls = ch_samplesheet
        .map { _meta, mzml -> mzml }
        .collect()

    //
    // MS1 feature detection, one task per run.
    //
    // This is the parallelism that matters. Quandenser runs Dinosaur in a
    // serial loop over inputs, each iteration a fresh JVM, and then skips any
    // run whose features file already exists. Computing them here and staging
    // them back in turns the longest stage from N sequential JVM invocations
    // into N concurrent ones, without changing a line of Quandenser.
    //
    if (params.parallel_feature_detection) {
        DINOSAUR( ch_samplesheet )
        ch_versions = ch_versions.mix(DINOSAUR.out.versions.first())
        ch_features = DINOSAUR.out.features.map { _meta, tsv -> tsv }.collect()
    }
    else {
        // Let the monolith do its own feature detection. Slower, but it is the
        // reference behaviour and a useful control when results are questioned.
        ch_features = channel.value([])
    }

    //
    // Clustering, alignment, match-between-runs, consensus spectra.
    //
    QUANDENSER( ch_mzmls, ch_batch_file, ch_features )
    ch_versions = ch_versions.mix(QUANDENSER.out.versions)

    //
    // Optional: identify the consensus spectra and quantify with Triqler.
    //
    ch_proteins = channel.empty()
    ch_triqler_input = channel.empty()

    if (params.run_identification) {
        CRUX_TIDE_INDEX( channel.fromPath(params.fasta, checkIfExists: true) )
        ch_versions = ch_versions.mix(CRUX_TIDE_INDEX.out.versions)

        CRUX_TIDE_SEARCH( QUANDENSER.out.consensus_spectra.collect(), CRUX_TIDE_INDEX.out.index )
        ch_versions = ch_versions.mix(CRUX_TIDE_SEARCH.out.versions)

        CRUX_PERCOLATOR( CRUX_TIDE_SEARCH.out.psms )
        ch_versions = ch_versions.mix(CRUX_PERCOLATOR.out.versions)

        TRIQLER_CONVERT(
            QUANDENSER.out.feature_groups,
            ch_batch_file,
            CRUX_PERCOLATOR.out.target_psms,
            CRUX_PERCOLATOR.out.decoy_psms
        )
        ch_versions = ch_versions.mix(TRIQLER_CONVERT.out.versions)
        ch_triqler_input = TRIQLER_CONVERT.out.triqler_input

        if (!params.skip_triqler) {
            TRIQLER( TRIQLER_CONVERT.out.triqler_input )
            ch_versions = ch_versions.mix(TRIQLER.out.versions)
            ch_proteins = TRIQLER.out.proteins
        }
    }

    ch_versions
        .unique()
        .collectFile(name: 'software_versions.yml', storeDir: "${params.outdir}/pipeline_info")

    emit:
    feature_groups    = QUANDENSER.out.feature_groups
    consensus_spectra = QUANDENSER.out.consensus_spectra
    triqler_input     = ch_triqler_input
    proteins          = ch_proteins
    versions          = ch_versions
}
