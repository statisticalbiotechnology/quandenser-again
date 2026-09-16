#!/usr/bin/env nextflow
/*
 * statisticalbiotechnology/quandenser-again
 *
 * Label-free MS1 quantification with match-between-runs.
 *
 *   nextflow run . -profile docker --input samplesheet.csv --outdir results
 */

nextflow.enable.dsl = 2

include { QUANDENSER_PIPELINE } from './workflows/quandenser'

//
// Read the samplesheet and fail early on the things that go wrong quietly.
//
def parseSamplesheet(csv_path) {
    channel
        .fromPath(csv_path, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            if (!row.sample || !row.mzml) {
                error("Samplesheet must have 'sample' and 'mzml' columns with non-empty values. Offending row: ${row}")
            }
            def meta = [
                id       : row.sample,
                condition: row.condition ?: 'A'
            ]
            [ meta, file(row.mzml, checkIfExists: true) ]
        }
}

workflow {

    if (!params.input) {
        error("No input given. Use --input <samplesheet.csv>. See docs/usage.md for the format.")
    }

    ch_samplesheet = parseSamplesheet(params.input)

    //
    // Quandenser identifies runs by the *stem* of the mzML file name and has a
    // standing TODO about two inputs sharing a basename in different folders
    // (src/Quandenser.cpp:555): both map to the same Dinosaur features file and
    // one silently overwrites the other. Catch it here rather than let it
    // corrupt the quantification.
    //
    ch_samplesheet
        .map { _meta, mzml -> mzml.baseName }
        .toList()
        .subscribe { names ->
            def dupes = names.countBy { name -> name }.findAll { _k, v -> v > 1 }.keySet()
            if (dupes) {
                error("Two or more input files share a base name: ${dupes.join(', ')}. " +
                      "Quandenser keys its intermediate files on the base name, so these would overwrite each other. Rename them.")
            }
            if (names.size() < 2) {
                error("Quandenser needs at least two runs to align; got ${names.size()}.")
            }
        }

    if (params.run_identification && !params.fasta) {
        error("--run_identification needs a protein database. Pass --fasta <file.fasta>.")
    }

    QUANDENSER_PIPELINE( ch_samplesheet )
}
