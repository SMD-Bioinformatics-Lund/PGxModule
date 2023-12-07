#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CHECK_INPUT                   } from '../subworkflows/local/create_meta'
include { PHARMACO_GENOMICS             } from '../subworkflows/local/pharmacoGenomics'
include { CUSTOM_DUMPSOFTWAREVERSIONS   } from '../modules/local/custom/dumpsoftwareversions/main'

csv = file(params.csv)

workflow PGX {

    take:
        input_pgx_target_beds
        input_pgx_target_rsids

    main:
        ch_versions = Channel.empty()

        CHECK_INPUT ( Channel.fromPath(csv) )

        PHARMACO_GENOMICS ( CHECK_INPUT.out.bam, input_pgx_target_beds, input_pgx_target_rsids )
        ch_versions = ch_versions.mix(PHARMACO_GENOMICS.out.versions)

        CUSTOM_DUMPSOFTWAREVERSIONS (
            ch_versions.unique().collectFile(name: 'collated_versions.yml'),
            CHECK_INPUT.out.meta
        )

    emit:
        report   = PHARMACO_GENOMICS.out.pgx_report
        versions = ch_versions
}