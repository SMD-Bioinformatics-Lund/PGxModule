#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CHECK_INPUT                   } from '../subworkflows/local/create_meta'
include { PHARMACO_GENOMICS             } from '../subworkflows/local/pharmacoGenomics'
include { CUSTOM_DUMPSOFTWAREVERSIONS   } from '../modules/local/custom/dumpsoftwareversions/main'

csv = file(params.csv)

workflow PGX {

    main:
        ch_versions = Channel.empty()

        CHECK_INPUT ( Channel.fromPath(csv) )

        PHARMACO_GENOMICS ( CHECK_INPUT.out.bam )
        ch_versions = ch_versions.mix(PHARMACO_GENOMICS.out.versions)

        
        samples =  CHECK_INPUT.out.meta.collect { it[0] }.map {concatenate_values(it) }

        CUSTOM_DUMPSOFTWAREVERSIONS (
            ch_versions.unique().collectFile(name: 'collated_versions.yml'),
            samples
        )

    emit:
        report   = PHARMACO_GENOMICS.out.pgx_report // channel: [ val(group), val(meta), file(pgx-report) ]
        versions = ch_versions                      // channel: [ file(versions) ]
}


def concatenate_values(channel) {
    channel.collect().join('_')
}