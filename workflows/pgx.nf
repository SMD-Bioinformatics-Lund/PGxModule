#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PHARMACO_GENOMICS } from '../subworkflows/local/pharmacoGenomics'

workflow PGX {

	take:
		input_bam
		input_pgx_target_beds
		input_pgx_target_rsids

	main:
		PHARMACO_GENOMICS ( input_bam, input_pgx_target_beds, input_pgx_target_rsids )
		.set { pgx_reports }

	emit:
        report = PHARMACO_GENOMICS.out.pgx_report
}