#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CHECK_INPUT                   } from '../subworkflows/local/create_meta'
include { PHARMACO_GENOMICS             } from '../subworkflows/local/pharmacoGenomics'
include { CUSTOM_DUMPSOFTWAREVERSIONS   } from '../modules/local/custom/dumpsoftwareversions/main'
include { MULTIQC                       } from '../modules/local/multiqc/main'

csv = file(params.csv)

// Config File Channels //

ch_multiqc_config                       = Channel.fromPath(params.multiqc_config, checkIfExists: true)
ch_multiqc_custom_config                = params.multiqc_extra_config           ? Channel.fromPath( params.multiqc_extra_config, checkIfExists: true )  : Channel.empty()
ch_multiqc_logo                         = params.multiqc_logo                   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true )          : Channel.empty()
ch_multiqc_custom_methods_description   = params.multiqc_methods_desc           ? file(params.multiqc_methods_desc, checkIfExists: true)                : Channel.empty()


workflow PGX {

    main:
        ch_versions = Channel.empty()

        //
        // SUBWORKFLOW: CHECK_INPUT
        //
        CHECK_INPUT ( Channel.fromPath(csv) )

        //
        // SUBWORKFLOW: PHARMACO_GENOMICS
        //
        PHARMACO_GENOMICS ( CHECK_INPUT.out.bam )
        ch_versions = ch_versions.mix(PHARMACO_GENOMICS.out.versions)

        //
        // MODULE: Software Versions
        //
        samples =  CHECK_INPUT.out.meta.collect { it[0] }.map {concatenate_values(it) }
        CUSTOM_DUMPSOFTWAREVERSIONS (
            ch_versions.unique().collectFile(name: 'collated_versions.yml'),
            samples
        )

        //
        // MODULE: MultiQC
        //

        ch_multiqc_files    = Channel.empty()
        ch_multiqc_files    = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
        // ch_multiqc_files    = ch_multiqc_files.mix(PHARMACO_GENOMICS.out.targets_depth.map{it[2]}.collect().ifEmpty([]))

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList(),
            PHARMACO_GENOMICS.out.targets_depth
        )
        multiqc_reports = MULTIQC.out.report.toList()


    emit:
        report                  = PHARMACO_GENOMICS.out.pgx_report              // channel: [ val(group), val(meta), file(pgx-report) ]
        // pharmcat_pgx_regions    = PHARMACO_GENOMICS.out.pharmcat_pgx_regions    // channel: [ tuple val(group), val(meta) file("pharmcat_preprocessed.vcf") ]
        pharmcat_preprocessed   = PHARMACO_GENOMICS.out.pharmcat_preprocessed   // channel: [ tuple val(group), val(meta) file("pharmcat_preprocessed.vcf") ]
        pharmcat_report         = PHARMACO_GENOMICS.out.pharmcat_report         // channel: [ tuple val(group), val(meta) file("pharmcat.html") ]
        pharmcat_pheno_json     = PHARMACO_GENOMICS.out.pharmcat_pheno_json     // channel: [ tuple val(group), val(meta) file("pharmcat.phenotype.json") ]
        pharmcat_match_json     = PHARMACO_GENOMICS.out.pharmcat_match_json     // channel: [ tuple val(group), val(meta) file("pharmcat.match.json") ]
        pharmcat_match_html     = PHARMACO_GENOMICS.out.pharmcat_match_html     // channel: [ tuple val(group), val(meta) file("pharmcat.match.html") ]
        pharmcat_report_json     = PHARMACO_GENOMICS.out.pharmcat_report_json   // channel: [ tuple val(group), val(meta) file("pharmcat.match.html") ]
        multiqc_report          = multiqc_reports                               // channel: [ tuple val(group), val(meta) file("multiqc.html") ]
        multiqc_data            = MULTIQC.out.data                              // channel: [ tuple val(group), val(meta) file("multiqc_data") ]
        multiqc_plots           = MULTIQC.out.plots                             // channel: [ tuple val(group), val(meta) file("multiqc_plots") ]
        versions                = ch_versions                                   // channel: [ file(versions) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow.onComplete {

	def msg = """\
		Pipeline execution summary
		---------------------------
		Completed at: ${workflow.complete}
		Duration    : ${workflow.duration}
		Success     : ${workflow.success}
		scriptFile  : ${workflow.scriptFile}
		workDir     : ${workflow.workDir}
		exit status : ${workflow.exitStatus}
		errorMessage: ${workflow.errorMessage}
		errorReport :
		"""
		.stripIndent()
	def error = """\
		${workflow.errorReport}
		"""
		.stripIndent()

	base = csv.getBaseName()
	logFile = file("${params.crondir}/logs/" + base + "pgx.complete")
	logFile.text = msg
	logFile.append(error)
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
def concatenate_values(channel) {
    channel.collect().join('_')
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/