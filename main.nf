#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PGX 							} from './workflows/pgx'

println(params.genome_file)

OUTDIR = params.outdir+'/'+params.subdir
CRONDIR = params.crondir

csv = file(params.csv)
println(csv)


Channel
    .fromPath("${params.pgx_target_regions}")
    .ifEmpty { exit 1, "PGX target regions bed file not found: ${params.pgx_target_regions}" }
    .set { input_pgx_target_beds }

Channel
    .fromPath("${params.pgx_target_rsid}")
    .ifEmpty { exit 1, "PGX target rsids bed file not found: ${params.pgx_target_rsid}" }
    .set { input_pgx_target_rsids }


workflow {

	PGX (
		input_pgx_target_beds, 
		input_pgx_target_rsids
	)

}


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
	logFile = file("${params.resultsdir}/cron/logs/" + base + "pgx.complete")
	logFile.text = msg
	logFile.append(error)
}
