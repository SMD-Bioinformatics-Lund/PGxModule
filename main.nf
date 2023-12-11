#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PGX 							} from './workflows/pgx'

println(params.genome_file)
csv = file(params.csv)
println(csv)


workflow {

	PGX ()
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
	logFile = file("${params.crondir}/logs/" + base + "pgx.complete")
	logFile.text = msg
	logFile.append(error)
}
