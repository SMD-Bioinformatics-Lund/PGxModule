#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PGX } from './workflows/pgx'

println(params.genome_file)
genome_file = file(params.genome_file)

OUTDIR = params.outdir+'/'+params.subdir
CRONDIR = params.crondir

params.csv = "/data/bnf/dev/ram/Pipelines/DSL2/PGxModule/test/data/test_sample_sheet.csv"
csv = file(params.csv)
println(csv)

params.vcf = "/data/bnf/dev/ram/Pipelines/DSL2/PGxModule/test/data/test_vcf_sheet.csv"
vcf = file(params.vcf)
println(vcf)


Channel
    .fromPath(params.csv).splitCsv(header:true)
    .map{ row-> tuple(row.group, row.id, row.type, file(row.bam), file(row.bai)) }
    .set { input_bam }

Channel
    .fromPath(params.vcf).splitCsv(header:true)
    .map{ row-> tuple(row.group, file(row.vcf)) }
    .set { input_vcfs }

Channel
    .fromPath("${params.pgx_target_regions}")
    .ifEmpty { exit 1, "PGX target regions bed file not found: ${params.pgx_target_regions}" }
    .set { input_pgx_target_beds }

Channel
    .fromPath("${params.pgx_target_rsid}")
    .ifEmpty { exit 1, "PGX target rsids bed file not found: ${params.pgx_target_rsid}" }
    .set { input_pgx_target_rsids }



workflow {

	PGX(input_bam, input_vcfs, input_pgx_target_beds, input_pgx_target_rsids)
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
	// logFile = file("/fs1/results/cron/logs/" + base + ".complete")
	logFile = file("/data/bnf/dev/ram/Pipelines/DSL2/PGxModule/logs/" + base + ".complete")
	logFile.text = msg
	logFile.append(error)
}