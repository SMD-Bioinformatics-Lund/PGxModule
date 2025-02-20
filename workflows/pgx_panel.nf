#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PREPROCESS_FASTQ              } from '../subworkflows/local/preprocess_fastq'
include { SAMPLE                        } from '../subworkflows/local/sample'
include { ALIGN_SENTIEON                } from '../subworkflows/local/alignment'
include { QC                            } from '../subworkflows/local/qc'
include { PADDED_INTERVALS              } from '../subworkflows/local/padded_intervals'
include { HAPLOTYPING                   } from '../subworkflows/local/haplotyping'
include { ONTARGET                      } from '../subworkflows/local/ontarget'
include { ANNOTATION                    } from '../subworkflows/local/annotations'
include { PHARMCAT                      } from '../subworkflows/local/pharmcat'
include { COVERAGE                      } from '../subworkflows/local/coverage'
include { CLINICAL_INFORMATION          } from '../subworkflows/local/clinical_information'
include { PGX_REPORT                    } from '../subworkflows/local/reports'
include { CUSTOM_DUMPSOFTWAREVERSIONS   } from '../modules/local/custom/dumpsoftwareversions/main'
// include { TEST_FASTQ   } from '../modules/local/test/main'
// include { TEST_BAM   } from '../modules/local/test/main'
include { RUN_MULTIQC                   } from '../subworkflows/local/multiqc'

csv = file(params.csv)

workflow PGX_PANEL {

    take:
        fastq_input                 // FASTQ CHANNEL for PGX_FULL
        bam_input
        samples


    main:
        // TEST_FASTQ ( fastq_input )
        // TEST_BAM ( bam_input )

        // Preprocessing Fastq (subsample and trim)
        PREPROCESS_FASTQ ( fastq_input ).set { fastq_processed }

        // Alignment with Sentieon
        ALIGN_SENTIEON ( fastq_processed.proccessed_fastq ).set { bam_aligned }

        // Basically it will either have input from csv or from align sub workflow
        bam = bam_input.mix(bam_aligned.bam)

    // emit:
    //     fastq_trimmed   =   TEST_FASTQ.out
    //     bam_trimmed     =   TEST_BAM.out
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