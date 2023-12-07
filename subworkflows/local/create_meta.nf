#!/usr/bin/env nextflow

// might need to add a check to csv? //
include { CSV_CHECK      } from '../../modules/local/check_input/main'

workflow CHECK_INPUT {
	take:
		csv		// file(csv)

	main:
		CSV_CHECK ( csv )
		checkedCsv = CSV_CHECK.out.csv.splitCsv( header:true, sep:',').set { csvmap }

		bam		= csvmap.map { create_bam_channel(it) }
		meta	= csvmap.map { create_samples_channel(it) }

	emit:
		bam        	// channel: [ val(meta), [ bam, bai ] ]
		meta        // channel: [ sample_id, sex, phenotype, paternal_id, maternal_id, case_id ]

}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_bam_channel(LinkedHashMap row) {
	// create meta map
	def meta = [:]
	meta.id					= row.id
	meta.group              = row.group
	meta.type               = row.type
	meta.clarity_sample_id  = (row.containsKey("clarity_sample_id") ? row.clarity_sample_id : None)
	meta.ffpe               = (row.containsKey("ffpe") ? row.ffpe : false)
	meta.purity             = (row.containsKey("purity") ? row.purity : false)
	// add path(s) of the fastq file(s) to the meta map
	def bam_meta = []
	// bam_meta = [row.group, row.id, meta, file(row.bam), file(row.bai) ]
	bam_meta = [row.group, row.id, row.type, file(row.bam), file(row.bai) ]

	return bam_meta
}

// Function to get a list of metadata (e.g. pedigree, case id) from the sample; [ meta ]
def create_samples_channel(LinkedHashMap row) {
	def meta                = [:]
	meta.id                 = row.id
	meta.group				= row.group
	meta.type               = row.type
	meta.clarity_sample_id  = (row.containsKey("clarity_sample_id") ? row.clarity_sample_id : None)
	meta.ffpe               = (row.containsKey("ffpe") ? row.ffpe : false)
	meta.purity             = (row.containsKey("purity") ? row.purity : false)
	def sample_meta = [row.group, meta]
	return sample_meta
}
