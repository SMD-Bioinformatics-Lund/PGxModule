process VARIANT_FILTRATION {
	publishDir "${params.outdir}/${params.subdir}/vcf", mode: 'copy', overwrite: true, pattern: "*.vcf"
	cpus 1
	time '1h'
	tag "$group"
	stageInMode 'copy'
	stageOutMode 'copy'
	container = "${params.python_image}"

	input:
		tuple val(group), file(vcf)

	output:
		tuple val(group), file("${group}.haplotypes.filtered.annotated.vcf"), emit: haplotypes_filtered 
	
	script:
	"""
    variant_filtration.py \
        --input_vcf=$vcf \
        --read_ratio=$params.read_ratio	\
        --depth=$params.DP \
        --output_file=${group}.haplotypes.filtered.annotated.vcf
	"""

	stub:
	"""
	touch ${group}.haplotypes.filtered.annotated.vcf
	"""

}