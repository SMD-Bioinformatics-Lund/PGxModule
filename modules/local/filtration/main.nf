process VARIANT_FILTRATION {
	cpus 1
	time '1h'
	tag "$group"
	stageInMode 'copy'
	stageOutMode 'copy'
	container = "${params.containers}/target_variants_python.simg"

	input:
		tuple val(group), file(vcf)

	output:
		tuple val(group), file("${group}.GATK.filtered.vcf"), emit: filtered_vcf
	
	script:
	"""
    variant_filtration.py \
        --input_vcf=$vcf \
        --read_ratio=$params.read_ratio	\
        --depth=$params.DP \
        --output_file=${group}.GATK.filtered.vcf
	"""

	stub:
	"""
    variant_filtration.py \
        --input_vcf=$vcf \
        --read_ratio=$params.read_ratio	\
        --depth=$params.DP \
        --output_file=${group}.GATK.filtered.vcf
	"""

}