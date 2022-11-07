process ANNOTATE_VARIANTS {
	publishDir "${params.outdir}/${params.subdir}/pgx/annotated", mode: 'copy', overwrite: true, pattern: "*.vcf"
	cpus 1
	time '1h'
	tag "$group"
	stageInMode 'copy'
	stageOutMode 'copy'
	container = "${params.containers}/target_variants_python.simg"

	input:
		tuple val(group), file(vcf)

	output:
		tuple val(group), file("${group}.pgx_annotated.vcf"), emit: annotated_vcf
	
	script:
	"""
    modifyVCF.py \
        --input $vcf \
        --addrsid True \
        --output ${group}.pgx_annotated.vcf
	"""

	stub:
	"""
    modifyVCF.py \
        --input $vcf \
        --addrsid True \
        --output ${group}.pgx_annotated.vcf
	"""

}