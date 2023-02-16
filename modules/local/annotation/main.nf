process ANNOTATE_VARIANTS {
	publishDir "${params.outdir}/${params.subdir}/pgx/annotated", mode: 'copy', overwrite: true, pattern: "*.vcf"
	cpus 1
	time '1h'
	tag "$group"
	stageInMode 'copy'
	stageOutMode 'copy'
	container = "${params.containers}/gatk4.simg"

	input:
		tuple val(group), file(bam), file(bai), file(vcf)

	output:
		tuple val(group), file("${group}.pgx_annotated.vcf"), emit: annotated_vcf
	
	script:
	"""
    gatk VariantAnnotator \
        -R $params.genome_file \
        -V $vcf \
        -I $bam \
        -O ${group}.pgx_annotated.vcf \
        --dbsnp $params.dbSNP
	"""

	stub:
	"""
    gatk VariantAnnotator \
        -R $params.genome_file \
        -V $vcf \
        -I $bam \
        -O ${group}.pgx_annotated.vcf \
        --dbsnp $params.dbSNP
	"""

}