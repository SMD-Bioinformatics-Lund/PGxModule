process BCFTOOLS_ANNOTATION {
	publishDir "${params.outdir}/${params.subdir}/vcf", mode: 'copy', overwrite: true, pattern: "*.vcf"
	cpus 30
	time '3h'
	tag "$group"
	stageInMode 'copy'
	stageOutMode 'copy'
	container = "${params.bcftools_image}"
	memory '25 GB'

	input:
		tuple val(group), file(vcf), file(tbi)

	output:
		tuple val(group), file("${group}.sentieon.haplotypes.anno.vcf"), emit: annotations
		tuple val(group), file("${group}.${task.process.split(':').last()}.versions.yaml"), emit: versions

	script:
	def processName = task.process.toString().split(':').last()
	"""
    bcftools annotate --threads ${task.cpus} -a $params.dbSNP -c ID -o ${group}".sentieon.haplotypes.anno.vcf" $vcf

	{
	echo -e "${processName}:"
	echo -e "\tBCFTools:"
	echo -e "\t\tversion: \$(bcftools --version | grep 'bcftools' 2>&1 | sed 's/^.*bcftools //')"
	echo -e "\t\tcontainer: ${task.container}"
	} > "${group}.${processName}.versions.yaml"

	"""

	stub:
	def processName = task.process.toString().split(':').last()
	"""
	touch ${group}".sentieon.haplotypes.anno.vcf"

	{
	echo -e "${processName}:"
	echo -e "\tBCFTools:"
	echo -e "\t\tversion: \$(bcftools --version | grep 'bcftools' 2>&1 | sed 's/^.*bcftools //')"
	echo -e "\t\tcontainer: ${task.container}"
	} > "${group}.${processName}.versions.yaml"
	"""
}