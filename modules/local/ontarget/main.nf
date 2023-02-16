process ONTARGET_BAM {
	cpus 2
	time '1h'
	tag "$id"
	stageInMode 'copy'
	stageOutMode 'copy'
	container = "${params.containers}/samtools.simg"

	input:
		tuple val(group), val(id), val(type), file(bam), file(bai), file(pgx_ontarget_padded_bed)

	output:
		tuple val(group), file("${group}.dedup.ontarget.bam"), file("${group}.dedup.ontarget.bam.bai"), emit: bam_ontarget

	script:
	"""
	samtools view -h -b $bam -L $pgx_ontarget_padded_bed -M > ${group}.dedup.ontarget.bam
    samtools index ${group}.dedup.ontarget.bam
	"""

	stub:
	"""
	samtools view -h -b $bam -L $pgx_ontarget_padded_bed -M > ${group}.dedup.ontarget.bam
    samtools index ${group}.dedup.ontarget.bam
	"""
}


process HAPLOTYPECALLING {
	cpus 20
	time '2h'
	tag "$group"
	stageInMode 'copy'
	stageOutMode 'copy'
	container = "${params.containers}/gatk4.simg"

	input:
		tuple val(group), file(bam), file(bai)

	output:
		tuple val(group), file("${group}.GATK.haplotype.vcf"), emit: ontarget_vcf

	script:
	"""
	gatk HaplotypeCaller -R $params.genome_file -I $bam --dbsnp $params.dbSNP -O ${group}.GATK.haplotype.vcf
	"""

	stub:
	"""
	gatk HaplotypeCaller -R $params.genome_file -I $bam --dbsnp $params.dbSNP -O ${group}.GATK.haplotype.vcf
	"""
}

