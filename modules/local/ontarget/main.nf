process ONTARGET_BAM {
	cpus 2
	time '1h'
	tag "$id"
	stageInMode 'copy'
	stageOutMode 'copy'
	container = "${params.containers}/samtools.simg"

	input:
		// tuple val(group), val(id), val(type), file("${id}.${type}.dedup.bam"), file("${id}.${type}.dedup.bam.bai")
		tuple val(group), val(id), val(type), file(bam), file(bai)
		file(pgx_ontarget_padded_bed)

	output:
		tuple val(group), val(id), val(type), file("${id}.${type}.dedup.ontarget.bam"), file("${id}.${type}.dedup.ontarget.bam.bai"), emit: bam_ontarget

	script:
	"""
	samtools view -h -b $bam -L $pgx_ontarget_padded_bed > ${id}.${type}.dedup.ontarget.bam
    samtools index ${id}.${type}.dedup.ontarget.bam
	"""

	stub:
	"""
	samtools view -h -b $bam -L $pgx_ontarget_padded_bed > ${id}.${type}.dedup.ontarget.bam
    samtools index ${id}.${type}.dedup.ontarget.bam
	"""
}

