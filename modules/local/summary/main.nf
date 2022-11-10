process SAMPLE_TARGET_LIST {
	publishDir "${params.outdir}/${params.subdir}/pgx/report/coverage/", mode: 'copy', overwrite: true, pattern: "*.list"
	cpus 1
	time '1h'
	tag "$group"
	stageInMode 'copy'
	stageOutMode 'copy'
	container = "${params.containers}/target_variants_python.simg"

	input:
		tuple val(group), file(detected_variants), file(target_rsids)
        

	output:
		tuple val(group), file("${group}.pgx_target_interval.list"), emit: pgx_target_interval_list

	script:
	"""
        reform_genomic_region.py \
            --target_bed=$target_rsids \
            --output_file=${group}.pgx_target_interval.list \
            --detected_variants=$detected_variants \
			--addchr=$params.addchr
	"""

	stub:
	"""
        reform_genomic_region.py \
            --target_bed=$target_rsids \
            --output_file=${group}.pgx_target_interval.list \
            --detected_variants=$detected_variants \
			--addchr=$params.addchr
	"""
}

process DEPTH_OF_TARGETS {
    // Get read depth of variant locations at wildtrype-called positions
	publishDir "${params.outdir}/${params.subdir}/pgx/report/coverage/", mode: 'copy', overwrite: true, pattern: "*.gdf"
	cpus 2
	time '1h'
	tag "$group"
	stageInMode 'copy'
	stageOutMode 'copy'
	container = "${params.containers}/gatk3.simg"

	input:
		tuple val(group), file(ontarget_bam), file(ontarget_bai), file(target_interval_list)

	output:
		tuple val(group), file("${group}.pgx_depth_at_missing.gdf"), emit: pgx_depth_at_missing

	script:
	"""
        java -jar /usr/GenomeAnalysisTK.jar -T DepthOfCoverage -R $params.genome_file -I $ontarget_bam -o ${group}.pgx_depth_at_missing.gdf -L $target_interval_list
	"""

	stub:
	"""
        java -jar /usr/GenomeAnalysisTK.jar -T DepthOfCoverage -R $params.genome_file -I $ontarget_bam -o ${group}.pgx_depth_at_missing.gdf -L $target_interval_list
	"""
}

process GET_PADDED_BAITS {
	publishDir "${params.outdir}/${params.subdir}/pgx/gdf/", mode: 'copy', overwrite: true, pattern: "*.list"
	cpus 1
	time '1h'
	tag "Get_padded_baits_list"
	stageInMode 'copy'
	stageOutMode 'copy'
	container = "${params.containers}/target_variants_python.simg"

	input:
        file(target_bed)

	output:
		path("padded_bait_interval.list"), emit: padded_baits_list

	script:
	"""
        reform_genomic_region.py \
            --target_bed=$target_bed \
            --output_file=padded_bait_interval.list \
            --padding=$params.padding \
			--addchr=$params.addchr
	"""

	stub:
	"""
        reform_genomic_region.py \
            --target_bed=$target_bed \
            --output_file=padded_bait_interval.list \
            --padding=$params.padding \
			--addchr=$params.addchr
	"""
}

process DEPTH_OF_BAITS {
    // Get read depth of baits
	publishDir "${params.outdir}/${params.subdir}/pgx/gdf/", mode: 'copy', overwrite: true, pattern: "*.gdf"
	cpus 2
	time '1h'
	tag "$group"
	stageInMode 'copy'
	stageOutMode 'copy'
	container = "${params.containers}/gatk3.simg"

	input:
        tuple val(group), file(ontarget_bam), file(ontarget_bai), file(padded_baits_interval_list)

	output:
		tuple val(group), file("${group}.pgx.gdf"), emit: padded_baits_list

	script:
	"""
        # NOTE: does not work with openjdk-11, openjdk-8 works
        java -jar /usr/GenomeAnalysisTK.jar -T DepthOfCoverage -R $params.genome_file -I $ontarget_bam -o ${group}.pgx.gdf -L $padded_baits_interval_list
	"""

	stub:
	"""
        # NOTE: does not work with openjdk-11, openjdk-8 works
        java -jar /usr/GenomeAnalysisTK.jar -T DepthOfCoverage -R $params.genome_file -I $ontarget_bam -o ${group}.pgx.gdf -L $padded_baits_interval_list
	"""
}

process PADDED_BED_INTERVALS {
	publishDir "${params.outdir}/${params.subdir}/pgx/bam/", mode: 'copy', overwrite: true, pattern: "*.gdf"
	cpus 1
	time '1h'
	tag "pgx_bed_padded_intervals"
	stageInMode 'copy'
	stageOutMode 'copy'
	container = "${params.containers}/target_variants_python.simg"

	input:
		file(target_bed)

	output:
		path('padded_bait_interval.bed'), emit: padded_bed_intervals
	
	script:
	"""
    reform_genomic_region.py \
        --target_bed=$target_bed \
        --output_file='padded_bait_interval.bed' \
        --padding=$params.padding \
        --format='bed' \
		--addchr=$params.addchr
	"""

	stub:
	"""
    reform_genomic_region.py \
        --target_bed=$target_bed \
        --output_file='padded_bait_interval.bed' \
        --padding=$params.padding \
        --format='bed' \
		--addchr=$params.addchr
	"""
}

process APPEND_ID_TO_GDF {
	//  Add variant id to appropriate location in gdf
	publishDir "${params.outdir}/${params.subdir}/pgx/report/coverage/", mode: 'copy', overwrite: true, pattern: "*.gdf"
	cpus 1
	time '1h'
	tag "$group"
	stageInMode 'copy'
	stageOutMode 'copy'
	container = "${params.containers}/target_variants_python.simg"

	input:
		tuple val(group), file(pgx_depth_at_missing_gdf), file(target_rsids)

	output:
		tuple val(group), file("${group}.pgx_depth_at_missing_annotated.gdf"),  emit: depth_at_missing_annotate_gdf

	script:
	"""
        append_rsid_to_gdf.py \
        	--input_gdf=$pgx_depth_at_missing_gdf \
            --target_bed=$target_rsids \
            --output_file=${group}.pgx_depth_at_missing_annotated.gdf \
			--addchr=$params.addchr
	"""

	stub:
	"""
        append_rsid_to_gdf.py \
        	--input_gdf=$pgx_depth_at_missing_gdf \
            --target_bed=$target_rsids \
            --output_file=${group}.pgx_depth_at_missing_annotated.gdf \
			--addchr=$params.addchr
	"""
}