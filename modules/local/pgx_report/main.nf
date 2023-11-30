process GET_CLIINICAL_GUIDELINES {
    // Given detected variants, get possible Haplotype combinations
    publishDir "${params.outdir}/${params.subdir}/report/possible_diploids/", mode: 'copy', overwrite: true, pattern: "*.tsv"
	cpus 1
	time '1h'
	tag "$group"
	stageInMode 'copy'
	stageOutMode 'copy'
	container = "${params.python_image}"

	input:
		tuple val(group), file(detected_variants)

	output:
		tuple val(group), file("${group}.possible_diplotypes.tsv"), emit: possible_diplotypes

    script:
	"""
	get_possible_diplotypes.py \
		--variant_csv $detected_variants \
		--haplotype_definitions $params.haplotype_definitions \
		--clinical_guidelines $params.clinical_guidelines \
		--haplotype_activity $params.haplotype_activity \
		--output ${group}.possible_diplotypes.tsv \
		--hidden_haplotypes $params.hidden_haplotypes
	"""

    stub:
	"""
	touch ${group}.possible_diplotypes.tsv
	"""
}

process GET_INTERACTION_GUIDELINES {
    // Given Haplotype Combinations, get possible interactions betweens these
    publishDir "${params.outdir}/${params.subdir}/report/possible_interactions/", mode: 'copy', overwrite: true, pattern: "*.tsv"
	cpus 1
	time '1h'
	tag "$group"
	stageInMode 'copy'
	stageOutMode 'copy'
	container = "${params.python_image}"

	input:
		tuple val(group), file(possible_diplotypes)

	output:
		tuple val(group), file("${group}.possible_interactions.tsv"), emit: possible_interactions

    script:
	"""
	get_interaction_guidelines.py \
		--diploids $possible_diplotypes \
		--interaction_guidelines $params.interacting_guidelines \
		--output ${group}.possible_interactions.tsv
	"""

    stub:
	"""
	touch ${group}.possible_interactions.tsv
	"""
}

process GET_PGX_REPORT {
    // Generates markdown report per sample
    publishDir "${params.outdir}/${params.subdir}/report/", mode: 'copy', overwrite: true, pattern: "*.html"
	cpus 3
	time '1h'
	tag "$group"
	stageInMode 'copy'
	stageOutMode 'copy'
	container = "${params.jinja_image}"

	input:
        tuple val(group), file(annotated_vcf), file(detected_variants), file(depth_at_missing_annotated_gdf), file(possible_diplotypes), file(depth_at_padded_baits), file(possible_interactions), file(target_bed), file(target_rsids)

	output:
		tuple val(group), file("${group}.pgx.html"), emit: pgx_html

    script:
	"""
	create_report.py \
		--group $group \
		--read_depth $params.DP \
		--detected_variants $detected_variants \
		--missing_annotated_depth $depth_at_missing_annotated_gdf \
		--haplotype_definitions $params.haplotype_definitions \
		--possible_diplotypes $possible_diplotypes \
		--possible_interactions $possible_interactions \
		--target_bed $target_bed \
		--padded_baits_depth $depth_at_padded_baits \
		--target_rsids $target_rsids \
		--annotated_vcf $annotated_vcf \
		--dbSNP_version $params.dbSNP_version \
		--genome_version $params.ref_version \
		--output ${group}.pgx.html \
		--report_template $params.report_template \
		--report_css $params.report_css \

	"""

    stub:
	"""
	touch ${group}.pgx.html
    """
}