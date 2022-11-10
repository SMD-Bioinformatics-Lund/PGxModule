process GET_CLIINICAL_GUIDELINES {
    // Given detected variants, get possible Haplotype combinations
    publishDir "${params.outdir}/${params.subdir}/pgx/report/possible_diploids/", mode: 'copy', overwrite: true, pattern: "*.csv"
	cpus 1
	time '1h'
	tag "$group"
	stageInMode 'copy'
	stageOutMode 'copy'
	container = "${params.containers}/target_variants_python.simg"

	input:
		tuple val(group), file(detected_variants)

	output:
		tuple val(group), file("${group}.possible_diplotypes.csv"), emit: possible_diplotypes

    script:
	"""
        get_possible_diplotypes.py \
            --variant_csv $detected_variants \
            --haplotype_definitions $params.haplotype_definitions \
            --clinical_guidelines $params.clinical_guidelines \
            --haplotype_activity $params.haplotype_activity \
            --output ${group}.possible_diplotypes.csv \
            --hidden_haplotypes $params.hidden_haplotypes
	"""

    stub:
	"""
        get_possible_diplotypes.py \
            --variant_csv $detected_variants \
            --haplotype_definitions $params.haplotype_definitions \
            --clinical_guidelines $params.clinical_guidelines \
            --haplotype_activity $params.haplotype_activity \
            --output ${group}.possible_diplotypes.csv \
            --hidden_haplotypes $params.hidden_haplotypes
	"""
}

process GET_INTERACTION_GUIDELINES {
    // Given Haplotype Combinations, get possible interactions betweens these
    publishDir "${params.outdir}/${params.subdir}/pgx/report/possible_interactions/", mode: 'copy', overwrite: true, pattern: "*.csv"
	cpus 1
	time '1h'
	tag "$group"
	stageInMode 'copy'
	stageOutMode 'copy'
	container = "${params.containers}/target_variants_python.simg"

	input:
		tuple val(group), file(possible_diplotypes)

	output:
		tuple val(group), file("${group}.possible_interactions.csv"), emit: possible_interactions

    script:
	"""
        get_interaction_guidelines.py \
            --diploids $possible_diplotypes \
            --interaction_guidelines $params.interacting_guidelines \
            --output ${group}.possible_interactions.csv
	"""

    stub:
	"""
        get_interaction_guidelines.py \
            --diploids $possible_diplotypes \
            --interaction_guidelines $params.interacting_guidelines \
            --output ${group}.possible_interactions.csv
	"""
}

process GET_PGX_REPORT {
    // Generates markdown report per sample
    publishDir "${params.outdir}/${params.subdir}/pgx/report/", mode: 'copy', overwrite: true, pattern: "*.html"
	cpus 3
	time '1h'
	tag "$group"
	stageInMode 'copy'
	stageOutMode 'copy'
	container = "${params.containers}/rmarkdown.sif"

	input:
        tuple val(group), file(annotated_vcf), file(detected_variants), file(depth_at_missing_annotated_gdf), file(possible_diplotypes), file(depth_at_padded_baits), file(possible_interactions), file(target_bed), file(target_rsids)

	output:
		tuple val(group), file("${group}.pgx.html"), emit: pgx_html

    script:
	"""
    bash runReport.sh $detected_variants $depth_at_missing_annotated_gdf $params.haplotype_definitions $possible_diplotypes $possible_interactions $group $target_bed $params.pgx_resources $depth_at_padded_baits $target_rsids \$PWD $params.scripts/report.Rmd ${group}.pgx.html $annotated_vcf
	"""

    stub:
	"""
    bash runReport.sh $detected_variants $depth_at_missing_annotated_gdf $params.haplotype_definitions $possible_diplotypes $possible_interactions $group $target_bed $params.pgx_resources $depth_at_padded_baits $target_rsids \$PWD $params.scripts/report.Rmd ${group}.pgx.html $annotated_vcf

    """
}