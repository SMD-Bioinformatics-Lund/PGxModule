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
        python3 $params.scripts/get_possible_diplotypes.py \
            --variant_csv $detected_variants \
            --haplotype_definitions $params.haplotype_definitions \
            --clinical_guidelines $params.clinical_guidelines \
            --haplotype_activity $params.haplotype_activity \
            --output ${group}.possible_diplotypes.csv \
            --hidden_haplotypes $params.hidden_haplotypes
	"""

    stub:
	"""
        python3 $params.scripts/get_possible_diplotypes.py \
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
        python3 $params.scripts/get_interaction_guidelines.py \
            --diploids $possible_diplotypes \
            --interaction_guidelines $params.interacting_guidelines \
            --output ${group}.possible_interactions.csv
	"""

    stub:
	"""
        python3 $params.scripts/get_interaction_guidelines.py \
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
		tuple val(group), file(detected_variants)
        tuple val(group), file(depth_at_missing_annotated_gdf)
        tuple val(group), file(possible_diplotypes)
        tuple val(group), file(depth_at_padded_baits)
        tuple val(group), file(possible_interactions)

	output:
		tuple val(group), file("${group}.pgx.html"), emit: pgx_html

    script:
	"""
        wkdir="\$PWD"  # Needed since Rscript will set wd to location of file not session
        intdir=\$(echo {output.html} | head -c -6)
        Rscript \
            -e ".libPaths('/lib/rlib'); library(rmdformats); rmarkdown::render('{$params.scripts/generate_sample_report.Rmd', output_file='${group}.pgx.html', output_format=c('readthedown'), intermediates_dir='$wkdir/$intdir')" \
            --args --title='Farmakogenomisk analys av $group' --author=Ram \
            --found_variants=$detected_variants \
            --missed_variants=$depth_at_missing_annotated_gdf  \
            --haplotype_definitions=$params.haplotype_definitions \
            --clinical_guidelines=$possible_diplotypes \
            --interaction_guidelines=$possible_interactions \
            --data_location={params.script_location}/data \
            --depth_file=$depth_at_padded_baits \
            --sample=$group \
            --seqid={wildcards.seqID} \
            --dbsnp=dbsnp147 \
            --ref=$params.genome_file \
            --name=Ram \
            --adress=Lund \
            --mail=Ram.Nanduri@skane.se \
            --phone=0123456789

            rmdir $wkdir/$intdir
	"""

    stub:
	"""
    touch ${group}.pgx.html
    """
}