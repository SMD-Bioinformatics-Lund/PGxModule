process DETECTED_VARIANTS {
	publishDir "${params.outdir}/${params.subdir}/pgx/report/detected_variants", mode: 'copy', overwrite: true, pattern: "*.csv"
	cpus 1
	time '0.5h'
	tag "$group"
	stageInMode 'copy'
	stageOutMode 'copy'
	container = "${params.containers}/target_variants_python.simg"

	input:
		tuple val(group), file(vcf), file(pgx_ontarget_rsids_bed)
		
	output:
		tuple val(group), file("${group}.detected_variants.csv"), emit: detected_csv

	script:
	"""
    get_target_variants.py \
        --target_bed $pgx_ontarget_rsids_bed \
        --vcf $vcf \
        --output ${group}.detected_variants.csv \
		--addchr $params.addchr
	"""

	stub:
	"""
    get_target_variants.py \
        --target_bed $pgx_ontarget_rsids_bed \
        --vcf $vcf \
        --output ${group}.detected_variants.csv \
		--addchr $params.addchr
	"""
}