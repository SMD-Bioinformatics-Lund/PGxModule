process DUMPSOFTWAREVERSIONS {
	publishDir "${params.outdir}/${params.subdir}/versions", mode: 'copy', overwrite: true, pattern: "*.versions.yml"
	cpus 1
	time '10m'
	tag "$group"
	stageInMode 'copy'
	stageOutMode 'copy'

    input:
    tuple val(group), file(ontarget_versions), file(haplotyper_versions), file(annotation_versions), file(depth_of_targets_versions), file(depth_of_baits_versions)

    output:
    tuple val(group), file("${group}.versions.yml"), emit: versions_yml

    script:
    """
    cat $ontarget_versions $haplotyper_versions $annotation_versions $depth_of_targets_versions $depth_of_baits_versions > ${group}.versions.yml
    """
}