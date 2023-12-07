process SAMPLE_TARGET_LIST {
    label 'process_single'
    label 'stage'
    tag "$group"

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
    touch ${group}.pgx_target_interval.list
    """
}

process DEPTH_OF_TARGETS {
    // Get read depth of variant locations at wildtrype-called positions
    label 'process_low'
    label 'stage'
    tag "$group"

    input:
        tuple val(group), file(ontarget_bam), file(ontarget_bai), file(target_interval_list)

    output:
        tuple val(group), file("${group}.pgx_depth_at_missing.gdf"), emit: pgx_depth_at_missing
        tuple val(group), file("${group}.${task.process.split(':').last()}.versions.yaml"), emit: versions

    script:
    def processName = task.process.toString().split(':').last()
    """
    java -jar /usr/GenomeAnalysisTK.jar -T DepthOfCoverage -allowPotentiallyMisencodedQuals -R $params.genome_file -I $ontarget_bam -o ${group}.pgx_depth_at_missing.gdf -L $target_interval_list

    {
        echo -e "${processName}:"
        echo -e "\tGATK DepthOfCoverage:"
        echo -e "\t\tversion: \$(java -jar /usr/GenomeAnalysisTK.jar --version)"
        echo -e "\t\tcontainer: ${task.container}"
    } > "${group}.${processName}.versions.yaml"
    """

    stub:
    def processName = task.process.toString().split(':').last()
    """
    touch ${group}.pgx_depth_at_missing.gdf

    {
        echo -e "${processName}:"
        echo -e "\tGATK DepthOfCoverage:"
        echo -e "\t\tversion: \$(java -jar /usr/GenomeAnalysisTK.jar --version)"
        echo -e "\t\tcontainer: ${task.container}"
    } > "${group}.${processName}.versions.yaml"
    """
}

process GET_PADDED_BAITS {
    label 'process_single'
    label 'stage'
    tag "Get_padded_baits_list"

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
    touch padded_bait_interval.list
    """
}

process DEPTH_OF_BAITS {
    // Get read depth of baits
    label 'process_low'
    label 'stage'
    tag "$group"

    input:
        tuple val(group), file(ontarget_bam), file(ontarget_bai), file(padded_baits_interval_list)

    output:
        tuple val(group), file("${group}.pgx.gdf"), emit: padded_baits_list
        tuple val(group), file("${group}.${task.process.split(':').last()}.versions.yaml"), emit: versions

    script:
    def processName = task.process.toString().split(':').last()
    """
    # NOTE: does not work with openjdk-11, openjdk-8 works
    java -jar /usr/GenomeAnalysisTK.jar -T DepthOfCoverage -allowPotentiallyMisencodedQuals -R $params.genome_file -I $ontarget_bam -o ${group}.pgx.gdf -L $padded_baits_interval_list

    {
        echo -e "${processName}:"
        echo -e "\tGATK DepthOfCoverage:"
        echo -e "\t\tversion: \$(java -jar /usr/GenomeAnalysisTK.jar --version)"
        echo -e "\t\tcontainer: ${task.container}"
    } > "${group}.${processName}.versions.yaml"
    """

    stub:
    def processName = task.process.toString().split(':').last()
    """
    touch ${group}.pgx.gdf

    {
        echo -e "${processName}:"
        echo -e "\tGATK DepthOfCoverage:"
        echo -e "\t\tversion: \$(java -jar /usr/GenomeAnalysisTK.jar --version)"
        echo -e "\t\tcontainer: ${task.container}"
    } > "${group}.${processName}.versions.yaml"

    """
}

process PADDED_BED_INTERVALS {
    label 'process_single'
    label 'stage'
    tag "pgx_bed_padded_intervals"

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
    touch padded_bait_interval.bed
    """
}

process APPEND_ID_TO_GDF {
    //  Add variant id to appropriate location in gdf
    label 'process_single'
    label 'stage'
    tag "$group"

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
    touch ${group}.pgx_depth_at_missing_annotated.gdf
    """
}