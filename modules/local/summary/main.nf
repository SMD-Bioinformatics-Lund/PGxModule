process SAMPLE_TARGET_LIST {
    label 'process_single'
    label 'stage'
    tag "$group"

    input:
        tuple val(group), file(detected_variants), file(target_rsids)

    output:
        tuple val(group), file("${group}.pgx_target_interval.list"),    emit: pgx_target_interval_list
        path "versions.yml",                                            emit: versions

    script:
        """
        reform_genomic_region.py \
            --target_bed=$target_rsids \
            --output_file=${group}.pgx_target_interval.list \
            --detected_variants=$detected_variants \
            --addchr=$params.addchr

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """

    stub:
        """
        touch ${group}.pgx_target_interval.list

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
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
        tuple val(group), file("${group}.pgx_depth_at_missing.gdf"),    emit: pgx_depth_at_missing
        path "versions.yml",                                            emit: versions

    script:
        """
        java -jar /usr/GenomeAnalysisTK.jar -T DepthOfCoverage -allowPotentiallyMisencodedQuals -R $params.genome_file -I $ontarget_bam -o ${group}.pgx_depth_at_missing.gdf -L $target_interval_list

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk: \$(java -jar /usr/GenomeAnalysisTK.jar --version)
        END_VERSIONS
        """

    stub:
        """
        touch ${group}.pgx_depth_at_missing.gdf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk: \$(java -jar /usr/GenomeAnalysisTK.jar --version)
        END_VERSIONS
        """
}

process GET_PADDED_BAITS {
    label 'process_single'
    label 'stage'
    tag "Get_padded_baits_list"

    input:
        file(target_bed)

    output:
        path("padded_bait_interval.list"),  emit: padded_baits_list
        path "versions.yml",                emit: versions

    script:
        """
        reform_genomic_region.py \
            --target_bed=$target_bed \
            --output_file=padded_bait_interval.list \
            --padding=$params.padding \
            --addchr=$params.addchr

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """

    stub:
        """
        touch padded_bait_interval.list

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
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
        path "versions.yml",                        emit: versions

    script:
        """
        # NOTE: does not work with openjdk-11, openjdk-8 works
        java -jar /usr/GenomeAnalysisTK.jar -T DepthOfCoverage -allowPotentiallyMisencodedQuals -R $params.genome_file -I $ontarget_bam -o ${group}.pgx.gdf -L $padded_baits_interval_list

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk: \$(java -jar /usr/GenomeAnalysisTK.jar --version)
        END_VERSIONS
        """

    stub:
        """
        touch ${group}.pgx.gdf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk: \$(java -jar /usr/GenomeAnalysisTK.jar --version)
        END_VERSIONS
        """
}

process PADDED_BED_INTERVALS {
    label 'process_single'
    label 'stage'
    tag "pgx_bed_padded_intervals"

    input:
        file(target_bed)

    output:
        path('padded_bait_interval.bed'),   emit: padded_bed_intervals
        path "versions.yml",                emit: versions
        
    script:
        """
        reform_genomic_region.py \
            --target_bed=$target_bed \
            --output_file='padded_bait_interval.bed' \
            --padding=$params.padding \
            --format='bed' \
            --addchr=$params.addchr

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """

    stub:
        """
        touch padded_bait_interval.bed

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
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
        path "versions.yml",                                                    emit: versions

    script:
        """
        append_rsid_to_gdf.py \
            --input_gdf=$pgx_depth_at_missing_gdf \
            --target_bed=$target_rsids \
            --output_file=${group}.pgx_depth_at_missing_annotated.gdf \
            --addchr=$params.addchr

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """

    stub:
        """
        touch ${group}.pgx_depth_at_missing_annotated.gdf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """
}