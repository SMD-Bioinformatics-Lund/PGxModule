process GET_CLIINICAL_GUIDELINES {
    // Given detected variants, get possible Haplotype combinations
    label 'process_single'
    label 'stage'
    tag "$meta.group"

    input:
        tuple val(group), val(meta), file(detected_variants)

    output:
        tuple val(group), val(meta), file("*.possible_diplotypes.tsv"), emit: possible_diplotypes
        path "versions.yml",                                            emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args   ?: ''
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """
        get_possible_diplotypes.py \
            --variant_csv $detected_variants \
            --output ${prefix}.possible_diplotypes.tsv  \
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """

    stub:
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """
        touch ${prefix}.possible_diplotypes.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """
}

process GET_INTERACTION_GUIDELINES {
    // Given Haplotype Combinations, get possible interactions betweens these
    label 'process_single'
    label 'stage'
    tag "$meta.group"

    input:
        tuple val(group), val(meta), file(possible_diplotypes)

    output:
        tuple val(group), val(meta), file("*.possible_interactions.tsv"),   emit: possible_interactions
        path "versions.yml",                                                emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args   ?: ''
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """
        get_interaction_guidelines.py \
            --diploids $possible_diplotypes \
            --output ${prefix}.possible_interactions.tsv \
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """

    stub:
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """
        touch ${prefix}.possible_interactions.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """
}

process GET_PGX_REPORT {
    // Generates markdown report per sample
    label 'process_single'
    label 'stage'
    tag "$meta.group"
    
    input:
        tuple val(group), val(meta), file(annotated_vcf), file(detected_variants), file(depth_at_missing_annotated_gdf), file(possible_diplotypes), file(depth_at_padded_baits), file(possible_interactions)

    output:
        tuple val(group), val(meta), file("*.pgx.html"),            emit: pgx_html
        tuple val(group), val(meta), file("*.targets.depth.tsv"),   emit: targets_depth
        path "versions.yml",                                        emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args   ?: ''
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """
        create_report.py \
            --group $group \
            --detected_variants $detected_variants \
            --missing_annotated_depth $depth_at_missing_annotated_gdf \
            --possible_diplotypes $possible_diplotypes \
            --possible_interactions $possible_interactions \
            --padded_baits_depth $depth_at_padded_baits \
            --annotated_vcf $annotated_vcf \
            --output ${prefix}.pgx.html \
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """

    stub:
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """
        touch ${prefix}.pgx.html
        touch ${prefix}.targets.depth.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """
    }