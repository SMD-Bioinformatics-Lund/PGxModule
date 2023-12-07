process VARIANT_FILTRATION {
    label 'process_single'
    label 'stage'
    tag "$group"

    input:
        tuple val(group), file(vcf)

    output:
        tuple val(group), file("${group}.haplotypes.filtered.annotated.vcf"),   emit: haplotypes_filtered
        path "versions.yml",                                                    emit: versions
        
    script:
        """
        variant_filtration.py \
            --input_vcf=$vcf \
            --read_ratio=$params.read_ratio \
            --depth=$params.DP \
            --output_file=${group}.haplotypes.filtered.annotated.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """

    stub:
        """
        touch ${group}.haplotypes.filtered.annotated.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """

}