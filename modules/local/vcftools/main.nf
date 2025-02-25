process VCFTOOLS_FILTER {
    tag "${meta.group}"

    input:
        tuple val(group), val(meta), file(vcf)

    output:
        tuple val(group), val(meta), file("*.filtered.vcf"),    emit: filtered_vcf
        path "versions.yml",                                    emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args ?: ''
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """ 
        vcffilter $args ${vcf} \\
        | vcffilter $args2 \\
        | vcfglxgt > ${prefix}.filtered.tagged.vcf

        vcftools --vcf ${prefix}.filtered.tagged.vcf --out ${prefix}.filtered.vcf $args3

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vcffilter: \$(echo \$( vcffilter -h 2>&1) | grep 'vcflib' | sed 's/ filter.*\$//g' | sed 's/.* //g' )
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.group}"
        """
        touch ${prefix}.filtered.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vcffilter: \$(echo \$( vcffilter -h 2>&1) | grep 'vcflib' | sed 's/ filter.*\$//g' | sed 's/.* //g' )
        END_VERSIONS
        """
}