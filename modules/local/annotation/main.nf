process BCFTOOLS_ANNOTATION {
    label 'process_very_high'
    label 'stage'
    tag "$meta.group"

    input:
        tuple val(group), val(meta), file(vcf), file(tbi)

    output:
        tuple val(group), val(meta), file("*.haplotypes.anno.vcf"), emit: annotations
        path "versions.yml",                                        emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args   ?: ''
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """
        bcftools annotate --threads ${task.cpus} $args -o ${prefix}".filtered.ontarget.haplotypes.anno.vcf" $vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(bcftools --version | grep 'bcftools' 2>&1 | sed 's/^.*bcftools //')
        END_VERSIONS
        """

    stub:
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """
        touch ${prefix}".filtered.ontarget.haplotypes.anno.vcf"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(bcftools --version | grep 'bcftools' 2>&1 | sed 's/^.*bcftools //')
        END_VERSIONS
        """
}