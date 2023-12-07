process BCFTOOLS_ANNOTATION {
    label 'process_very_high'
    label 'stage'
    tag "$group"

    input:
        tuple val(group), file(vcf), file(tbi)

    output:
        tuple val(group), file("${group}.sentieon.haplotypes.anno.vcf"),    emit: annotations
        path "versions.yml",                                                emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        bcftools annotate --threads ${task.cpus} -a $params.dbSNP -c ID -o ${group}".sentieon.haplotypes.anno.vcf" $vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(bcftools --version | grep 'bcftools' 2>&1 | sed 's/^.*bcftools //')
        END_VERSIONS
        """

    stub:
        """
        touch ${group}".sentieon.haplotypes.anno.vcf"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(bcftools --version | grep 'bcftools' 2>&1 | sed 's/^.*bcftools //')
        END_VERSIONS
        """
}