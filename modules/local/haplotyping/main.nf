process GATK_HAPLOTYPING {
    label 'process_medium_cpus'
    label 'process_medium_memory'
    label 'stage'
    tag "$meta.group"

    input:
        tuple val(group), val(meta), file(bam), file(bai)

    output:
        tuple val(group), val(meta), file("*.haplotypes.vcf"),      emit: haplotypes
        path "versions.yml",                                        emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args        = task.ext.args   ?: ''
        def prefix      = task.ext.prefix ?: "${meta.group}.GATK"
        def avail_mem   = 12288
        if (!task.memory) {
            log.info '[GATK CollectAllelicCounts] Available memory not known - defaulting to 12GB. Specify process memory requirements to change this.'
        } else {
            avail_mem = (task.memory.mega*0.8).intValue()
        }
        """
        gatk --java-options "-Xmx${avail_mem}M" HaplotypeCaller $args -I $bam -O ${prefix}.haplotypes.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        def prefix      = task.ext.prefix ?: "${meta.group}.GATK"
        """
        touch ${prefix}.haplotypes.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """
}


process SENTIEON_HAPLOTYPING {
    label 'process_medium_cpus'
    label 'process_medium_memory'
    label 'stage'
    tag "$meta.group"

    input:
        tuple val(group), val(meta), file(bam), file(bai)

    output:
        tuple val(group), val(meta), file("*.haplotypes.vcf"),  emit: haplotypes 
        path "versions.yml",                                    emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args        = task.ext.args   ?: ''
        def args2       = task.ext.args2  ?: ''
        def prefix      = task.ext.prefix ?: "${meta.group}.sentieon"
        """
        sentieon driver -t ${task.cpus} $args -i $bam --algo Haplotyper $args2 ${prefix}.haplotypes.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.group}.sentieon"
        """
        touch ${prefix}.haplotypes.vcf ${prefix}.haplotypes.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        END_VERSIONS
        """
}
