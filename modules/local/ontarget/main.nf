process ONTARGET_BAM {
    label 'process_low'
    label 'stage'
    tag "$meta.group"

    input:
        tuple val(group), val(meta), file(bam), file(bai) 
        path pgx_ontarget_padded_bed

    output:
        tuple val(group),  val(meta), file("*.dedup.ontarget.pgx.bam"), file("*.dedup.ontarget.pgx.bam.bai"),   emit: bam_ontarget
        path "versions.yml",                                                                                    emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args   ?: ''
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """
        samtools view $args -b $bam -L $pgx_ontarget_padded_bed > ${prefix}.dedup.ontarget.pgx.bam
        samtools index ${prefix}.dedup.ontarget.pgx.bam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(samtools --version 2>&1 | grep 'samtools' | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """

    stub:
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """
        touch ${prefix}.dedup.ontarget.pgx.bam ${prefix}.dedup.ontarget.pgx.bam.bai

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(samtools --version 2>&1 | grep 'samtools' | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
}


process GATK_HAPLOTYPING {
    label 'process_high'
    label 'stage'
    tag "$meta.group"

    input:
        tuple val(group), val(meta), file(bam), file(bai)

    output:
        tuple val(group), val(meta), file("*.haplotypes.vcf.gz"), file("*.haplotypes.vcf.gz.tbi"),  emit: haplotypes
        path "versions.yml",                                                                        emit: versions

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
        bgzip -c ${prefix}.haplotypes.vcf > ${prefix}.haplotypes.vcf.gz
        tabix ${prefix}.haplotypes.vcf.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
            bgzip: \$(bgzip --v | grep 'bgzip' | sed 's/.* //g')
            tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        def prefix      = task.ext.prefix ?: "${meta.group}.GATK"
        """
        touch ${prefix}.haplotypes.vcf.gz ${prefix}.haplotypes.vcf.gz.tbi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
            bgzip: \$(bgzip --v | grep 'bgzip' | sed 's/.* //g')
            tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        END_VERSIONS
        """
}


process SENTIEON_HAPLOTYPING {
    label 'process_high'
    label 'stage'
    tag "$meta.group"

    input:
        tuple val(group), val(meta), file(bam), file(bai)

    output:
        tuple val(group), val(meta), file("*.haplotypes.vcf.gz"), file("*.haplotypes.vcf.gz.tbi"),  emit: haplotypes
        path "versions.yml",                                                                        emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args        = task.ext.args   ?: ''
        def args2       = task.ext.args2  ?: ''
        def prefix      = task.ext.prefix ?: "${meta.group}.sentieon"
        """
        sentieon driver -t ${task.cpus} $args -i $bam --algo Haplotyper $args2 ${prefix}.haplotypes.vcf
        bgzip -c ${prefix}.haplotypes.vcf > ${prefix}.haplotypes.vcf.gz
        tabix ${prefix}.haplotypes.vcf.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
            bgzip: \$(bgzip --v | grep 'bgzip' | sed 's/.* //g')
            tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.group}.sentieon"
        """
        touch ${prefix}.haplotypes.vcf.gz ${prefix}.haplotypes.vcf.gz.tbi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
            bgzip: \$(bgzip --v | grep 'bgzip' | sed 's/.* //g')
            tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        END_VERSIONS
        """
}
