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

process ONTARGET_VCF {
    label 'process_low'
    label 'stage'
    tag "$meta.group"

    input:
        tuple val(group), val(meta), file(vcf), file(tbi) 
        path pgx_ontarget_padded_bed

    output:
        tuple val(group),  val(meta), file("*.ontarget.filtered.haplotypes.vcf.gz"), file("*.ontarget.filtered.haplotypes.vcf.gz.tbi"), emit: vcf_ontarget
        path "versions.yml",                                                                                                            emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args   ?: ''
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """
        tabix -h -R $pgx_ontarget_padded_bed $vcf | bgzip -c > ${prefix}.ontarget.filtered.haplotypes.vcf.gz
        tabix -p vcf ${prefix}.ontarget.filtered.haplotypes.vcf.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            tabix: \$(tabix --help 2>&1 | grep 'Version' | sed 's/^.*Version: //'')
        END_VERSIONS
        """

    stub:
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """
        touch ${prefix}.ontarget.filtered.haplotypes.vcf.gz ${prefix}.ontarget.filtered.haplotypes.vcf.gz.tbi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            tabix: \$(tabix --help 2>&1 | grep 'Version' | sed 's/^.*Version: //'')
        END_VERSIONS
        """
}


