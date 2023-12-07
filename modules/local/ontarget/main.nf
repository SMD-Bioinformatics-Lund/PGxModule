process ONTARGET_BAM {
    label 'process_low'
    label 'stage'
    tag "$group"

    input:
        tuple val(group), val(id), val(type), file(bam), file(bai), file(pgx_ontarget_padded_bed)

    output:
        tuple val(group), file("${group}.dedup.ontarget.pgx.bam"), file("${group}.dedup.ontarget.pgx.bam.bai"), emit: bam_ontarget
        path "versions.yml",                                                                                    emit: versions

    script:
        """
        samtools view -h -b $bam -L $pgx_ontarget_padded_bed -M > ${group}.dedup.ontarget.pgx.bam
        samtools index ${group}.dedup.ontarget.pgx.bam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(samtools --version 2>&1 | grep 'samtools' | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """

    stub:
        """
        touch ${group}.dedup.ontarget.pgx.bam ${group}.dedup.ontarget.pgx.bam.bai

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(samtools --version 2>&1 | grep 'samtools' | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
}


process GATK_HAPLOTYPING {
    label 'process_high'
    label 'stage'
    tag "$group"

    input:
        tuple val(group), file(bam), file(bai)

    output:
        tuple val(group), file("${group}.GATK.haplotypes.vcf.gz"), file("${group}.GATK.haplotypes.vcf.gz.tbi"), emit: haplotypes
        path "versions.yml",                                                                                    emit: versions

    script:
        """
        gatk HaplotypeCaller -R $params.genome_file -I $bam -O ${group}.GATK.haplotypes.vcf
        bgzip -c ${group}.sentieon.haplotypes.vcf > ${group}.sentieon.haplotypes.vcf.gz
        tabix ${group}.sentieon.haplotypes.vcf.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
            bgzip: \$(bgzip --v | grep 'bgzip' | sed 's/.* //g')
            tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        """
        touch ${group}.sentieon.haplotypes.vcf.gz ${group}.sentieon.haplotypes.vcf.gz.tbi

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
    tag "$group"

    input:
        tuple val(group), file(bam), file(bai)

    output:
        tuple val(group), file("${group}.sentieon.haplotypes.vcf.gz"), file("${group}.sentieon.haplotypes.vcf.gz.tbi"), emit: haplotypes
        path "versions.yml",                                                                                            emit: versions
        
    script:
        """
        sentieon driver -t ${task.cpus} -r $params.genome_file -i $bam --algo Haplotyper --emit_mode confident ${group}.sentieon.haplotypes.vcf
        bgzip -c ${group}.sentieon.haplotypes.vcf > ${group}.sentieon.haplotypes.vcf.gz
        tabix ${group}.sentieon.haplotypes.vcf.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
            bgzip: \$(bgzip --v | grep 'bgzip' | sed 's/.* //g')
            tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        """
        touch ${group}.sentieon.haplotypes.vcf.gz ${group}.sentieon.haplotypes.vcf.gz.tbi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
            bgzip: \$(bgzip --v | grep 'bgzip' | sed 's/.* //g')
            tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        END_VERSIONS
        """
}
