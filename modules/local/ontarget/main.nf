process ONTARGET_BAM {
    label 'process_low'
    label 'stage'
    tag "$group"

    input:
        tuple val(group), val(id), val(type), file(bam), file(bai), file(pgx_ontarget_padded_bed)

    output:
        tuple val(group), file("${group}.dedup.ontarget.pgx.bam"), file("${group}.dedup.ontarget.pgx.bam.bai"), emit: bam_ontarget
        tuple val(group), file("${group}.${task.process.split(':').last()}.versions.yaml"), emit: versions

    script:
    def processName = task.process.toString().split(':').last()
    """
    samtools view -h -b $bam -L $pgx_ontarget_padded_bed -M > ${group}.dedup.ontarget.pgx.bam
    samtools index ${group}.dedup.ontarget.pgx.bam
    
    {
        echo -e "${processName}:"
        echo -e "\tSAMtools:"
        echo -e "\t\tversion: \$(samtools --version 2>&1 | grep 'samtools' | sed 's/^.*samtools //; s/Using.*\$//')"
        echo -e "\t\tcontainer: ${task.container}"
    } > "${group}.${processName}.versions.yaml"
    """

    stub:
    def processName = task.process.toString().split(':').last()
    """
    touch ${group}.dedup.ontarget.pgx.bam ${group}.dedup.ontarget.pgx.bam.bai

    {
        echo -e "${processName}:"
        echo -e "\tSAMtools:"
        echo -e "\t\tversion: \$(samtools --version 2>&1 | grep 'samtools' | sed 's/^.*samtools //; s/Using.*\$//')"
        echo -e "\t\tcontainer: ${task.container}"
    } > "${group}.${processName}.versions.yaml"
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
        tuple val(group), file("${group}.${task.process.split(':').last()}.versions.yaml"), emit: versions

    script:
    def processName = task.process.toString().split(':').last()
    """
    gatk HaplotypeCaller -R $params.genome_file -I $bam -O ${group}.GATK.haplotypes.vcf
    bgzip -c ${group}.sentieon.haplotypes.vcf > ${group}.sentieon.haplotypes.vcf.gz
    tabix ${group}.sentieon.haplotypes.vcf.gz
    {
        echo -e "${processName}:"
        echo -e "\tGATK HaplotypeCaller:"
        echo -e "\t\tversion: \$(gatk --version 2>&1 | grep 'The Genome Analysis Toolkit (GATK)' | sed -e 's/The Genome Analysis Toolkit (GATK) //g')"
        echo -e "\t\tcontainer: ${task.container}"
        echo -e "\tbgzip:"
        echo -e "\t\tversion: \$(bgzip --version 2>&1 | grep 'bgzip' | sed 's/.* //g')"
        echo -e "\t\tcontainer: ${task.container}"
        echo -e "\ttabix:"
        echo -e "\t\tversion: \$(tabix --version 2>&1 | grep 'tabix' | sed 's/.* //g')"
        echo -e "\t\tcontainer: ${task.container}"
    } > "${group}.${processName}.versions.yaml"
    """

    stub:
    def processName = task.process.toString().split(':').last()
    """
    touch ${group}.sentieon.haplotypes.vcf.gz ${group}.sentieon.haplotypes.vcf.gz.tbi
    {
        echo -e "${processName}:"
        echo -e "\tGATK HaplotypeCaller:"
        echo -e "\t\tversion: \$(gatk --version 2>&1 | grep 'The Genome Analysis Toolkit (GATK)' | sed -e 's/The Genome Analysis Toolkit (GATK) //g')"
        echo -e "\t\tcontainer: ${task.container}"
        echo -e "\tbgzip:"
        echo -e "\t\tversion: \$(bgzip --version 2>&1 | grep 'bgzip' | sed 's/.* //g')"
        echo -e "\t\tcontainer: ${task.container}"
        echo -e "\ttabix:"
        echo -e "\t\tversion: \$(tabix --version 2>&1 | grep 'tabix' | sed 's/.* //g')"
        echo -e "\t\tcontainer: ${task.container}"
    } > "${group}.${processName}.versions.yaml"
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
        tuple val(group), file("${group}.${task.process.split(':').last()}.versions.yaml"), emit: versions
        
    script:
    def processName = task.process.toString().split(':').last()
    """
    sentieon driver -t ${task.cpus} -r $params.genome_file -i $bam --algo Haplotyper --emit_mode confident ${group}.sentieon.haplotypes.vcf
    bgzip -c ${group}.sentieon.haplotypes.vcf > ${group}.sentieon.haplotypes.vcf.gz
    tabix ${group}.sentieon.haplotypes.vcf.gz

    {
        echo -e "${processName}:"
        echo -e "\tSentieon Haplotyper:"
        echo -e "\t\tversion: \$(sentieon driver --version 2>&1 | sed -e 's/sentieon-genomics-//g')"
        echo -e "\t\tcontainer: ${task.container}"
        echo -e "\tbgzip:"
        echo -e "\t\tversion: \$(bgzip --version 2>&1 | grep 'bgzip' | sed 's/.* //g')"
        echo -e "\t\tcontainer: ${task.container}"
        echo -e "\ttabix:"
        echo -e "\t\tversion: \$(tabix --version 2>&1 | grep 'tabix' | sed 's/.* //g')"
        echo -e "\t\tcontainer: ${task.container}"
    } > "${group}.${processName}.versions.yaml"

    """

    stub:
    def processName = task.process.toString().split(':').last()
    """
    touch ${group}.sentieon.haplotypes.vcf.gz ${group}.sentieon.haplotypes.vcf.gz.tbi
    {
        echo -e "${processName}:"
        echo -e "\tSentieon Haplotyper:"
        echo -e "\t\tversion: \$(sentieon driver --version 2>&1 | sed -e 's/sentieon-genomics-//g')"
        echo -e "\t\tcontainer: ${task.container}"
        echo -e "\tbgzip:"
        echo -e "\t\tversion: \$(bgzip --version 2>&1 | grep 'bgzip' | sed 's/.* //g')"
        echo -e "\t\tcontainer: ${task.container}"
        echo -e "\ttabix:"
        echo -e "\t\tversion: \$(tabix --version 2>&1 | grep 'tabix' | sed 's/.* //g')"
        echo -e "\t\tcontainer: ${task.container}"
    } > "${group}.${processName}.versions.yaml"
    """

}
