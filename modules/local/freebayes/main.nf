process FREEBAYES {
    label "process_low"
    tag "$group"
    
    input:
        tuple val(group), val(meta), file(bam), file(bai), file(bqsr)

    output:
        tuple val(group), val(meta), file("*.freebayes.vcf"),   emit: freebayes_vcf
        path "versions.yml",                                    emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args        = task.ext.args                ?: ''
        def args2       = task.ext.args2               ?: ''
        def args3       = task.ext.args3               ?: ''
        def prefix      = task.ext.prefix ?: "${meta.group}.freebayes"
        """
        freebayes \\
        $args \\
        -F ${params.fb_var_freq_cutoff_up} \\
        $bam > ${prefix}.vcf.raw

        vcffilter $args2 ${prefix}.vcf.raw \\
        | vcffilter $args3 \\
        | vcfglxgt > ${prefix}.filt1.vcf

        filter_freebayes_unpaired.pl freebayes.filt1.vcf > ${prefix}.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            freebayes: \$(echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g' )
            vcffilter: \$(echo \$( vcffilter -h 2>&1) | grep 'vcflib' | sed 's/ filter.*\$//g' | sed 's/.* //g' )
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """
    stub:
        def prefix      = task.ext.prefix ?: "${meta.group}.freebayes"
        """
        echo tumor:$bams
        touch ${prefix}.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            freebayes: \$(echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g' )
            vcffilter: \$(echo \$( vcffilter -h 2>&1) | grep 'vcflib' | sed 's/ filter.*\$//g' | sed 's/.* //g' )
            perl: \$( echo \$(perl -v 2>&1) |sed 's/.*(v//; s/).*//')
        END_VERSIONS
        """
}