process PHARMCAT_PREPROCESSING {
    label 'process_single'
    label 'stage'
    tag "$meta.group"

    input:
        tuple val(group), val(meta), file(vcf) 

    output:
        tuple val(group), val(meta), file("*.pharmcat.preprocessed.vcf"),    emit: pharmcat_preprocessed
        path "versions.yml",                                                emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args   ?: ''
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """
        /pharmcat/pharmcat_vcf_preprocessor.py -vcf $vcf --base-filename ${prefix}.pharmcat $args
        gunzip -c ${prefix}.pharmcat.preprocessed.vcf.bgz > ${prefix}.pharmcat.preprocessed.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            pharmcat_vcf_preprocessor: \$(/pharmcat/pharmcat_vcf_preprocessor.py -V 2>&1 | sed -e 's/PharmCAT VCF Preprocessor //g')
        END_VERSIONS
        """

    stub:
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """
        touch ${prefix}.pharmcat.preprocessed.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            pharmcat_vcf_preprocessor: \$(/pharmcat/pharmcat_vcf_preprocessor.py -V 2>&1 | sed -e 's/PharmCAT VCF Preprocessor //g')
        END_VERSIONS
        """

}


process PHARMCAT {
    label 'process_single'
    label 'stage'
    tag "$meta.group"

    input:
        tuple val(group), val(meta), file(vcf)

    output:
        tuple val(group), val(meta), file("*.phenotype.json"),  emit: pharmcat_pheno_json
        tuple val(group), val(meta), file("*.match.json"),      emit: pharmcat_macth_json
        tuple val(group), val(meta), file("*.report.html"),     emit: pharmcat_report
        path "versions.yml",                                    emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args   ?: ''
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """
        /pharmcat/pharmcat -vcf  $vcf --base-filename ${prefix}.pharmcat $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            PharmCAT: \$(/pharmcat/pharmcat --version 2>&1 | sed -e 's/PharmCAT //g')
        END_VERSIONS
        """

    stub:
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """
        touch ${prefix}.pharmcat.phenotype.json
        touch ${prefix}.match.json
        touch ${prefix}.report.html

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            PharmCAT: \$(/pharmcat/pharmcat --version 2>&1 | sed -e 's/PharmCAT //g')
        END_VERSIONS
        """

}