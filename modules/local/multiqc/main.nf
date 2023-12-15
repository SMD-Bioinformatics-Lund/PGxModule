process MULTIQC {
    label 'process_single'

    input:
        path multiqc_files, stageAs: "?/*"
        path(multiqc_config)
        path(extra_multiqc_config)
        path(multiqc_logo)
        val(sample_name)

    output:
        path "*multiqc_report.html", emit: report
        path "*_data"              , emit: data
        path "*_plots"             , optional:true, emit: plots
        path "versions.yml"        , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def config = multiqc_config ? "--config $multiqc_config" : ''
        def extra_config = extra_multiqc_config ? "--config $extra_multiqc_config" : ''
        def logo = multiqc_logo ? /--cl-config 'custom_logo: "${multiqc_logo}"'/ : ''
        """
        multiqc \\
            --force \\
            --filename ${sample_name}_multiqc_report.html \\
            $args \\
            $extra_config \\
            $config \\
            $logo . 

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
        END_VERSIONS
        """

    stub:
        """
        mkdir ${sample_name}_multiqc_data
        mkdir ${sample_name}_multiqc_plots
        touch ${sample_name}_multiqc_report.html

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
        END_VERSIONS
        """
}