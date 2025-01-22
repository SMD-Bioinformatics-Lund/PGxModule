#!/usr/bin/env nextflow

include { PHARMCAT_PREPROCESSING                } from '../../modules/local/pharmcat/main'
include { PHARMCAT_RUN                          } from '../../modules/local/pharmcat/main'

workflow PHARMCAT {

    take:
        filtered_haplotypes             // channel: [ tuple val(group), val(meta) file("haplotypes_filtered.vcf.gz") file("haplotypes_filtered.vcf.gz.tbi") ]

    main:
        ch_versions = Channel.empty()

        // Preprocess the pharmcat
        PHARMCAT_PREPROCESSING ( filtered_haplotypes )
        ch_versions = ch_versions.mix(PHARMCAT_PREPROCESSING.out.versions)

        // Run the pharmcat
        PHARMCAT_RUN ( PHARMCAT_PREPROCESSING.out.pharmcat_preprocessed_vcf )
        ch_versions = ch_versions.mix(PHARMCAT_RUN.out.versions)


    emit:
        pharmcat_preprocessed   = PHARMCAT_PREPROCESSING.out.pharmcat_preprocessed_vcf       // channel: [ tuple val(group), val(meta) file("pharmcat_preprocessed.vcf") ]
        pharmcat_report         = PHARMCAT_RUN.out.pharmcat_report                           // channel: [ tuple val(group), val(meta) file("pharmcat.html") ]
        pharmcat_pheno_json     = PHARMCAT_RUN.out.pharmcat_pheno_json                       // channel: [ tuple val(group), val(meta) file("pharmcat.phenotype.json") ]
        pharmcat_match_json     = PHARMCAT_RUN.out.pharmcat_match_json                       // channel: [ tuple val(group), val(meta) file("pharmcat.match.json") ]
        pharmcat_match_html     = PHARMCAT_RUN.out.pharmcat_match_html                       // channel: [ tuple val(group), val(meta) file("pharmcat.match.html") ]
        pharmcat_report_json    = PHARMCAT_RUN.out.pharmcat_report_json                      // channel: [ tuple val(group), val(meta) file("pharmcat.report.json") ]
        versions                = ch_versions                                                // channel: [ path(versions.yml) ]
}