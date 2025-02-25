include { SENTIEON_HAPLOTYPING                      } from '../../modules/local/sentieon/main'
include { GATK_HAPLOTYPING                          } from '../../modules/local/gatk/main'
include { FREEBAYES                                 } from '../../modules/local/freebayes/main'
include { VT_DECOMPOSE_NORMALIZE as NORM_SENTIEON   } from '../../modules/local/vt/main'
include { VT_DECOMPOSE_NORMALIZE as NORM_GATK       } from '../../modules/local/vt/main'
include { VT_DECOMPOSE_NORMALIZE as NORM_FREEBAYES  } from '../../modules/local/vt/main'
include { AGGREGATE_VCFS                            } from '../../modules/local/aggregate_vcfs/main'



workflow VARIANT_CALLING {

    take:
        bam_ch

    main:
        ch_versions = Channel.empty()

        // Sentieon haplotyping
        SENTIEON_HAPLOTYPING ( bam_ch )
        ch_versions = ch_versions.mix(SENTIEON_HAPLOTYPING.out.versions)

        NORM_SENTIEON ( SENTIEON_HAPLOTYPING.out.haplotypes, "sentieon" )
        ch_versions = ch_versions.mix(NORM_SENTIEON.out.versions)
        
        // GATK haplotyping
        GATK_HAPLOTYPING ( bam_ch )
        ch_versions = ch_versions.mix(GATK_HAPLOTYPING.out.versions)

        NORM_GATK ( GATK_HAPLOTYPING.out.haplotypes, "gatk" )
        ch_versions = ch_versions.mix(NORM_GATK.out.versions)

        // Freebayes haplotyping
        FREEBAYES ( bam_ch )
        ch_versions = ch_versions.mix(FREEBAYES.out.versions)

        NORM_FREEBAYES ( FREEBAYES.out.freebayes_vcf, "freebayes" )
        ch_versions = ch_versions.mix(NORM_FREEBAYES.out.versions)

        // Aggregate all callers to one VCF
        AGGREGATE_VCFS ( NORM_SENTIEON.out.decomposed_normalized_vcfs.mix(NORM_GATK.out.decomposed_normalized_vcfs, NORM_FREEBAYES.out.decomposed_normalized_vcfs) )
        // ch_versions = ch_versions.mix(AGGREGATE_VCFS.out.versions)



    emit:
        sentieon_vcf        =   NORM_SENTIEON.out.decomposed_normalized_vcfs
        gatk_vcf            =   NORM_GATK.out.decomposed_normalized_vcfs
        freebayes_vcf       =   NORM_FREEBAYES.out.decomposed_normalized_vcfs
        aggregate_vcf       =   AGGREGATE_VCFS.out.vcf_agg
        aggregate_vcf_tbi   =   AGGREGATE_VCFS.out.vcf_agg_tbi
        versions            =   ch_versions

}