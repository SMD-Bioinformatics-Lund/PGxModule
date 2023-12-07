#!/usr/bin/env nextflow

if (params.haplotyper == 'GATK') {

    include { GATK_HAPLOTYPING as HAPLOTYPER   } from '../../modules/local/ontarget/main'

} else if (params.haplotyper == 'SENTIEON') {

    include { SENTIEON_HAPLOTYPING as HAPLOTYPER  } from '../../modules/local/ontarget/main'
}

include { ONTARGET_BAM                          } from '../../modules/local/ontarget/main'
include { BCFTOOLS_ANNOTATION                   } from '../../modules/local/annotation/main'
include { VARIANT_FILTRATION                    } from '../../modules/local/filtration/main'
include { DETECTED_VARIANTS                     } from '../../modules/local/variant_detection/main'
include { SAMPLE_TARGET_LIST                    } from '../../modules/local/summary/main'
include { DEPTH_OF_TARGETS                      } from '../../modules/local/summary/main'
include { GET_PADDED_BAITS                      } from '../../modules/local/summary/main'
include { DEPTH_OF_BAITS                        } from '../../modules/local/summary/main'
include { APPEND_ID_TO_GDF                      } from '../../modules/local/summary/main'
include { PADDED_BED_INTERVALS                  } from '../../modules/local/summary/main'
include { GET_CLIINICAL_GUIDELINES              } from '../../modules/local/pgx_report/main'
include { GET_INTERACTION_GUIDELINES            } from '../../modules/local/pgx_report/main'
include { GET_PGX_REPORT                        } from '../../modules/local/pgx_report/main'

workflow PHARMACO_GENOMICS {

    take: 
        bam_input
        pgx_targets_bed
        pgx_target_rsids
    
    main:
        ch_versions = Channel.empty()

        PADDED_BED_INTERVALS ( pgx_targets_bed )
        ch_versions = ch_versions.mix(PADDED_BED_INTERVALS.out.versions)
        
        GET_PADDED_BAITS ( pgx_targets_bed )
        ch_versions = ch_versions.mix(GET_PADDED_BAITS.out.versions)

        ONTARGET_BAM ( bam_input.combine(PADDED_BED_INTERVALS.out.padded_bed_intervals) )
        ch_versions = ch_versions.mix(ONTARGET_BAM.out.versions)

        HAPLOTYPER ( ONTARGET_BAM.out.bam_ontarget )
        ch_versions = ch_versions.mix(HAPLOTYPER.out.versions)

        BCFTOOLS_ANNOTATION ( HAPLOTYPER.out.haplotypes )
        ch_versions = ch_versions.mix(BCFTOOLS_ANNOTATION.out.versions)

        VARIANT_FILTRATION ( BCFTOOLS_ANNOTATION.out.annotations )        
        ch_versions = ch_versions.mix(VARIANT_FILTRATION.out.versions)

        DETECTED_VARIANTS ( VARIANT_FILTRATION.out.haplotypes_filtered.combine(pgx_target_rsids) )
        ch_versions = ch_versions.mix(DETECTED_VARIANTS.out.versions)

        SAMPLE_TARGET_LIST ( DETECTED_VARIANTS.out.detected_tsv.combine(pgx_target_rsids) )
        ch_versions = ch_versions.mix(SAMPLE_TARGET_LIST.out.versions)

        DEPTH_OF_TARGETS ( ONTARGET_BAM.out.bam_ontarget.join(SAMPLE_TARGET_LIST.out.pgx_target_interval_list) )
        ch_versions = ch_versions.mix(DEPTH_OF_TARGETS.out.versions)

        DEPTH_OF_BAITS ( ONTARGET_BAM.out.bam_ontarget.combine(GET_PADDED_BAITS.out.padded_baits_list) )
        ch_versions = ch_versions.mix(DEPTH_OF_BAITS.out.versions)

        APPEND_ID_TO_GDF ( DEPTH_OF_TARGETS.out.pgx_depth_at_missing.combine(pgx_target_rsids) )
        ch_versions = ch_versions.mix(APPEND_ID_TO_GDF.out.versions)

        GET_CLIINICAL_GUIDELINES ( DETECTED_VARIANTS.out.detected_tsv )
        ch_versions = ch_versions.mix(GET_CLIINICAL_GUIDELINES.out.versions)

        GET_INTERACTION_GUIDELINES ( GET_CLIINICAL_GUIDELINES.out.possible_diplotypes )
        ch_versions = ch_versions.mix(GET_INTERACTION_GUIDELINES.out.versions)

        GET_PGX_REPORT ( VARIANT_FILTRATION.out.haplotypes_filtered.join(DETECTED_VARIANTS.out.detected_tsv).join(APPEND_ID_TO_GDF.out.depth_at_missing_annotate_gdf).join( GET_CLIINICAL_GUIDELINES.out.possible_diplotypes).join
        (DEPTH_OF_BAITS.out.padded_baits_list).join(GET_INTERACTION_GUIDELINES.out.possible_interactions).combine(pgx_targets_bed).combine(pgx_target_rsids) )
        ch_versions = ch_versions.mix(GET_PGX_REPORT.out.versions)

    emit:
        pgx_report = GET_PGX_REPORT.out.pgx_html
        versions   = ch_versions   // channel: [ path(versions.yml) ]

}