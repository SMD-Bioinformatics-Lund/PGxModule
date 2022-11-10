#!/usr/bin/env nextflow

include { ONTARGET_BAM                  } from '../../modules/local/ontarget/main'
include { DETECTED_VARIANTS             } from '../../modules/local/detected_variants/main'
include { ANNOTATE_VARIANTS             } from '../../modules/local/annotation/main'
include { SAMPLE_TARGET_LIST            } from '../../modules/local/summary/main'
include { DEPTH_OF_TARGETS              } from '../../modules/local/summary/main'
include { GET_PADDED_BAITS              } from '../../modules/local/summary/main'
include { DEPTH_OF_BAITS                } from '../../modules/local/summary/main'
include { APPEND_ID_TO_GDF              } from '../../modules/local/summary/main'
include { PADDED_BED_INTERVALS          } from '../../modules/local/summary/main'
include { GET_CLIINICAL_GUIDELINES      } from '../../modules/local/pgx_report/main'
include { GET_INTERACTION_GUIDELINES    } from '../../modules/local/pgx_report/main'
include { GET_PGX_REPORT                } from '../../modules/local/pgx_report/main'

workflow PHARMACO_GENOMICS {
    take: 
        bam_input
        vcfs_input
        pgx_targets_bed
        pgx_target_rsids
    
    main:
        PADDED_BED_INTERVALS ( pgx_targets_bed )
        GET_PADDED_BAITS ( pgx_targets_bed )
        ONTARGET_BAM ( bam_input.combine(PADDED_BED_INTERVALS.out.padded_bed_intervals) )
        ANNOTATE_VARIANTS ( vcfs_input )
        DETECTED_VARIANTS ( ANNOTATE_VARIANTS.out.annotated_vcf.combine(pgx_target_rsids) )
        SAMPLE_TARGET_LIST ( DETECTED_VARIANTS.out.detected_csv.combine(pgx_target_rsids) )
        DEPTH_OF_TARGETS ( ONTARGET_BAM.out.bam_ontarget.join(SAMPLE_TARGET_LIST.out.pgx_target_interval_list) )
        DEPTH_OF_BAITS ( ONTARGET_BAM.out.bam_ontarget.combine(GET_PADDED_BAITS.out.padded_baits_list) )
        APPEND_ID_TO_GDF ( DEPTH_OF_TARGETS.out.pgx_depth_at_missing.combine(pgx_target_rsids) )
        GET_CLIINICAL_GUIDELINES ( DETECTED_VARIANTS.out.detected_csv )
        GET_INTERACTION_GUIDELINES ( GET_CLIINICAL_GUIDELINES.out.possible_diplotypes )
        GET_PGX_REPORT ( ANNOTATE_VARIANTS.out.annotated_vcf.join(DETECTED_VARIANTS.out.detected_csv).join(APPEND_ID_TO_GDF.out.depth_at_missing_annotate_gdf).join( GET_CLIINICAL_GUIDELINES.out.possible_diplotypes).join(DEPTH_OF_BAITS.out.padded_baits_list).join(GET_INTERACTION_GUIDELINES.out.possible_interactions).combine(pgx_targets_bed).combine(pgx_target_rsids) )
    emit:
        pgx_report = GET_PGX_REPORT.out.pgx_html

}