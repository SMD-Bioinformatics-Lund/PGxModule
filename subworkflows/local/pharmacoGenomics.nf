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
include { DUMPSOFTWAREVERSIONS                  } from '../../modules/local/dumpsoftwareversions/main'

workflow PHARMACO_GENOMICS {

    take: 
        bam_input
        pgx_targets_bed
        pgx_target_rsids
    
    main:
        PADDED_BED_INTERVALS ( pgx_targets_bed )
        GET_PADDED_BAITS ( pgx_targets_bed )
        ONTARGET_BAM ( bam_input.combine(PADDED_BED_INTERVALS.out.padded_bed_intervals) )
        HAPLOTYPER ( ONTARGET_BAM.out.bam_ontarget )
        BCFTOOLS_ANNOTATION ( HAPLOTYPER.out.haplotypes )
        VARIANT_FILTRATION ( BCFTOOLS_ANNOTATION.out.annotations )        
        DETECTED_VARIANTS ( VARIANT_FILTRATION.out.haplotypes_filtered.combine(pgx_target_rsids) )
        SAMPLE_TARGET_LIST ( DETECTED_VARIANTS.out.detected_tsv.combine(pgx_target_rsids) )
        DEPTH_OF_TARGETS ( ONTARGET_BAM.out.bam_ontarget.join(SAMPLE_TARGET_LIST.out.pgx_target_interval_list) )
        DEPTH_OF_BAITS ( ONTARGET_BAM.out.bam_ontarget.combine(GET_PADDED_BAITS.out.padded_baits_list) )
        APPEND_ID_TO_GDF ( DEPTH_OF_TARGETS.out.pgx_depth_at_missing.combine(pgx_target_rsids) )
        GET_CLIINICAL_GUIDELINES ( DETECTED_VARIANTS.out.detected_tsv )
        GET_INTERACTION_GUIDELINES ( GET_CLIINICAL_GUIDELINES.out.possible_diplotypes )
        GET_PGX_REPORT ( VARIANT_FILTRATION.out.haplotypes_filtered.join(DETECTED_VARIANTS.out.detected_tsv).join(APPEND_ID_TO_GDF.out.depth_at_missing_annotate_gdf).join( GET_CLIINICAL_GUIDELINES.out.possible_diplotypes).join(DEPTH_OF_BAITS.out.padded_baits_list).join(GET_INTERACTION_GUIDELINES.out.possible_interactions).combine(pgx_targets_bed).combine(pgx_target_rsids) )

        DUMPSOFTWAREVERSIONS (ONTARGET_BAM.out.versions.join(
                                            HAPLOTYPER.out.versions).join(
                                            BCFTOOLS_ANNOTATION.out.versions).join(
                                            DEPTH_OF_TARGETS.out.versions).join(
                                            DEPTH_OF_BAITS.out.versions)
                                            )
    emit:
        pgx_report = GET_PGX_REPORT.out.pgx_html
        versions   = DUMPSOFTWAREVERSIONS.out.versions_yml   // channel: [ path(versions.yml) ]

}