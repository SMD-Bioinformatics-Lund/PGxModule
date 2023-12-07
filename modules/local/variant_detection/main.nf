process DETECTED_VARIANTS {
    label 'process_single'
    label 'stage'
    tag "$group"

    input:
        tuple val(group), file(vcf), file(pgx_ontarget_rsids_bed)
        
    output:
        tuple val(group), file("${group}.detected_variants.tsv"), emit: detected_tsv

    script:
    """
    get_target_variants.py \
        --target_bed $pgx_ontarget_rsids_bed \
        --vcf $vcf \
        --output ${group}.detected_variants.tsv \
        --addchr $params.addchr
    """

    stub:
    """
    touch ${group}.detected_variants.tsv
    """
}