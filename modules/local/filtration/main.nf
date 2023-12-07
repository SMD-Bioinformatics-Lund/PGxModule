process VARIANT_FILTRATION {
    label 'process_single'
    label 'stage'
    tag "$group"

    input:
        tuple val(group), file(vcf)

    output:
        tuple val(group), file("${group}.haplotypes.filtered.annotated.vcf"), emit: haplotypes_filtered 
        
    script:
    """
    variant_filtration.py \
        --input_vcf=$vcf \
        --read_ratio=$params.read_ratio \
        --depth=$params.DP \
        --output_file=${group}.haplotypes.filtered.annotated.vcf
    """

    stub:
    """
    touch ${group}.haplotypes.filtered.annotated.vcf
    """

}