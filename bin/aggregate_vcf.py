import argparse
import gzip
import os
import sys
from collections import defaultdict
from typing import List

SUPPORTED_CALLERS: set[str] = {'freebayes', 'mutect2', 'tnscope', 'vardict', 'pindel', 'samtools', 'gatk-haplotyper', 'sentieon-haplotyper'}


def parse_args():
    parser = argparse.ArgumentParser(description="Aggregate multiple VCFs into one.")
    parser.add_argument('--vcfs', required=True, help='Comma separated list of VCF files to aggregate')
    parser.add_argument('--tumor-id', help='Tumor sample ID')
    parser.add_argument('--normal-id', help='Normal sample ID')
    parser.add_argument('--fluffify-pindel', action='store_true', help='Modify Pindel REF/ALT fields')
    parser.add_argument('--sample-order', help='Comma separated list of sample order')
    return parser.parse_args()


def read_vcf(file):
    metadata = []
    variants = []
    with (gzip.open(file, 'rt') if file.endswith('.gz') else open(file)) as f:
        for line in f:
            if line.startswith('##'):
                metadata.append(line.strip())
            elif line.startswith('#CHROM'):
                header = line.strip().split('\t')
            else:
                variants.append(line.strip().split('\t'))
    return metadata, header, variants


def which_variantcaller(metadata):
    source = None
    for line in metadata:
        if line.startswith('##source='):
            source = line.split('##source=')[1].split(' ')[0].split('_')[0].lower()
            break
        elif line.startswith('##GATKCommandLine'):
            source = "gatk-haplotyper"
            break
        elif line.startswith('##SentieonCommandLine.Haplotyper'):
            source = "sentieon-haplotyper"
            break

    for caller in SUPPORTED_CALLERS:
        if source and caller in source:
            return caller
    return "unknown"


def summarize_filters(filters):
    non_pass = [f for f in filters if f not in {"PASS", "."}]
    return "PASS" if not non_pass else ";".join(non_pass)


def is_weird_freebayes(var):
    if isinstance(var, list) and len(var) > 8:
        format_fields = var[8].split(':')
        if "AD" in format_fields:
            ad_index = format_fields.index("AD")
            sample_data = var[9].split(':') if len(var) > 9 else []
            return len(sample_data) <= ad_index or not sample_data[ad_index]
    return True


def aggregate_headers(headers):
    return headers[0]


def add_info_field(vcf_str, data):
    parts = vcf_str.split('\t')
    parts[7] += f";{data}"
    return '\t'.join(parts)


def excl_prefix(prefix, name):
    return name if prefix == "DO_NOTHING" else name.lstrip(prefix)


def add_info(var, key, val):
    var['info'][key] = val


def add_gt(var, sample, key, val):
    var['format'] = var.get('format', [])
    if key not in var['format']:
        var['format'].append(key)
    for gt in var['samples']:
        if gt['_sample_id'] == sample:
            gt[key] = val


def fix_gt(var, caller):
    for sample in var['samples']:
        gt_data = sample.split(':')
        ref_dp, alt_dp = 0, 0
        if caller in {'mutect2', 'tnscope', 'vardict', 'pindel'}:
            if len(gt_data) > 1 and ',' in gt_data[1]:
                ref_dp, alt_dp = map(int, gt_data[1].split(','))
            af = alt_dp / (alt_dp + ref_dp) if alt_dp + ref_dp > 0 else 0
            var['info']['VAF'] = str(af)
            var['info']['VD'] = str(alt_dp)
            var['info']['DP'] = str(alt_dp + ref_dp)


def aggregate_vcfs(vcf_files):
    aggregated = {}
    filters = defaultdict(set)
    all_filters = set()
    headers = []

    for vcf_file in vcf_files:
        metadata, header, variants = read_vcf(vcf_file)
        headers.append(metadata)
        caller = which_variantcaller(metadata)

        for var in variants:
            chrom, pos, var_id, ref, alt, qual, filt, info, fmt, *samples = var
            key = f"{chrom}_{pos}_{ref}_{alt}"
            
            if caller == "freebayes" and is_weird_freebayes(var):
                continue
            
            if key in aggregated:
                aggregated[key]['info']['variant_callers'] += f"|{caller}"
            else:
                aggregated[key] = {
                    'chrom': chrom, 'pos': pos, 'id': var_id, 'ref': ref, 'alt': alt,
                    'qual': qual, 'filter': filt, 'info': {'variant_callers': caller},
                    'format': fmt, 'samples': samples
                }
                
                if caller == "MELT":
                    add_info(aggregated[key], "custom", info)
                
                fix_gt(aggregated[key], caller)
            
            filters[key].update(filt.split(';'))
            all_filters.update(filt.split(';'))

    for key in aggregated:
        aggregated[key]['filter'] = summarize_filters(filters[key])
    
    return aggregated, headers, sorted(all_filters)


def print_header(filters, vcf_files):
    print("##fileformat=VCFv4.2")
    print(f"##origin={','.join(vcf_files)}")
    print('##INFO=<ID=variant_callers,Number=.,Type=String,Description="List of variant callers">')
    print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    print('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">')
    print('##FORMAT=<ID=VAF,Number=1,Type=Float,Description="ALT allele observation fraction">')
    print('##FORMAT=<ID=VD,Number=1,Type=Integer,Description="ALT allele observation count">')
    for filt in filters:
        print(f'##FILTER=<ID={filt},Description="{filt}">')
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLES")


def main():
    args = parse_args()
    vcf_files = args.vcfs.split(',')
    
    aggregated_vcfs, headers, filters = aggregate_vcfs(vcf_files)
    print_header(filters, vcf_files)
    
    for key, record in aggregated_vcfs.items():
        info_str = ';'.join(f"{k}={v}" for k, v in record['info'].items())
        sample_str = '\t'.join(record['samples'])
        print(f"{record['chrom']}\t{record['pos']}\t{record['id']}\t{record['ref']}\t{record['alt']}\t"
              f"{record['qual']}\t{record['filter']}\t{info_str}\t{record['format']}\t{sample_str}")


if __name__ == "__main__":
    main()
