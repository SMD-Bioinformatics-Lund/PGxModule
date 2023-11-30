#!/usr/bin/env python3

import os
import sys
import argparse
from report_class import Report
import argparse


# DEFAULTS
TEMPLATE_FOLDER = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "resources", "templates"
)

# Create a Pharmacogenomics (PGx) Report
parser = argparse.ArgumentParser(description="Generate a Pharmacogenomics (PGx) report")

# Required Arguments
parser.add_argument("--group", type=str, required=True, help="Specify the sample group")
parser.add_argument(
    "--read_depth",
    type=int,
    default=100,
    help="Read depth threshold for target regions (default: 100)",
)
parser.add_argument(
    "--report_template",
    type=str,
    default=os.path.join(
        TEMPLATE_FOLDER,
        "report.html",
    ),
    help=f"Path to the report template file (default: {os.path.join(TEMPLATE_FOLDER,'report.html')})",
)
parser.add_argument(
    "--report_css",
    type=str,
    default=os.path.join(TEMPLATE_FOLDER, "report.css"),
    help=f"Path to the report css file (default: {os.path.join(TEMPLATE_FOLDER,'report.css')})",
)
parser.add_argument(
    "--detected_variants",
    type=str,
    required=True,
    help="Path to detected variants file",
)
parser.add_argument(
    "--missing_annotated_depth",
    type=str,
    required=True,
    help="Path to missing annotated depth file",
)
parser.add_argument(
    "--haplotype_definitions",
    type=str,
    required=True,
    help="Path to haplotype definitions file",
)
parser.add_argument(
    "--possible_diplotypes",
    type=str,
    required=True,
    help="Path to possible diplotypes file",
)
parser.add_argument(
    "--possible_interactions",
    type=str,
    required=True,
    help="Path to possible interactions file",
)
parser.add_argument(
    "--target_bed", type=str, required=True, help="Path to the target BED file"
)
parser.add_argument(
    "--padded_baits_depth",
    type=str,
    required=True,
    help="Path to padded baits depth file",
)
parser.add_argument(
    "--target_rsids", type=str, required=True, help="Path to target rsIDs file"
)
parser.add_argument(
    "--annotated_vcf", type=str, required=True, help="Path to annotated VCF file"
)
parser.add_argument(
    "--dbSNP_version", type=str, required=True, help="Specify the dbSNP version"
)
parser.add_argument(
    "--genome_version", type=str, required=True, help="Specify the Genome version"
)
parser.add_argument(
    "--output",
    type=str,
    required=True,
    help="Specify the output file path for the PGx report",
)


if __name__ == "__main__":
    args = parser.parse_args()

    # Example usage:
    report_instance = Report(
        group=args.group,
        read_depth=args.read_depth,
        detected_variants=args.detected_variants,
        missing_annotated_depth=args.missing_annotated_depth,
        haplotype_definitions=args.haplotype_definitions,
        possible_diplotypes=args.possible_diplotypes,
        possible_interactions=args.possible_interactions,
        target_bed=args.target_bed,
        padded_baits_depth=args.padded_baits_depth,
        target_rsids=args.target_rsids,
        annotated_vcf=args.annotated_vcf,
        dbSNP_version=args.dbSNP_version,
        genome_version=args.genome_version,
        output=args.output,
        report_template=args.report_template,
        report_css=args.report_css,
    )

    report_instance.create_report()

"""sumary_line
# example usage:

python /data/bnf/dev/ram/Pipelines/DSL2/PGxModule/bin/create_report.py --group 1DPYD-PGx-231113 --read_depth 100 --detected_variants 1DPYD-PGx-231113.detected_variants.tsv --missing_annotated_depth 1DPYD-PGx-231113.pgx_depth_at_missing_annotated.gdf --haplotype_definitions /data/bnf/dev/ram/Pipelines/DSL2/PGxModule/resources/haplotypes/haplotype_definitions.csv --possible_diplotypes 1DPYD-PGx-231113.possible_diplotypes.tsv --possible_interactions 1DPYD-PGx-231113.possible_interactions.tsv --target_bed exons_variants_pharmacogenomics_18_06_2019_ex_cyp2d6_hg38.bed --padded_baits_depth 1DPYD-PGx-231113.pgx.gdf --target_rsids target_rsid_hg38.bed --annotated_vcf 1DPYD-PGx-231113.haplotypes.filtered.annotated.vcf --dbSNP_version 151 --output example_output.html --genome_version hg38 --report_template /data/bnf/dev/ram/Pipelines/DSL2/PGxModule/resources/templates/report.html


singularity run --bind /fs1 --bind /data /data/bnf/dev/ram/Pipelines/DSL2/PGxModule/envs/jinja_report.sif python3 /data/bnf/dev/ram/Pipelines/DSL2/PGxModule/bin/create_report.py --group 1DPYD-PGx-231113 --read_depth 100 --detected_variants 1DPYD-PGx-231113.detected_variants.tsv --missing_annotated_depth 1DPYD-PGx-231113.pgx_depth_at_missing_annotated.gdf --haplotype_definitions /data/bnf/dev/ram/Pipelines/DSL2/PGxModule/resources/haplotypes/haplotype_definitions.csv --possible_diplotypes 1DPYD-PGx-231113.possible_diplotypes.tsv --possible_interactions 1DPYD-PGx-231113.possible_interactions.tsv --target_bed exons_variants_pharmacogenomics_18_06_2019_ex_cyp2d6_hg38.bed --padded_baits_depth 1DPYD-PGx-231113.pgx.gdf --target_rsids target_rsid_hg38.bed --annotated_vcf 1DPYD-PGx-231113.haplotypes.filtered.annotated.vcf --dbSNP_version 151 --output example_output.html --genome_version hg38 --report_template /data/bnf/dev/ram/Pipelines/DSL2/PGxModule/resources/templates/report.html --report_css /data/bnf/dev/ram/Pipelines/DSL2/PGxModule/resources/templates/report.css

"""
