#!/usr/bin/bash

set -x

#Create variables for report
found_variants=$1
missed_variants=$2
haplotype_definitions=$3
clinical_guidelines=$4
interaction_guidelines=$5
sample=$6
target_regions_bed=$7
data_location=$8
depth_file=$9
target_rsid=${10}
workdir=${11}
rmd_html=${12}
output_file_html=${workdir}/${13}
annotated_vcf=${14}
dbSNP=${15}
assembly_genome=${16}


tmp_dir=$workdir/report_tmp
mkdir -p $tmp_dir
cp $workdir/* $tmp_dir

R -e "found_variants='$found_variants';missed_variants='$missed_variants';haplotype_definitions='$haplotype_definitions';clinical_guidelines='$clinical_guidelines';interaction_guidelines='$interaction_guidelines';sample='$sample';target_regions_bed='$target_regions_bed';data_location='$data_location';depth_file='$depth_file';target_rsid='$target_rsid';dbSNP='$dbSNP';assembly_genome='$assembly_genome';rmarkdown::render('$rmd_html',output_file='$output_file_html',intermediates_dir='$tmp_dir', knit_root_dir='$tmp_dir', clean=F)"