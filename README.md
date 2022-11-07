# PGxModule

This workflow is designed in nextflow DSL2 to generate sample-specific report with clinical guidelines dependent on variants detected in a Genomic Medicine Sweden sequencing panel. (Originally taken from https://github.com/JoelAAs/pgx_module)

It is designed to be coupled to an existing pipeline and need to be given the location of analysis ready bam files, and VCF with VEP annotations.

### Setup
Set paths in `nextflow.config` and `workflows/pgx_hg38.config` 
Add `bam_location` and `bam_location` to main.nf

If using singularity with nextflow please run script `envs/get_containers.sh` 

### This was tested using:
+ `N E X T F L O W version 22.10.1 build 5828`
+ _Hg38_ reference genome 
+ _dbsnp_ build 154