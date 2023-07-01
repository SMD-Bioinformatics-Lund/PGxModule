# PGxModule

This workflow is designed in nextflow DSL2 to generate sample-specific report with clinical guidelines dependent on variants detected in a Genomic Medicine Sweden sequencing panel. (Originally taken from <https://github.com/JoelAAs/pgx_module>)

It is designed to be coupled to an existing pipeline and need to be given the location of analysis ready bam files, and VCF with VEP annotations.

### Setup

Set paths in `nextflow.config` and `workflows/pgx_hg38.config`
Add `bam_location` and `bam_location` to main.nf

If using singularity with nextflow please run script `envs/get_containers.sh`

### Run Sentieon

To use Sentieon you have to:

1. Ensure that Sentieon executable is in path.
2. If necessary, store license information as a secret in Nextflow's local store.
3. Set environmental variable `NXF_ENABLE_SECRETS` to an appropriate value.

To elaborate more on item #2 in the list above, if you are running nf-core/raredisease with Sentieon on AWS or on a local machine and you do not want other users to know your license details, we recommend that you use [Nextflow's secrets feature](https://www.nextflow.io/docs/latest/secrets.html) to pass the that information. To do that run the command below after replacing LICENSE with the value corresponding to your license server (for example, 1.1.1.1:4000)

```
nextflow secrets set SENTIEON_LICENSE_BASE64 <LICENSE>  
```  

If you are using Nextflow secrets, you have to set the environment variable `NXF_ENABLE_SECRETS` to true. This will see to it that the pipeline can retrieve the secret from Nextflow's secrets store during the pipeline execution. Keep in mind that versions of Nextflow Version 22.09.2-edge and onwards have NXF_ENABLE_SECRETS to true by default. If you are not using secrets set `NXF_ENABLE_SECRETS` to false, but make sure that the environment variable [`SENTIEON_LICENSE`](`NXF_ENABLE_SECRETS`) is set to reflect the value of your license server on your machine.  

### This was tested using  

+ `N E X T F L O W version 22.10.1 build 5828`
+ _Hg38_ reference genome
+ _dbsnp_ build 154
