#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PGX 							} from './workflows/pgx'

println(params.genome_file)
csv = file(params.csv)
println(csv)



workflow {

	PGX ()
}
