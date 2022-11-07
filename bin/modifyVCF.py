#!/usr/bin/env python3

import os
import sys
import argparse
from extractVCF import vcfRecord, get_rsIDs, reconstructAnnotation_withrsID


parser = argparse.ArgumentParser(
    description="Finds selected RSIDs form bed file in input VCF")

parser.add_argument("--input", type=str,required=True)
parser.add_argument("--addrsid", type=bool,default=False,choices=[True,False])
parser.add_argument("--output", type=str,required='--addrsid' in sys.argv)
args = parser.parse_args(sys.argv[1:])

inVCF = args.input
outVCF = args.output
addRSID = args.addrsid


infoFormat = "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|UNIPROT_ISOFORM|REFSEQ_MATCH|SOURCE|REFSEQ_OFFSET|GIVEN_REF|USED_REF|BAM_EDIT|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|AA_AF|EA_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|CADD_PHRED|CADD_RAW|LoFtool|gnomADg|gnomADg_AF_popmax|gnomADg_AF|gnomADg_popmax|COSMIC|COSMIC_CNT".split('|')

if addRSID:
    try:
        os.remove(outVCF)
    except:
        pass

    outVCF_handle = open(outVCF, 'a')

with open(inVCF) as vcf:
    lines = vcf.readlines()

for vcf.line in lines:
    if vcf.line.startswith("#"):
        if addRSID:
            outVCF_handle.write(f'{vcf.line}')
        continue
    else:
        vcftuple = vcfRecord(vcf.line.strip(), infoFormat)
        vcfDict = vcftuple[0]
        transcriptCount = vcftuple[1]

        if addRSID:
            rsIDs = get_rsIDs(vcfDict['INFO']['CSQ']['Existing_variation'])
            newVCF_entry = reconstructAnnotation_withrsID(vcfDict,transcriptCount,rsIDs)
            outVCF_handle.write(f'{newVCF_entry}\n')
        
if addRSID:
    outVCF_handle.close()