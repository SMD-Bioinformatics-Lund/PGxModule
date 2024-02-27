#!/usr/bin/python
import sys
import os
import shutil
import glob
import pandas as pd


def main():

    # remove the formatted directory if it exists
    if os.path.exists("formatted"):
        shutil.rmtree("formatted")

    os.makedirs("formatted")

    # dictionary of gene names and their corresponding chromosome number
    Gene_chromosome_number = {
        "1": ["NRAS", "CACNA1S", "DPYD", "F5", "MTHFR"],
        "2": ["UGT1A1", "ALK"],
        "3": [],
        "4": ["ABCG2", "KIT"],
        "5": ["ADRB2"],
        "6": ["HLA-B", "TPMT"],
        "7": ["ABCB1", "CFTR", "CYP3A4", "CYP3A5", "BRAF", "EGFR"],
        "8": ["NAT2"],
        "9": ["ABL1"],
        "10": ["ADRB1", "CYP2C19", "CYP2C8", "CYP2C9"],
        "11": ["DRD2", "GSTP1"],
        "12": ["SLCO1B1", "KRAS"],
        "13": ["NUDT15"],
        "14": [],
        "15": [],
        "16": ["VKORC1"],
        "17": ["ACE", "ERBB2"],
        "18": ["TYMS"],
        "19": ["CYP2A6", "CYP2B6", "CYP4F2", "IFNL3", "RYR1"],
        "20": [],
        "21": ["SLC19A1"],
        "22": ["COMT", "CYP2D6", "BCR"],
        "X": ["G6PD"],
        "Y": [],
        "MT": ["MT-RNR1"],
    }

    def get_chromosome(gene):
        for chromosome, genes in Gene_chromosome_number.items():
            if gene in genes:
                return chromosome
        return None

    # load all the csv files
    csv_files = glob.glob("raw/*.csv")

    for csv in csv_files:
        # read the csv file
        df = pd.read_csv(csv)
        df = df.fillna("")

        lines = df.values.tolist()

        # get the gene name
        gene_name = csv.split(".")[0]
        chrom = get_chromosome(gene_name)
        formatted_lines = [
            "No.,Exon / Intron,Chrom,Start,End,Start Phase,End Phase,Length,Sequence",
        ]

        bed_lines = []

        # create a bed file for each exon
        for line in lines:
            # format the line and separate the columns with tabs
            exon_number = line[0]
            region_name = line[1]
            start = (
                str(int(line[2])).replace(",", "")
                if "." in str(line[2])
                else str(line[2]).replace(",", "")
            )
            end = str(line[3]).replace(",", "")
            start_phase = str(line[4]).replace(",", "")
            end_phase = str(line[5]).replace(",", "")
            length = (
                str(int(line[6])).replace(",", "")
                if "." in str(line[6])
                else str(line[6]).replace(",", "")
            )
            sequence = line[-1].split("}")[-1].replace('"', "").strip()
            formatted_lines.append(
                f"{exon_number},{region_name},{chrom},{start},{end},{start_phase},{end_phase},{length},{sequence}"
            )

            if start == "" or end == "":
                continue

            if region_name.startswith("ENS"):
                bed_entry_name = f"Exon_{exon_number}_{gene_name}".replace(
                    " ", "_"
                ).replace("__", "_")
            else:
                bed_entry_name = (
                    f"{region_name}_{exon_number}_{gene_name}".replace(" ", "_")
                    .replace("__", "_")
                    .replace("__", "_")
                )

            if int(start) > int(end):
                # print(gene_name)
                bed_lines.append(f"{chrom}\t{end}\t{start}\t{bed_entry_name}")
                # print(f"{chrom}\t{start}\t{end}\t{bed_entry_name}")
                # break

            else:
                bed_lines.append(f"{chrom}\t{start}\t{end}\t{bed_entry_name}")
                # print(gene_name)
                # print(f"{chrom}\t{end}\t{start}\t{bed_entry_name}")
                # break

        # write the formatted lines to a new file
        with open(f"formatted/{gene_name}.csv", "w") as f:
            f.write("\n".join(formatted_lines))

        # write the bed lines to a new file
        with open(f"formatted/{gene_name}.bed", "w") as f:
            f.write("\n".join(bed_lines))

        # combine all the bed files into one
        with open(f"formatted/35_Tier1_9_cancer_PGx_HH38.bed", "a") as f:
            f.write("\n".join(bed_lines) + "\n")

    # Sort all genes bed file
    os.system(
        "sort -k1,1n -k2,2n -k3,3n formatted/35_Tier1_9_cancer_PGx_HH38.bed > formatted/35_Tier1_9_cancer_PGx_HH38_sorted.bed"
    )

    os.system(
        "grep -v '^X\|^Y\|^MT' formatted/35_Tier1_9_cancer_PGx_HH38_sorted.bed > formatted/35_Tier1_9_cancer_PGx_HH38_sorted_no_X_Y_MT.bed"
    )

    os.system(
        "grep '^X' formatted/35_Tier1_9_cancer_PGx_HH38_sorted.bed > formatted/35_Tier1_9_cancer_PGx_HH38_sorted_X.bed"
    )

    os.system(
        "grep '^Y' formatted/35_Tier1_9_cancer_PGx_HH38_sorted.bed > formatted/35_Tier1_9_cancer_PGx_HH38_sorted_Y.bed"
    )

    os.system(
        "grep '^MT' formatted/35_Tier1_9_cancer_PGx_HH38_sorted.bed > formatted/35_Tier1_9_cancer_PGx_HH38_sorted_MT.bed"
    )

    os.system(
        "cat formatted/35_Tier1_9_cancer_PGx_HH38_sorted_no_X_Y_MT.bed formatted/35_Tier1_9_cancer_PGx_HH38_sorted_X.bed formatted/35_Tier1_9_cancer_PGx_HH38_sorted_Y.bed formatted/35_Tier1_9_cancer_PGx_HH38_sorted_MT.bed > formatted/35_Tier1_9_cancer_PGx_HH38_sorted.bed"
    )

    os.system(
        "rm formatted/35_Tier1_9_cancer_PGx_HH38_sorted_no_X_Y_MT.bed formatted/35_Tier1_9_cancer_PGx_HH38_sorted_X.bed formatted/35_Tier1_9_cancer_PGx_HH38_sorted_Y.bed formatted/35_Tier1_9_cancer_PGx_HH38_sorted_MT.bed"
    )


if __name__ == "__main__":
    main()
