import argparse
import os
import pandas as pd
import numpy as np
import logging
from datetime import datetime
import matplotlib.pyplot as plt

# Read input CNV file
def read_cnv_file(cnv_file):
    """
    Reads a CNV file into a Pandas DataFrame.

    Parameters:
    cnv_file (str): Path to the CNV file.

    Returns:
    pd.DataFrame: DataFrame containing the CNV data.
    """
    return pd.read_csv(cnv_file)

def gene_starts_and_ends(bed_file):
    """
    Reads a BED file and extracts the start and end positions of unique genes.

    Parameters:
    bed_file (str): Path to the BED file.

    Returns:
    pd.DataFrame: DataFrame with columns ["chr", "start", "end", "gene"].
    """
    # Read the BED file (assuming tab-separated values with columns: chr, start, end, gene)
    df_bed = pd.read_csv(bed_file, sep="\t", names=["chr", "start", "end", "gene"])

    gene_limits = []
    this_gene = None
    this_start = None

    for index, row in df_bed.iterrows():
        gene_name_parts = str(row["gene"]).split("_")

        # Handling gene naming convention
        if "/" in str(row["gene"]) and "ENS" not in str(row["gene"]):
            this_line_gene_name = f"{gene_name_parts[0]}_{gene_name_parts[1]}"
        else:
            this_line_gene_name = gene_name_parts[0]

        # If this is the first gene encountered
        if this_gene is None:
            this_gene = this_line_gene_name
            this_start = row["start"] + 1  # Convert 0-based start to 1-based

        # If a new gene is encountered, store previous gene details
        if this_line_gene_name != this_gene:
            gene_limits.append([row["chr"], this_start, df_bed.loc[index - 1, "end"], this_gene])

            # Reset values for the new gene
            this_gene = this_line_gene_name
            this_start = row["start"] + 1

        # If it's the last row, store the final gene details
        if index == len(df_bed) - 1:
            gene_limits.append([row["chr"], this_start, row["end"], this_gene])

    return pd.DataFrame(gene_limits, columns=["chr", "start", "end", "gene"])

def plot_complete_panel(data_frame, gene_lims, sample_name, z_threshold, cn_threshold, output_dir):
    """
    Generates plots for CNV analysis including Z-score, Copy Number, and Read Depth.
    
    Parameters:
    - data_frame: Pandas DataFrame containing CNV data (columns: Chr, Position, Depth, prop, mean, std, expDepth, zscore, copynumber)
    - gene_lims: DataFrame containing chromosome start and end points for genes
    - sample_name: Name of the sample for output file naming
    - z_threshold: Z-score threshold
    - cn_threshold: Copy number threshold
    - output_dir: Output directory to save the plots
    """
    
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Get unique chromosomes and first row indices where they appear
    panel_chrs = data_frame["Chr"].unique()
    panel_chrs_row_no = {chr_: data_frame.index[data_frame["Chr"] == chr_][0] for chr_ in panel_chrs}

    # Extract series for plotting
    y_series_z = data_frame["zscore"]
    y_series_cn = data_frame["copynumber"]
    y_series_depth = data_frame["Depth"]
    y_series_expDepth = data_frame["expDepth"]
    y_series_stdev = np.abs((data_frame["Depth"] / data_frame["prop"]) * data_frame["std"])

    # Identify points above thresholds (for red markers)
    threshold_series = (np.abs(y_series_z) > z_threshold) & (np.abs(y_series_cn - 2) > cn_threshold)

    # Create figure with 3 subplots
    fig, axes = plt.subplots(3, 1, figsize=(15, 10), sharex=True)
    
    # Read Depth Plot
    axes[0].errorbar(range(len(y_series_depth)), y_series_expDepth, yerr=y_series_stdev, fmt='o', color="gray", label="Expected Depth")
    axes[0].scatter(range(len(y_series_depth)), y_series_depth, c=np.where(threshold_series, "red", "black"), s=5, label="Observed Depth")
    axes[0].axhline(y=20, linestyle="dashed", color="gray")  # Threshold line
    axes[0].set_ylabel("Read Depth")
    axes[0].legend()

    # Copy Number Plot
    axes[1].scatter(range(len(y_series_cn)), y_series_cn, c=np.where(threshold_series, "red", "black"), s=5, label="Copy Number")
    axes[1].axhline(y=2, linestyle="dashed", color="gray")
    axes[1].axhline(y=2 + cn_threshold, linestyle="dotted", color="gray")
    axes[1].axhline(y=2 - cn_threshold, linestyle="dotted", color="gray")
    axes[1].set_ylabel("Copy Number")

    # Z-score Plot
    axes[2].scatter(range(len(y_series_z)), y_series_z, c=np.where(threshold_series, "red", "black"), s=5, label="Z-score")
    axes[2].axhline(y=0, linestyle="dashed", color="gray")
    axes[2].axhline(y=z_threshold, linestyle="dotted", color="gray")
    axes[2].axhline(y=-z_threshold, linestyle="dotted", color="gray")
    axes[2].set_ylabel("Z-score")
    axes[2].set_xlabel("Positions covered by panel")

    # Add vertical lines for chromosomes
    for chr_, row_no in panel_chrs_row_no.items():
        for ax in axes:
            ax.axvline(x=row_no, linestyle="dotted", color="gray")
            ax.text(row_no, ax.get_ylim()[1] * 0.9, chr_, rotation=90, verticalalignment="top")

    # Title and save plot
    plt.suptitle(f"{sample_name} - Complete Panel CNV Analysis", fontsize=14)
    plt.tight_layout()
    output_path = os.path.join(output_dir, f"{sample_name}_cnv_complete_panel.png")
    plt.savefig(output_path)
    plt.close()

    print(f"Plot saved at {output_path}")

def cnv_report(data_frame, gene_lims, sample_name, z_threshold, cn_threshold, prop_threshold, bed_file, flag_bed_file, output_dir):
    """
    Generates a CNV report based on z-score and copy number thresholds.

    Parameters:
    - data_frame: Pandas DataFrame with columns: Chr, Position, Depth, prop, mean, std, expDepth, zscore, copynumber
    - gene_lims: Pandas DataFrame containing chromosome start and end positions for genes
    - sample_name: Sample name prefix for output files
    - z_threshold: Threshold for the z-score
    - cn_threshold: Threshold for the copy number deviation from 2
    - prop_threshold: Proportion threshold for calling CNVs
    - bed_file: Path to BED file containing gene annotations
    - flag_bed_file: Path to flag BED file with problematic genomic regions
    - output_dir: Directory to save the CNV report

    Returns:
    - DataFrame of CNV regions
    """

    print(f"Creating CNV report for {sample_name}...")

    report_file = os.path.join(output_dir, f"{sample_name}_cnv_report_all.txt")
    filtered_report_file = os.path.join(output_dir, f"{sample_name}_cnv_report_filtered.txt")

    # Read BED and flag BED files
    bed_df = pd.read_csv(bed_file, sep="\t", names=["chr", "start_pos", "end_pos", "gene"])
    flag_bed_df = pd.read_csv(flag_bed_file, sep="\t", names=["chr", "start_pos", "end_pos", "annot"])

    # Identify positions exceeding both CNV and Z-score thresholds
    threshold_series = (data_frame["zscore"].abs() > z_threshold) & ((data_frame["copynumber"] - 2).abs() > cn_threshold)

    cnv_records = []

    if threshold_series.any():
        # Process each gene region from BED file
        for _, gene_row in bed_df.iterrows():
            chr_region = gene_row["chr"]
            start = gene_row["start_pos"]
            end = gene_row["end_pos"]
            gene_name = gene_row["gene"]

            # Filter CNV data for the current gene region
            gene_cnv_df = data_frame[(data_frame["Chr"] == chr_region) & 
                                     (data_frame["Position"] >= start) & 
                                     (data_frame["Position"] <= end)]
            
            if not gene_cnv_df.empty:
                # Calculate mean CN and Z-score
                mean_cn = gene_cnv_df["copynumber"].mean()
                mean_z = gene_cnv_df["zscore"].mean()
                mean_depth = gene_cnv_df["Depth"].mean()
                exp_mean_depth = gene_cnv_df["expDepth"].mean()
                
                # Determine CNV type
                if mean_cn > 2:
                    cnv_type = "DUP"
                elif mean_cn < 2:
                    cnv_type = "DEL"
                else:
                    cnv_type = "UNDETERMINED"

                # Count positions above threshold
                reg_length = len(gene_cnv_df)
                reg_pos_above = threshold_series.loc[gene_cnv_df.index].sum()
                reg_percent_above = reg_pos_above / reg_length

                # Flagging CNV
                reg_flag = "PASS" if reg_percent_above >= prop_threshold else "FAIL"
                reg_comment = "FAIL - Low proportion above threshold. " if reg_flag == "FAIL" else ""

                # Check for low mean expected depth
                if exp_mean_depth < 50:
                    reg_flag = "FAIL"
                    reg_comment += "FAIL - Cohort mean depth too low. "

                # Check for overlap with problematic regions
                flag_overlaps = flag_bed_df[(flag_bed_df["chr"] == chr_region) &
                                            ((flag_bed_df["start_pos"] <= end) & 
                                             (flag_bed_df["end_pos"] >= start))]

                if not flag_overlaps.empty:
                    reg_comment += "WARNING - Overlaps with problematic region: "
                    for _, flag_row in flag_overlaps.iterrows():
                        reg_comment += f"{flag_row['chr']}:{flag_row['start_pos']}-{flag_row['end_pos']} {flag_row['annot']} "

                # Save CNV data
                cnv_records.append([chr_region, start, end, gene_name, cnv_type, mean_depth, exp_mean_depth, 
                                    mean_cn, mean_z, reg_length, reg_pos_above, reg_percent_above, reg_flag, reg_comment.strip()])

    # Convert to DataFrame
    cnv_df = pd.DataFrame(cnv_records, columns=[
        "chr", "start_pos", "end_pos", "gene", "cnv", "mean_depth", "exp_mean_depth", "mean_cn", 
        "mean_z-score", "reg.length", "pos.above", "%above", "flag", "comment"
    ])

    # Save full CNV report
    cnv_df.to_csv(report_file, sep="\t", index=False)
    print(f"Report saved: {report_file}")

    # Save only PASS CNVs
    cnv_df_filtered = cnv_df[cnv_df["flag"] == "PASS"]
    cnv_df_filtered.to_csv(filtered_report_file, sep="\t", index=False)
    print(f"Filtered report saved: {filtered_report_file}")

    return cnv_df

def plot_complete_panel_all_pos(data_frame, gene_lims, sample_name, chr_length_file, z_threshold, cn_threshold, output_dir):
    """
    Plots CNV data across all chromosome positions.

    Parameters:
    - data_frame: Pandas DataFrame containing columns: Chr, Position, Depth, prop, mean, std, expDepth, zscore, copynumber
    - gene_lims: Pandas DataFrame with gene start and end positions
    - sample_name: Name of the sample for plot title and output file naming
    - chr_length_file: Path to chromosome length file
    - z_threshold: Threshold for Z-score
    - cn_threshold: Threshold for Copy Number deviation from 2
    - output_dir: Directory to save the plots

    Saves:
    - A combined genome-wide CNV plot with Depth, CNV, and Z-score across all chromosomes.
    """

    # Load chromosome length file
    chr_lengths = pd.read_csv(chr_length_file, sep="\t", names=["chr", "length", "cumulative_length"])

    # Compute adjusted genome-wide positions
    adj_pos = data_frame.apply(lambda row: row["Position"] + 
                               chr_lengths.loc[chr_lengths["chr"] == row["Chr"], "cumulative_length"].values[0], axis=1)

    # Identify positions exceeding both CNV and Z-score thresholds
    threshold_series = (data_frame["zscore"].abs() > z_threshold) & ((data_frame["copynumber"] - 2).abs() > cn_threshold)

    # Prepare data for plotting
    x_series = adj_pos
    y_series_z = data_frame["zscore"]
    y_series_cn = data_frame["copynumber"]
    y_series_depth = data_frame["Depth"]
    y_series_expDepth = data_frame["expDepth"]
    y_series_stdev = np.abs((data_frame["Depth"] / data_frame["prop"]) * data_frame["std"])

    # Setup plot layout
    fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(12, 8), sharex=True)

    # Plot Read Depth
    axes[0].errorbar(x_series, y_series_expDepth, yerr=y_series_stdev, fmt='o', markersize=2, color='gray', label="Expected Depth")
    axes[0].scatter(x_series, y_series_depth, c=threshold_series.map({True: "red", False: "black"}), s=3)
    axes[0].axhline(y=20, linestyle='dashed', color='gray')
    axes[0].set_ylabel("Read Depth")
    axes[0].legend()

    # Plot Copy Number
    axes[1].scatter(x_series, y_series_cn, c=threshold_series.map({True: "red", False: "black"}), s=3)
    axes[1].axhline(y=2, linestyle='dashed', color='gray')
    axes[1].axhline(y=2 + cn_threshold, linestyle='dotted', color='gray')
    axes[1].axhline(y=2 - cn_threshold, linestyle='dotted', color='gray')
    axes[1].set_ylabel("Copy Number")
    axes[1].set_ylim(0, max(y_series_cn.max(), 4))

    # Plot Z-score
    axes[2].scatter(x_series, y_series_z, c=threshold_series.map({True: "red", False: "black"}), s=3)
    axes[2].axhline(y=0, linestyle='dashed', color='gray')
    axes[2].axhline(y=z_threshold, linestyle='dotted', color='gray')
    axes[2].axhline(y=-z_threshold, linestyle='dotted', color='gray')
    axes[2].set_ylabel("Z-score")
    axes[2].set_xlabel("Genomic Position")

    # Add chromosome separation lines
    for _, row in chr_lengths.iloc[1:].iterrows():
        for ax in axes:
            ax.axvline(x=row["cumulative_length"], linestyle='dotted', color='gray')

    # Add chromosome labels
    for i, row in chr_lengths.iloc[1:].iterrows():
        mid_chr_pos = row["cumulative_length"] - (row["length"] / 2)
        axes[-1].text(mid_chr_pos, z_threshold * 1.2, row["chr"], ha='center', fontsize=8)

    # Final plot settings
    fig.suptitle(f"{sample_name} - Complete Panel")
    plt.tight_layout(rect=[0, 0.03, 1, 0.97])

    # Save plot
    output_file = os.path.join(output_dir, f"{sample_name}_cnv_whole_genome.png")
    plt.savefig(output_file, dpi=300)
    print(f"Plot saved as {output_file}")
    plt.close()

def plot_chrs(data_frame, gene_lims, sample_name, chr_list, z_threshold, cn_threshold, output_dir):
    """
    Plots CNV data for each chromosome separately.

    Parameters:
    - data_frame: Pandas DataFrame containing columns: Chr, Position, Depth, prop, mean, std, expDepth, zscore, copynumber
    - gene_lims: Pandas DataFrame with gene start and end positions (chr, start, end, gene)
    - sample_name: Name of the sample for plot title and output file naming
    - chr_list: List of chromosomes to plot
    - z_threshold: Threshold for Z-score
    - cn_threshold: Threshold for Copy Number deviation from 2
    - output_dir: Directory to save the plots

    Saves:
    - PNG plots for each chromosome.
    """

    print("Plotting chromosomes...")

    # Iterate through each chromosome in the list
    for chrom in chr_list:
        print(f"Processing {chrom}...")

        # Filter data for the current chromosome
        this_chr_df = data_frame[data_frame["Chr"] == chrom].copy()
        this_chr_gene_lims = gene_lims[gene_lims["chr"] == chrom].copy()

        if this_chr_df.empty:
            print(f"No data for chromosome {chrom}. Skipping...")
            continue

        # Create vertical line positions for genes
        gene_starts = []
        gene_ends = []
        gene_labels = []

        for _, row in this_chr_gene_lims.iterrows():
            if row["start"] in this_chr_df["Position"].values and row["end"] in this_chr_df["Position"].values:
                gene_starts.append(row["start"])
                gene_ends.append(row["end"])
                gene_labels.append(row["gene"])

        # Identify CNV and Z-score threshold exceedances
        threshold_series = (this_chr_df["zscore"].abs() > z_threshold) & ((this_chr_df["copynumber"] - 2).abs() > cn_threshold)

        # Extract values for plotting
        x_series = this_chr_df["Position"]
        y_series_z = this_chr_df["zscore"]
        y_series_cn = this_chr_df["copynumber"]
        y_series_depth = this_chr_df["Depth"]
        y_series_expDepth = this_chr_df["expDepth"]
        y_series_stdev = np.abs((this_chr_df["Depth"] / this_chr_df["prop"]) * this_chr_df["std"])

        # Setup plot layout
        fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(12, 8), sharex=True)

        # Read Depth Plot
        axes[0].errorbar(x_series, y_series_expDepth, yerr=y_series_stdev, fmt='o', markersize=2, color='gray', label="Expected Depth")
        axes[0].scatter(x_series, y_series_depth, c=threshold_series.map({True: "red", False: "black"}), s=3)
        axes[0].axhline(y=20, linestyle='dashed', color='gray')
        axes[0].set_ylabel("Read Depth")
        axes[0].legend()

        # Copy Number Plot
        axes[1].scatter(x_series, y_series_cn, c=threshold_series.map({True: "red", False: "black"}), s=3)
        axes[1].axhline(y=2, linestyle='dashed', color='gray')
        axes[1].axhline(y=2 + cn_threshold, linestyle='dotted', color='gray')
        axes[1].axhline(y=2 - cn_threshold, linestyle='dotted', color='gray')
        axes[1].set_ylabel("Copy Number")
        axes[1].set_ylim(0, max(y_series_cn.max(), 4))

        # Z-score Plot
        axes[2].scatter(x_series, y_series_z, c=threshold_series.map({True: "red", False: "black"}), s=3)
        axes[2].axhline(y=0, linestyle='dashed', color='gray')
        axes[2].axhline(y=z_threshold, linestyle='dotted', color='gray')
        axes[2].axhline(y=-z_threshold, linestyle='dotted', color='gray')
        axes[2].set_ylabel("Z-score")
        axes[2].set_xlabel("Genomic Position")

        # Add vertical gene lines
        for ax in axes:
            for start, end in zip(gene_starts, gene_ends):
                ax.axvline(x=start, linestyle='solid', color='black')
                ax.axvline(x=end, linestyle='solid', color='black')

        # Add gene annotations
        for start, gene in zip(gene_starts, gene_labels):
            axes[2].text(start, z_threshold * 1.2, gene, ha='center', fontsize=8, rotation=45)

        # Final plot settings
        fig.suptitle(f"{sample_name} - Chromosome {chrom}")
        plt.tight_layout(rect=[0, 0.03, 1, 0.97])

        # Save plot
        output_file = os.path.join(output_dir, f"{sample_name}_cnv_{chrom}.png")
        plt.savefig(output_file, dpi=300)
        print(f"Plot saved as {output_file}")
        plt.close()

def plot_region(data_frame, gene_lims, sample_name, z_threshold, cn_threshold, chrom, start_pos, end_pos, bed_file, output_dir):
    """
    Plots CNV data for a specific genomic region.

    Parameters:
    - data_frame: Pandas DataFrame containing columns: Chr, Position, Depth, prop, mean, std, expDepth, zscore, copynumber
    - gene_lims: DataFrame with gene start and end positions (chr, start, end, gene)
    - sample_name: Sample name for labeling and file naming
    - z_threshold: Z-score threshold
    - cn_threshold: Copy Number threshold
    - chrom: Chromosome of the region
    - start_pos: Start position of the region
    - end_pos: End position of the region
    - bed_file: Path to BED file
    - output_dir: Directory to save plots

    Saves:
    - PNG plot of the specified genomic region
    """
    print(f"Creating plot for region {chrom}:{start_pos}-{end_pos}")

    chr = chrom.replace("chr","")   

    # Load BED file containing annotation regions
    bed_df = pd.read_csv(bed_file, sep='\t', names=["chr", "start_pos", "end_pos", "annot"])

    print(data_frame.head())

    # Filter CNV data for the given region
    this_reg_df = data_frame[
        (data_frame["Chr"] == chr) & 
        (data_frame["Position"] >= start_pos) & 
        (data_frame["Position"] <= end_pos)
    ].copy()


    if this_reg_df.empty:
        print(f"No data available for region {chrom}:{start_pos}-{end_pos}. Skipping...")
        return

    # Filter gene limits for the given region
    this_reg_gene_lims = gene_lims[
        (gene_lims["chr"] == chrom) & 
        (gene_lims["start"] < end_pos) & 
        (gene_lims["end"] > start_pos)
    ].copy()

    # Identify CNV and Z-score threshold exceedances
    threshold_series = (this_reg_df["zscore"].abs() > z_threshold) & ((this_reg_df["copynumber"] - 2).abs() > cn_threshold)

    # Extract values for plotting
    x_series = this_reg_df["Position"]
    y_series_z = this_reg_df["zscore"]
    y_series_cn = this_reg_df["copynumber"]
    y_series_depth = this_reg_df["Depth"]
    y_series_expDepth = this_reg_df["expDepth"]
    y_series_stdev = np.abs((this_reg_df["Depth"] / this_reg_df["prop"]) * this_reg_df["std"])

    # Setup plot layout
    fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(12, 8), sharex=True)

    # Read Depth Plot
    axes[0].errorbar(x_series, y_series_expDepth, yerr=y_series_stdev, fmt='o', markersize=2, color='gray', label="Expected Depth")
    axes[0].scatter(x_series, y_series_depth, c=threshold_series.map({True: "red", False: "black"}), s=3)
    axes[0].axhline(y=20, linestyle='dashed', color='gray')
    axes[0].set_ylabel("Read Depth")
    axes[0].legend()

    # Copy Number Plot
    axes[1].scatter(x_series, y_series_cn, c=threshold_series.map({True: "red", False: "black"}), s=3)
    axes[1].axhline(y=2, linestyle='dashed', color='gray')
    axes[1].axhline(y=2 + cn_threshold, linestyle='dotted', color='gray')
    axes[1].axhline(y=2 - cn_threshold, linestyle='dotted', color='gray')
    axes[1].set_ylabel("Copy Number")
    axes[1].set_ylim(0, max(y_series_cn.max(), 4))

    # Z-score Plot
    axes[2].scatter(x_series, y_series_z, c=threshold_series.map({True: "red", False: "black"}), s=3)
    axes[2].axhline(y=0, linestyle='dashed', color='gray')
    axes[2].axhline(y=z_threshold, linestyle='dotted', color='gray')
    axes[2].axhline(y=-z_threshold, linestyle='dotted', color='gray')
    axes[2].set_ylabel("Z-score")
    axes[2].set_xlabel("Genomic Position")

    # Final plot settings
    fig.suptitle(f"{sample_name} - {chrom}:{start_pos}-{end_pos}")
    plt.tight_layout(rect=[0, 0.03, 1, 0.97])

    # Save plot
    output_file = os.path.join(output_dir, f"{sample_name}_{chrom}_{start_pos}_{end_pos}.png")
    plt.savefig(output_file, dpi=300)
    print(f"Plot saved as {output_file}")
    plt.close()

def plot_cnv_regions(data_frame, gene_lims, sample_name, z_threshold, cn_threshold, cnv_regions, bed_file, output_dir):
    """
    Plots CNV regions by generating zoomed-in plots for significant CNVs.

    Parameters:
    - data_frame: Pandas DataFrame containing CNV data
    - gene_lims: DataFrame containing gene start/end positions
    - sample_name: Sample name for labeling
    - z_threshold: Z-score threshold
    - cn_threshold: Copy Number threshold
    - cnv_regions: DataFrame with detected CNVs (chr, start_pos, end_pos, gene, etc.)
    - bed_file: Path to BED file
    - output_dir: Directory to save plots
    """
    print("Generating CNV region plots...")

    # Load BED file containing annotation regions
    bed_df = pd.read_csv(bed_file, sep='\t', names=["chr", "start_pos", "end_pos", "gene"])

    # Get unique chromosomes with CNVs
    cnv_chrs = cnv_regions["chr"].unique()

    # Generate chromosome-wide plots
    for chrom in cnv_chrs:
        plot_region(data_frame, gene_lims, sample_name, z_threshold, cn_threshold, chrom, 0, float("inf"), bed_file, output_dir)

    # Determine zoomed-in CNV plot regions
    cnv_plot_regions = []
    
    for _, row in cnv_regions.iterrows():
        chrom = row["chr"]
        start_pos = row["start_pos"]
        end_pos = row["end_pos"]

        # Determine start and end positions for zoomed-in plot
        cnv_length = end_pos - start_pos
        cnv_length = max(cnv_length, 40)

        # Locate bed file regions covering the CNV
        bed_start_idx = bed_df[
            (bed_df["chr"] == chrom) & 
            (bed_df["start_pos"] <= start_pos) & 
            (bed_df["end_pos"] >= start_pos)
        ].index.min()

        bed_end_idx = bed_df[
            (bed_df["chr"] == chrom) & 
            (bed_df["start_pos"] <= end_pos) & 
            (bed_df["end_pos"] >= end_pos)
        ].index.max()

        # Adjust start and end positions
        plot_start_pos = bed_df.iloc[max(0, bed_start_idx - 3)]["start_pos"] if bed_start_idx is not None else start_pos - 2 * cnv_length
        plot_end_pos = bed_df.iloc[min(len(bed_df) - 1, bed_end_idx + 3)]["end_pos"] if bed_end_idx is not None else end_pos + 2 * cnv_length

        # Ensure plot regions do not overlap excessively
        if not cnv_plot_regions or (cnv_plot_regions[-1][1] < plot_start_pos):
            cnv_plot_regions.append((chrom, plot_start_pos, plot_end_pos))

    # Generate zoomed-in CNV plots
    for chrom, start_pos, end_pos in cnv_plot_regions:
        plot_region(data_frame, gene_lims, sample_name, z_threshold, cn_threshold, chrom, start_pos, end_pos, bed_file, output_dir)

def cyp2d6_report(data_frame, gene_lims, sample_name, z_threshold, cn_threshold, prop_threshold, bed_file, flag_bed_file, output_dir):
    """
    Generates a CYP2D6 CNV report.

    Parameters:
    - data_frame: Pandas DataFrame containing CNV data (Chr, Position, Depth, prop, mean, std, expDepth, zscore, copynumber)
    - gene_lims: DataFrame with gene start and end positions (chr, start, end, gene)
    - sample_name: Sample name for labeling
    - z_threshold: Z-score threshold
    - cn_threshold: Copy Number threshold
    - prop_threshold: Proportion threshold for CNV detection
    - bed_file: Path to BED file
    - flag_bed_file: Path to BED file containing flagged regions
    - output_dir: Directory to save reports and plots
    """
    print(f"Creating CYP2D6 CNV report for {sample_name}...")

    # Load BED and flagged regions data
    bed_df = pd.read_csv(bed_file, sep='\t', names=["chr", "start_pos", "end_pos", "gene"])
    flag_bed_df = pd.read_csv(flag_bed_file, sep='\t', names=["chr", "start_pos", "end_pos", "annot"])

    print(bed_df.head())
    print(flag_bed_df.head())
    print (data_frame.head())

    # Define the CYP2D6-D8 region boundaries
    cyp2d6_series = (data_frame["Chr"] == 22) & (data_frame["Position"] > 42123142) & (data_frame["Position"] < 42155251)

    # Identify CNV threshold exceedances
    threshold_series = (data_frame["zscore"].abs() > z_threshold) & ((data_frame["copynumber"] - 2).abs() > cn_threshold)

    report_file = os.path.join(output_dir, f"{sample_name}_CYP2D6_report_all.txt")
    with open(report_file, "w") as f:
        f.write("chr\tstart_pos\tend_pos\tgene\tcnv\tmean_depth\texp_mean_depth\tmean_cn\tmean_z-score\treg.length\tpos.above\t%above\tflag\tcomment\n")

        # Process CNVs within the CYP2D6 region
        for idx, row in data_frame.loc[cyp2d6_series].iterrows():
            chr_, pos = row["Chr"], row["Position"]

            # Find corresponding gene region
            gene_region = bed_df[(bed_df["chr"] == chr_) & (bed_df["start_pos"] < pos) & (bed_df["end_pos"] >= pos)]
            if gene_region.empty:
                continue

            start, end, gene_name = gene_region.iloc[0][["start_pos", "end_pos", "gene"]]

            # Extract relevant CNV data
            region_mask = (data_frame["Chr"] == chr_) & (data_frame["Position"] > start) & (data_frame["Position"] <= end)
            region_df = data_frame.loc[region_mask]

            mean_cn = region_df["copynumber"].mean()
            mean_z = region_df["zscore"].mean()
            mean_depth = region_df["Depth"].mean()
            mean_exp_depth = region_df["expDepth"].mean()

            cnv_type = "DUP" if mean_cn > 2 else "DEL" if mean_cn < 2 else "UNDETERMINED"
            region_length = len(region_df)
            pos_above = region_df.loc[threshold_series & region_mask].shape[0]
            percent_above = pos_above / region_length if region_length else 0

            flag = "PASS" if percent_above >= prop_threshold else "FAIL"
            comment = "FAIL - Low proportion of positions above threshold. " if flag == "FAIL" else ""

            # Flag low expected depth
            if mean_exp_depth < 50:
                flag = "FAIL"
                comment += "FAIL - Cohort mean depth too low. "

            # Flag overlaps with flagged regions
            overlap_flags = flag_bed_df[
                (flag_bed_df["chr"] == chr_) &
                (flag_bed_df["start_pos"] < end) &
                (flag_bed_df["end_pos"] > start)
            ]
            if not overlap_flags.empty:
                comment += "WARNING - Overlap with problematic region: " + "; ".join(
                    [f"{r['chr']}:{r['start_pos']}-{r['end_pos']} {r['annot']}" for _, r in overlap_flags.iterrows()]
                )

            # Write to report file
            f.write(f"{chr_}\t{start}\t{end}\t{gene_name}\t{cnv_type}\t{mean_depth:.2f}\t{mean_exp_depth:.2f}\t{mean_cn:.2f}\t"
                    f"{mean_z:.2f}\t{region_length}\t{pos_above}\t{percent_above:.2f}\t{flag}\t{comment}\n")

    print(f"CYP2D6 CNV report saved: {report_file}")

    # Load the report into a DataFrame
    cnv_regions = pd.read_csv(report_file, sep='\t')

    # Save a filtered report with only PASS flags
    filtered_report_file = os.path.join(output_dir, f"{sample_name}_CYP2D6_report_filtered.txt")
    cnv_regions_pass = cnv_regions[cnv_regions["flag"] == "PASS"]
    cnv_regions_pass.to_csv(filtered_report_file, sep='\t', index=False)
    print(f"Filtered report saved: {filtered_report_file}")

    # Create CYP2D6 vs CYP2D7 table
    element_df = cnv_regions.copy()
    element_df["gene"] = element_df["gene"].apply(lambda x: x.split("_")[0])
    element_df["element"] = element_df["gene"].apply(lambda x: x.split("_")[-1])

    elements = ["ex1", "ex2", "ex3", "ex4", "ex5", "ex6", "ex7", "ex8", "ex9", "REP"]
    cn_df = pd.DataFrame({"element": elements, "CYP2D6": 0.0, "CYP2D7": 0.0, "total": 0.0, "diff": 0.0, "ratio": 0.0, "CYP2D6_DEL_DUP": ""})

    for i, element in enumerate(elements):
        d6_values = element_df[(element_df["element"] == element) & (element_df["gene"] == "CYP2D6")]["mean_cn"]
        d7_values = element_df[(element_df["element"] == element) & (element_df["gene"] == "CYP2D7")]["mean_cn"]

        cn_df.at[i, "CYP2D6"] = d6_values.mean() if not d6_values.empty else 0.0
        cn_df.at[i, "CYP2D7"] = d7_values.mean() if not d7_values.empty else 0.0
        cn_df.at[i, "total"] = cn_df.at[i, "CYP2D6"] + cn_df.at[i, "CYP2D7"]
        cn_df.at[i, "diff"] = abs(cn_df.at[i, "CYP2D6"] - cn_df.at[i, "CYP2D7"])
        cn_df.at[i, "ratio"] = cn_df.at[i, "CYP2D6"] / cn_df.at[i, "CYP2D7"] if cn_df.at[i, "CYP2D7"] != 0 else np.nan

    # Save CN table
    cn_table_file = os.path.join(output_dir, f"{sample_name}_CYP2D6_vs_CYP2D7_table.txt")
    cn_df.to_csv(cn_table_file, sep='\t', index=False)
    print(f"CYP2D6 vs CYP2D7 table saved: {cn_table_file}")

    return cnv_regions

def plot_all_genes(data_frame, gene_lims, sample_name, z_threshold, cn_threshold, bed_file, output_dir):
    """
    Generates plots for all genes in the dataset.

    Parameters:
    - data_frame: Pandas DataFrame containing CNV data (Chr, Position, Depth, prop, mean, std, expDepth, zscore, copynumber)
    - gene_lims: DataFrame with gene start and end positions (chr, start, end, gene)
    - sample_name: Sample name for labeling
    - z_threshold: Z-score threshold
    - cn_threshold: Copy Number threshold
    - bed_file: Path to BED file
    - output_dir: Directory to save plots
    """
    print("Generating plots for all genes...")
    

    for index, row in gene_lims.iterrows():
        chr_, start, end, gene = row["chr"], row["start"], row["end"], row["gene"]
        print(f"Plotting gene: {gene} ({chr_}:{start}-{end})")
        
        # Call the plot_region function for each gene
        plot_region(data_frame, gene_lims, sample_name, z_threshold, cn_threshold, chr_, start, end, bed_file, output_dir)

def setup_logging(log_file):
    logging.basicConfig(
        filename=log_file,
        filemode='w',
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
    )
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    console.setFormatter(formatter)
    logging.getLogger().addHandler(console)

def main():
    """
    Main entry point for the script.

    This script takes in a sample name, input CNV file, input BED file, input flag BED file, z-score threshold, copy number threshold, proportion threshold, and chromosome lengths file as input and generates a report in the specified output folder.

    The report includes a CNV report, a plot of the CYP2D6 region, and plots for all genes in the dataset.

    The report directory is created if it does not exist, otherwise it is used as is.

    A log file is created in the report directory with the same name as the sample name.

    The script prints out the start time, the command line, and the end time to the log file.

    The script prints out the analysis steps to the log file.

    The script prints out the end time to the log file.
    """
    parser = argparse.ArgumentParser(description="Generate CYP2D6 CNV report from input files")
    parser.add_argument("input_sample_name", type=str, help="Prefix for output report and plots")
    parser.add_argument("input_cnv", type=str, help="Path to the input CNV file")
    parser.add_argument("input_bed", type=str, help="Path to the input BED file")
    parser.add_argument("input_flag_bed", type=str, help="Path to the input flag BED file")
    parser.add_argument("input_z_threshold", type=float, help="Threshold for the z-score")
    parser.add_argument("input_cn_threshold", type=float, help="Threshold for the copy number deviation from 2")
    parser.add_argument("input_prop_threshold", type=float, help="Proportion threshold")
    parser.add_argument("chr_lengths", type=str, help="Path to chromosome lengths file")
    parser.add_argument("output_folder", type=str, help="Folder to store the report output")
    
    args = parser.parse_args()

    print("****************************************************")
    print(f"Running Python CYP2D6 CNV Report Script with arguments: {args}")

    # Create report directory
    report_dir = os.path.join(args.output_folder, f"{args.input_sample_name}_CYP2D6_report")
    if not os.path.exists(report_dir):
        print(f"\nCreating report directory: {report_dir}")
        os.makedirs(report_dir)
    else:
        print(f"\nUsing existing report directory: {report_dir}")

    log_file = os.path.join(report_dir, f"{args.input_sample_name}_CYP2D6_report_log.txt")
    print(f"\nCreating log file: {log_file}")

    # Setup logging
    setup_logging(log_file)

    logging.info("****************************************************")
    logging.info(f"Running Python generate_CYP2D6_CNV_report.py {args.input_sample_name} {args.input_cnv} {args.input_bed} {args.input_flag_bed} {args.input_z_threshold} {args.input_cn_threshold} {args.input_prop_threshold} {args.chr_lengths} {args.output_folder}")
    logging.info(f"Start time: {datetime.now()}")

    # Read input CNV file
    logging.info(f"\nReading input CNV file: {args.input_cnv}")
    df_cnv = read_cnv_file(args.input_cnv)

    # Read BED file
    logging.info(f"\nFinding gene starts and ends from input BED file: {args.input_bed}")
    ##file_name = args.input_bed
    ##df_bed = gene_starts_and_ends(file_name)
    

    # Extract gene limits
    df_gene_limits = gene_starts_and_ends(args.input_bed)

    # Create CYP2D6 report
    logging.info("\nCreating CYP2D6 report...")
    cyp2d6_regions = cyp2d6_report(df_cnv, df_gene_limits, args.input_sample_name, args.input_z_threshold, args.input_cn_threshold, args.input_prop_threshold, args.input_bed, args.input_flag_bed, report_dir)

    # Create CYP2D6 region plot
    logging.info("\nCreating CYP2D6 region plot...")
    plot_region(df_cnv, df_gene_limits, args.input_sample_name, args.input_z_threshold, args.input_cn_threshold, "22", 42123143, 42155250, args.input_bed, report_dir)

    # Plot all genes
    logging.info("\nPlotting all genes...")
    plot_all_genes(df_cnv, df_gene_limits, args.input_sample_name, args.input_z_threshold, args.input_cn_threshold, args.input_bed, report_dir)

    # Final log statements
    logging.info("Analysis complete.")
    logging.info(f"End time: {datetime.now()}")

if __name__ == "__main__":
    main()
