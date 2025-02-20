import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import os


def read_bed_file(bed_file: str) -> pd.DataFrame:
    """Read BED file and return as a DataFrame."""
    return pd.read_csv(
        bed_file,
        sep="\t",
        header=None,
        names=["Chr", "Start", "End", "Annot"],
        dtype={"Chr": str},
    )


def read_depth_file(depth_file: str) -> pd.DataFrame:
    """Read depth file and return as a DataFrame."""
    return pd.read_csv(
        depth_file,
        sep="\t",
        header=None,
        names=["Chr", "Position", "Depth"],
        dtype={"Chr": str},
    )


def read_snp_bed_file(snp_bed_file: str) -> pd.DataFrame:
    """Read SNP BED file and return as a DataFrame."""
    return pd.read_csv(
        snp_bed_file,
        sep="\t",
        header=None,
        names=["Chr", "Start", "End", "SNP_ID"],
        dtype={"Chr": str},
    )


def read_pc_bed_file(pc_bed_file: str) -> pd.DataFrame:
    """Read PC BED file and return as a DataFrame."""
    return pd.read_csv(
        pc_bed_file,
        sep="\t",
        header=None,
        names=["Chr", "Start", "ID", "Ref", "Alt"],
        dtype={"Chr": str},
    )


def annotate_depths(depth_df: pd.DataFrame, bed_df: pd.DataFrame) -> pd.DataFrame:
    """Annotate depth file with gene names from BED file."""
    bed_df = bed_df.sort_values(by=["Chr", "Start"])
    annotated_df = depth_df.copy()
    annotated_df["Gene"] = "ND"
    annotated_df["Element"] = "ND"

    for idx, row in annotated_df.iterrows():
        match = bed_df[
            (bed_df["Chr"] == row["Chr"])
            & (bed_df["Start"] <= row["Position"])
            & (bed_df["End"] >= row["Position"])
        ].copy()
        if not match.empty:
            annotated_df.at[idx, "Gene"] = match.iloc[0]["Annot"].split("/")[0]
            annotated_df.at[idx, "Element"] = match.iloc[0]["Annot"].split("/")[-1]

    return annotated_df


def plot_gene(
    annot_df: pd.DataFrame,
    sample_name: str,
    gene_name: str,
    depth_threshold: int,
    outfolder: str,
):
    """Generate a gene-specific coverage plot."""
    gene_df = annot_df[annot_df["Gene"] == gene_name]
    plt.figure(figsize=(12, 6))
    plt.scatter(
        gene_df["Position"],
        gene_df["Depth"],
        c=(gene_df["Depth"] < depth_threshold),
        cmap="coolwarm",
        s=5,
    )
    plt.axhline(
        y=depth_threshold,
        color="gray",
        linestyle="--",
        label=f"Threshold {depth_threshold}",
    )
    plt.xlabel("Genomic Positions")
    plt.ylabel("Read Depth")
    plt.title(f"{sample_name} - {gene_name}")
    plt.legend()
    plt.grid(True, linestyle="--", linewidth=0.5)
    plt.savefig(f"{outfolder}/{sample_name}_depth_{gene_name}.png")
    plt.close()


def generate_html_report(sample_id: str, coverage_stats: dict, output_dir: str):
    """Generate an HTML report summarizing the coverage statistics."""
    html_file = os.path.join(output_dir, f"{sample_id}_coverage_report.html")
    with open(html_file, "w") as f:
        f.write("<html><head><title>Coverage Report</title></head><body>")
        f.write(f"<h1>Coverage Report for {sample_id}</h1>")
        f.write("<h2>Overall Coverage</h2>")
        f.write("<table border='1'><tr><th>Metric</th><th>Value</th></tr>")
        for key, value in coverage_stats.items():
            f.write(
                f"<tr><td>{key.replace('_', ' ').capitalize()}</td><td>{value}</td></tr>"
            )
        f.write("</table>")
        f.write("</body></html>")
    print(f"HTML report saved: {html_file}")


def calculate_coverage(
    annot_df: pd.DataFrame, snp_df: pd.DataFrame, pc_df: pd.DataFrame, d_threshold: int
) -> dict:
    """Calculate coverage statistics above a given threshold, including SNP and PharmCat sites."""
    thres_levels = [80, 50, 20, 10]
    total_positions = len(annot_df)
    coverage_stats = {
        "total_positions": total_positions,
        "above_threshold": (annot_df["Depth"] >= d_threshold).sum(),
    }
    for thres in thres_levels:
        coverage_stats[f"above_{thres}"] = (annot_df["Depth"] >= thres).sum()

    # Calculate SNP coverage
    snp_coverage = annot_df[annot_df["Position"].isin(snp_df["Start"])]
    coverage_stats["snp_covered"] = (snp_coverage["Depth"] >= d_threshold).sum()

    # Calculate PharmCat site coverage
    pc_coverage = annot_df[annot_df["Position"].isin(pc_df["Start"])]
    coverage_stats["pharmcat_covered"] = (pc_coverage["Depth"] >= d_threshold).sum()

    coverage_stats = {
        k: round(100 * v / total_positions, 2) if "above" in k else v
        for k, v in coverage_stats.items()
    }
    return coverage_stats


def main():
    parser = argparse.ArgumentParser(description="Coverage Analysis Tool")
    parser.add_argument("--sample_id", type=str, help="Sample Identifier")
    parser.add_argument("--depth_file", type=str, help="Depth file path")
    parser.add_argument("--bed_file", type=str, help="BED file path")
    parser.add_argument("--snp_bed_file", type=str, help="SNP BED file path")
    parser.add_argument("--pc_bed_file", type=str, help="PharmCat BED file path")
    parser.add_argument("--threshold", type=int, help="Coverage threshold")
    parser.add_argument("--output_dir", type=str, help="Output directory for results")
    args = parser.parse_args()

    # Read input files
    depth_df = read_depth_file(args.depth_file)
    bed_df = read_bed_file(args.bed_file)
    snp_bed_df = read_snp_bed_file(args.snp_bed_file)
    pc_bed_df = read_pc_bed_file(args.pc_bed_file)

    # Annotate depths
    annotated_df = annotate_depths(depth_df, bed_df)

    # Calculate coverage statistics
    coverage_stats = calculate_coverage(
        annotated_df, snp_bed_df, pc_bed_df, args.threshold
    )

    # Generate HTML report
    generate_html_report(args.sample_id, coverage_stats, args.output_dir)

    # Generate gene coverage plots
    unique_genes = annotated_df["Gene"].unique()
    for gene in unique_genes:
        if gene != "ND":
            plot_gene(
                annotated_df, args.sample_id, gene, args.threshold, args.output_dir
            )

    # Save annotated depths
    output_file = f"{args.output_dir}/{args.sample_id}_depth_annotated.csv"
    annotated_df.to_csv(output_file, index=False)
    print(f"Annotated depth file saved: {output_file}")

    # Save coverage statistics
    stats_file = f"{args.output_dir}/{args.sample_id}_coverage_stats.txt"
    with open(stats_file, "w") as f:
        f.write(f"Coverage Statistics for {args.sample_id}\n")
        for key, value in coverage_stats.items():
            f.write(f"{key.replace('_', ' ').capitalize()}: {value}\n")
    print(f"Coverage statistics saved: {stats_file}")


if __name__ == "__main__":
    main()
