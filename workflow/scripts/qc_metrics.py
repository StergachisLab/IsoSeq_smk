#!/usr/bin/env python

import sys
import os
import pysam
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def extract_qc_metrics(bam_file, threads=8):
    """Extract QC metrics from a given FLNC BAM file."""
    try:
        bam = pysam.AlignmentFile(bam_file, "rb", threads=threads)

        read_lengths = []
        isoform_counts = {}
        gene_counts = {}
        transcript_counts = {}
        sb_counts = {}
        sc_counts = {}
        hp_counts = {}
        mapq_scores = []

        for read in bam:
            read_lengths.append(read.query_length)

            isoform = read.get_tag("in") if read.has_tag("in") else None
            gene = read.get_tag("gn") if read.has_tag("gn") else None
            transcript = read.get_tag("tn") if read.has_tag("tn") else None

            if isoform:
                isoform_counts[isoform] = isoform_counts.get(isoform, 0) + 1
            if gene:
                gene_counts[gene] = gene_counts.get(gene, 0) + 1
            if transcript:
                transcript_counts[transcript] = transcript_counts.get(transcript, 0) + 1

            sb = read.get_tag("sb") if read.has_tag("sb") else None
            sc = read.get_tag("sc") if read.has_tag("sc") else None
            hp = read.get_tag("HP") if read.has_tag("HP") else None

            if sb and sb != "0":
                sb_counts[sb] = sb_counts.get(sb, 0) + 1
            if sc and sc != "0":
                sc_counts[sc] = sc_counts.get(sc, 0) + 1
            if hp:
                hp_counts[hp] = hp_counts.get(hp, 0) + 1

            mapq_scores.append(read.mapping_quality)

        bam.close()

        metrics = {
            "Total_FLNC_Reads": len(read_lengths),
            "Total_Unique_Isoforms": len(isoform_counts),
            "Total_Unique_Genes": len(gene_counts),
            "Total_Unique_Transcripts": len(transcript_counts),
        }

        return (
            metrics,
            read_lengths,
            isoform_counts,
            gene_counts,
            transcript_counts,
            sb_counts,
            sc_counts,
            hp_counts,
            mapq_scores,
        )

    except Exception as e:
        print(f"Error processing BAM file {bam_file}: {e}")
        sys.exit(1)


def plot_histogram(
    data,
    xlabel,
    ylabel,
    title,
    output_file,
    bins=50,
    log_scale=False,
    percentiles=None,
    color="blue",
):
    plt.figure(figsize=(8, 5))
    sns.histplot(data, bins=bins, kde=True, color=color, label=xlabel)
    if percentiles and len(data) > 0:
        for pct, col in percentiles.items():
            value = np.percentile(data, pct)
            plt.axvline(
                value, color=col, linestyle="dashed", label=f"{pct}%: {value:.0f}"
            )
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    if log_scale:
        plt.yscale("log")
    plt.savefig(output_file, format="pdf")
    plt.close()


def plot_bar_chart(
    df, x_label, y_label, title, output_file, log_scale=False, color="blue"
):
    plt.figure(figsize=(10, 5))
    if df.empty:
        print(f"Skipping {output_file} (No data to plot)")
        return
    sns.barplot(x=np.arange(len(df)), y=df["Count"], color=color, label=title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    plt.legend()
    if log_scale:
        plt.yscale("log")
    plt.savefig(output_file, format="pdf")
    plt.close()


def plot_read_length_by_category(length_dict, category_name, output_file):
    plt.figure(figsize=(10, 6))
    for category, lengths in length_dict.items():
        if len(lengths) > 1:
            sns.kdeplot(lengths, label=category, fill=True, alpha=0.4)
    plt.xlabel("Read Length")
    plt.ylabel("Density")
    plt.title(f"Read Length Distribution by {category_name}")
    plt.legend(title=category_name)
    plt.tight_layout()
    plt.savefig(output_file, format="pdf")
    plt.close()


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python qc_metrics.py <input_bam> <output_tsv> <output_plots_dir> <threads>")
        sys.exit(1)

    bam_file = sys.argv[1]
    output_tsv = sys.argv[2]
    output_dir = sys.argv[3]
    threads = int(sys.argv[4]) if len(sys.argv) > 4 else 8

    print(f"Processing BAM: {bam_file}")
    print(f"Saving output to: {output_tsv}")
    os.makedirs(output_dir, exist_ok=True)

    (
        qc_metrics,
        read_lengths,
        isoform_counts,
        gene_counts,
        transcript_counts,
        sb_counts,
        sc_counts,
        hp_counts,
        mapq_scores,
    ) = extract_qc_metrics(bam_file, threads)

    qc_df = pd.DataFrame(qc_metrics.items(), columns=["Metric", "Value"])

    if read_lengths:
        qc_df_extra = pd.DataFrame(
            {
                "Metric": [
                    "Median_FLNC_Read_Length",
                    "Percentile_10_FLNC_Length",
                    "Percentile_90_FLNC_Length",
                ],
                "Value": [
                    int(np.median(read_lengths)),
                    int(np.percentile(read_lengths, 10)),
                    int(np.percentile(read_lengths, 90)),
                ],
            }
        )
        qc_df = pd.concat([qc_df, qc_df_extra], ignore_index=True)

    if qc_df.empty:
        print(f"Warning: No QC data for {bam_file}. Writing placeholder.")
        qc_df = pd.DataFrame({"Metric": ["No Data"], "Value": ["N/A"]})

    qc_df.to_csv(output_tsv, sep="\t", index=False)

    percentiles = {10: "green", 50: "red", 90: "purple"}

    if read_lengths:
        plot_histogram(
            read_lengths,
            "FLNC Read Length",
            "Count",
            "FLNC Read Length Distribution",
            os.path.join(
                output_dir,
                os.path.basename(bam_file).replace(".bam", "_read_length.pdf"),
            ),
            percentiles=percentiles,
            color="blue",
        )

    if mapq_scores:
        plot_histogram(
            mapq_scores,
            "MAPQ Score",
            "Count",
            "MAPQ Score Distribution",
            os.path.join(
                output_dir, os.path.basename(bam_file).replace(".bam", "_mapq.pdf")
            ),
            bins=30,
            color="black",
        )

    sc_length_dict = {}
    sb_length_dict = {}

    bam = pysam.AlignmentFile(bam_file, "rb")
    for read in bam:
        read_length = read.query_length
        if read_length is None:
            continue
        if read.has_tag("sc"):
            sc = read.get_tag("sc")
            if sc and sc != "0":
                sc_length_dict.setdefault(sc, []).append(read_length)
        if read.has_tag("sb"):
            sb = read.get_tag("sb")
            if sb and sb != "0":
                sb_length_dict.setdefault(sb, []).append(read_length)
    bam.close()

    plot_read_length_by_category(
        sc_length_dict,
        "Structural Category (sc)",
        os.path.join(
            output_dir,
            os.path.basename(bam_file).replace(".bam", "_read_length_by_sc.pdf"),
        ),
    )

    plot_read_length_by_category(
        sb_length_dict,
        "Subcategory (sb)",
        os.path.join(
            output_dir,
            os.path.basename(bam_file).replace(".bam", "_read_length_by_sb.pdf"),
        ),
    )

    sc_counts_df = pd.DataFrame(
        [
            {
                "structural_category": k,
                "read_count": len(v),
                "median_read_length": int(np.median(v)) if v else 0,
            }
            for k, v in sc_length_dict.items()
        ]
    )
    if not sc_counts_df.empty:
        sc_counts_df["percentage"] = (
            sc_counts_df["read_count"] / sc_counts_df["read_count"].sum() * 100
        ).round(2)

    sc_counts_df.to_csv(
        os.path.join(
            output_dir, os.path.basename(bam_file).replace(".bam", "_sc_counts.tsv")
        ),
        sep="\t",
        index=False,
    )

    sb_counts_df = pd.DataFrame(
        [
            {
                "subcategory": k,
                "read_count": len(v),
                "median_read_length": int(np.median(v)) if v else 0,
            }
            for k, v in sb_length_dict.items()
        ]
    )
    if not sb_counts_df.empty:
        sb_counts_df["percentage"] = (
            sb_counts_df["read_count"] / sb_counts_df["read_count"].sum() * 100
        ).round(2)

    sb_counts_df.to_csv(
        os.path.join(
            output_dir, os.path.basename(bam_file).replace(".bam", "_sb_counts.tsv")
        ),
        sep="\t",
        index=False,
    )

    print(f"Finished processing {bam_file}")
