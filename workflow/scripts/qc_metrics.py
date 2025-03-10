import sys
import os
import pysam
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def extract_qc_metrics(bam_file):
    """Extract QC metrics from a given FLNC BAM file."""
    try:
        bam = pysam.AlignmentFile(bam_file, "rb")

        # Storage for metrics
        read_lengths = []
        isoform_counts = {}
        gene_counts = {}
        transcript_counts = {}
        sb_counts = {}
        st_counts = {}
        hp_counts = {}
        mapq_scores = []

        # Iterate through BAM reads
        for read in bam:
            read_lengths.append(read.query_length)

            # Extract tags (isoform, gene, transcript)
            isoform = read.get_tag("in") if read.has_tag("in") else None
            gene = read.get_tag("gn") if read.has_tag("gn") else None
            transcript = read.get_tag("tn") if read.has_tag("tn") else None

            if isoform:
                isoform_counts[isoform] = isoform_counts.get(isoform, 0) + 1
            if gene:
                gene_counts[gene] = gene_counts.get(gene, 0) + 1
            if transcript:
                transcript_counts[transcript] = transcript_counts.get(transcript, 0) + 1

            # Extract categorical tags (sb, st, HP)
            sb = read.get_tag("sb") if read.has_tag("sb") else None
            st = read.get_tag("st") if read.has_tag("st") else None
            hp = read.get_tag("HP") if read.has_tag("HP") else None

            if sb:
                sb_counts[sb] = sb_counts.get(sb, 0) + 1
            if st:
                st_counts[st] = st_counts.get(st, 0) + 1
            if hp:
                hp_counts[hp] = hp_counts.get(hp, 0) + 1

            # MAPQ Score
            mapq_scores.append(read.mapping_quality)

        bam.close()

        # Summarize QC metrics
        metrics = {
            "Total_FLNC_Reads": len(read_lengths),
            "Total_Unique_Isoforms": len(isoform_counts),
            "Total_Unique_Genes": len(gene_counts),
            "Total_Unique_Transcripts": len(transcript_counts),
        }

        return metrics, read_lengths, isoform_counts, gene_counts, transcript_counts, sb_counts, st_counts, hp_counts, mapq_scores

    except Exception as e:
        print(f"Error processing BAM file {bam_file}: {e}")
        sys.exit(1)


def plot_histogram(data, xlabel, ylabel, title, output_file, bins=50, log_scale=False, percentiles=None, color="blue"):
    """Generate histograms and save as PDF."""
    plt.figure(figsize=(8, 5))
    sns.histplot(data, bins=bins, kde=True, color=color, label=xlabel)

    if percentiles and len(data) > 0:
        for pct, col in percentiles.items():
            value = np.percentile(data, pct)
            plt.axvline(value, color=col, linestyle="dashed", label=f"{pct}%: {value:.0f}")

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()  # Ensuring the legend has labels

    if log_scale:
        plt.yscale("log")

    plt.savefig(output_file, format="pdf")
    plt.close()


def plot_bar_chart(df, x_label, y_label, title, output_file, log_scale=False, color="blue"):
    """Generate bar charts and save as PDF."""
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


if __name__ == "__main__":
    # Capture arguments
    if len(sys.argv) < 4:
        print("Usage: python qc_metrics.py <input_bam> <output_tsv> <output_plots_dir>")
        sys.exit(1)

    bam_file = sys.argv[1]
    output_tsv = sys.argv[2]
    output_dir = sys.argv[3]

    print(f"Processing BAM: {bam_file}")
    print(f"Saving output to: {output_tsv}")

    os.makedirs(output_dir, exist_ok=True)

    # Extract QC metrics
    qc_metrics, read_lengths, isoform_counts, gene_counts, transcript_counts, sb_counts, st_counts, hp_counts, mapq_scores = extract_qc_metrics(bam_file)

    # Save QC summary as TSV (avoid empty files)
    qc_df = pd.DataFrame(qc_metrics.items(), columns=["Metric", "Value"])
    if qc_df.empty:
        print(f"Warning: No QC data for {bam_file}. Writing placeholder.")
        qc_df = pd.DataFrame({"Metric": ["No Data"], "Value": ["N/A"]})

    qc_df.to_csv(output_tsv, sep="\t", index=False)

    # Compute percentiles for read length histogram
    percentiles = {10: "green", 50: "red", 90: "purple"}

    # Generate and save plots
    if read_lengths:
        plot_histogram(read_lengths, "FLNC Read Length", "Count", "FLNC Read Length Distribution",
                       os.path.join(output_dir, f"{os.path.basename(bam_file).replace('.bam', '_read_length.pdf')}"),
                       percentiles=percentiles, color="blue")

    if mapq_scores:
        plot_histogram(mapq_scores, "MAPQ Score", "Count", "MAPQ Score Distribution",
                       os.path.join(output_dir, f"{os.path.basename(bam_file).replace('.bam', '_mapq.pdf')}"),
                       bins=30, color="black")

    print(f"Finished processing {bam_file}")
