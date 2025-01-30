#!/usr/bin/env python3
# Optimized conversion of R script to Python with Dask for Snakemake compatibility

import os
import sys
import dask.dataframe as dd
import pandas as pd
import numpy as np
import glob
import re
from dask import delayed, compute
from dask.diagnostics import ProgressBar

# -----------------------
# 1) Helper Functions
# -----------------------

## Read the collapsed file
def process_collapsed(filepath):
    return dd.read_csv(filepath, sep="\t", names=["read_id", "isoform"], header=0, dtype=str)

## Read the classification file
def process_classification(filepath):
    cols = ["isoform", "structural_category", "associated_gene", "associated_transcript", "subcategory"]
    return dd.read_csv(filepath, sep="\t", usecols=cols, header=0, dtype=str)

## Process the haplotagged files
def process_haplotagged(filepath):
    df = dd.read_csv(filepath, sep="\t", dtype=str)  # Lazy loading

    # Extract sample name from ReadName
    df["sample"] = df["ReadName"].str.extract(r"^(.*?)_(?:treated|untreated)_.*")
    
    # Assign condition
    df["condition"] = dd.np.where(df["ReadName"].str.contains("_cyclo"), "cyclo",
                       dd.np.where(df["ReadName"].str.contains("_non-cyclo"), "non-cyclo", None))

    # Extract identifier
    df["identifier"] = df["ReadName"].str.extract(r"(PS[^_]+)")
    
    return df

## Merge collapsed + classification + haplotagged data
def create_dictionary_sample(collapsed_data, classification_data, hap_data):
    collapsed_sub = collapsed_data.merge(hap_data, left_on="read_id", right_on="ReadName", how="inner")
    classification_sub = classification_data[classification_data["isoform"].isin(collapsed_sub["isoform"])]

    dictionary = (
        collapsed_sub
        .merge(classification_sub, on="isoform", how="left")
        .merge(hap_data, left_on="read_id", right_on="ReadName", how="left")
        .fillna(0)  # Replace NaNs with 0
    )

    # Format the columns using vectorized operations
    dictionary["isoform"] = dictionary["isoform"].str.replace("PB.", "in:Z:")
    dictionary["structural_category"] = "sc:Z:" + dictionary["structural_category"].astype(str)
    dictionary["associated_gene"] = "gn:Z:" + dictionary["associated_gene"].astype(str)
    dictionary["associated_transcript"] = "tn:Z:" + dictionary["associated_transcript"].astype(str)
    dictionary["subcategory"] = "sb:Z:" + dictionary["subcategory"].astype(str)

    dictionary = dictionary.drop_duplicates()
    dictionary = dictionary[["read_id", "haplotype", "isoform", "condition", "sample",
                             "structural_category", "associated_gene", "associated_transcript", "subcategory"]]

    # Aggregate counts
    counts_hap = dictionary.groupby(["sample", "isoform", "condition", "haplotype"]).size().reset_index(name="n_hap_cond")
    counts = dictionary.groupby(["sample", "isoform", "condition"]).size().reset_index(name="n_cond")

    counts_hap_wide = counts_hap.pivot_table(index=["sample", "isoform", "haplotype"],
                                             columns="condition",
                                             values="n_hap_cond", fill_value=0).reset_index()
    counts_wide = counts.pivot_table(index=["sample", "isoform"],
                                     columns="condition",
                                     values="n_cond", fill_value=0).reset_index()

    # Rename for final output
    for df in [counts_hap_wide, counts_wide]:
        df["iso_hap_noncyclo_counts"] = "hn:i:" + df.get("non-cyclo", 0).astype(str)
        df["iso_hap_cyclo_counts"] = "hc:i:" + df.get("cyclo", 0).astype(str)
        df.drop(columns=["cyclo", "non-cyclo"], errors="ignore", inplace=True)

    dictionary = dictionary.merge(counts_hap_wide, on=["sample", "isoform", "haplotype"], how="left")\
                           .merge(counts_wide, on=["sample", "isoform"], how="left")

    return dictionary


# -----------------------
# 2) Main Processing
# -----------------------

def process_sample(sample_name, files, collapsed_data, classification_data, output_dir):
    print(f"Processing sample: {sample_name}")
    
    # Read all haplotag files for the sample using Dask
    hap_data = dd.concat([process_haplotagged(file) for file in files])

    # Create dictionary
    dictionary = create_dictionary_sample(collapsed_data, classification_data, hap_data)

    # Compute and write output
    output_file = os.path.join(output_dir, f"dictionary_{sample_name}.txt")
    dictionary.compute().to_csv(output_file, sep="\t", index=False)
    print(f"Saved: {output_file}")

def main():
    # Get command line arguments
    collapsed_file = sys.argv[1]
    haplotag_dir = sys.argv[2]
    classification_file = sys.argv[3]

    print(f"Collapsed file: {collapsed_file}")
    print(f"Haplotag directory: {haplotag_dir}")
    print(f"Classification file: {classification_file}")

    # Load fixed input files with Dask
    collapsed_data = process_collapsed(collapsed_file)
    classification_data = process_classification(classification_file)

    # Create output folder if needed
    output_dir = "tag"
    os.makedirs(output_dir, exist_ok=True)

    # Find haplotagged files
    haplotag_files = glob.glob(os.path.join(haplotag_dir, "*haplotagged.txt"))
    print(f"Found {len(haplotag_files)} haplotag files.")

    # Group files by sample
    df_files = pd.DataFrame({"file_path": haplotag_files})
    df_files["sample"] = df_files["file_path"].apply(lambda x: re.sub(r"^(.*)_(?:treated|untreated).*", r"\1", os.path.basename(x)))

    files_by_sample = df_files.groupby("sample")["file_path"].apply(list).reset_index()

    # Use Dask's Delayed for parallel processing
    tasks = [delayed(process_sample)(row["sample"], row["file_path"], collapsed_data, classification_data, output_dir) for _, row in files_by_sample.iterrows()]

    with ProgressBar():
        compute(*tasks)

    print("Done!")

if __name__ == "__main__":
    main()
