#!/usr/bin/env python3

import sys
import os
import re
import glob
import pandas as pd

# ------------------
# 1) Helper functions
# ------------------

def process_collapsed(filepath):
    """
    Replicates R's process_collapsed():
      read.table(filepath, header=TRUE, col.names=c("read_id","isoform"))
    """
    # If the file truly has headers, you can use header=True and remove the names=...
    # Adjust if your input file is formatted differently.
    collapsed_data = pd.read_csv(filepath, sep="\t", header=0,
                                 names=["read_id", "isoform"])
    return collapsed_data


def process_classification(filepath):
    """
    Replicates R's process_classification():
      read.table(..., header=TRUE, sep="\t") %>%
      select(isoform, structural_category, associated_gene, associated_transcript, subcategory)
    """
    classification_data = pd.read_csv(filepath, sep="\t", header=1)
    classification_data = classification_data[
        ["isoform", "structural_category",
         "associated_gene", "associated_transcript", "subcategory"]
    ]
    return classification_data


def create_dictionary_sample(collapsed_data, classification_data, hap_data):
    """
    Replicates R's create_dictionary_sample() steps:
      1) semi_join collapsed/hap_data
      2) semi_join classification/collapsed_sub
      3) left_join everything
      4) fill NA with 0
      5) rename_with(tolower)
      6) string replacements
      7) distinct()
      8) final select(...)
      9) pivot/wider for condition stats
      10) final merges with pivoted stats
    """

    # 1) collapsed_sub <- collapsed_data %>%
    #      semi_join(hap_data, by=join_by("read_id"=="ReadName"))
    # In Python, replicate semi_join via .isin():
    collapsed_sub = collapsed_data[
        collapsed_data["read_id"].isin(hap_data["ReadName"])
    ]

    # 2) classification_sub <- classification_data %>%
    #      semi_join(collapsed_sub, by="isoform")
    classification_sub = classification_data[
        classification_data["isoform"].isin(collapsed_sub["isoform"])
    ]

    # 3) dictionary <- collapsed_sub %>%
    #      left_join(classification_sub, by="isoform") %>%
    #      left_join(hap_data, by=join_by("read_id"=="ReadName"))
    dictionary = pd.merge(
        collapsed_sub,
        classification_sub,
        on="isoform",
        how="left"
    )
    dictionary = pd.merge(
        dictionary,
        hap_data,
        left_on="read_id",
        right_on="ReadName",
        how="left"
    )

    # 4) fill NA with 0 (replicate mutate_all(~ replace(., is.na(.), 0)))
    dictionary = dictionary.fillna(0)

    # 5) rename_with(tolower)
    dictionary.columns = [col.lower() for col in dictionary.columns]

    # 6) string replacements
    dictionary["isoform"] = dictionary["isoform"].astype(str).str.replace("PB.", "in:Z:")
    dictionary["structural_category"] = "sc:Z:" + dictionary["structural_category"].astype(str)
    dictionary["associated_gene"] = "gn:Z:" + dictionary["associated_gene"].astype(str)
    dictionary["associated_transcript"] = "tn:Z:" + dictionary["associated_transcript"].astype(str)
    dictionary["subcategory"] = "sb:Z:" + dictionary["subcategory"].astype(str)

    # 7) distinct()
    dictionary = dictionary.drop_duplicates()

    # 8) final select(...)
    # NOTE: If your haplotag data does NOT have "haplotype" or "condition" or "sample",
    # you must create them or remove them here.
    dictionary = dictionary[
        ["read_id", "haplotype", "isoform", "condition", "sample",
         "structural_category", "associated_gene", "associated_transcript", "subcategory"]
    ]

    # 9) Build additional stats (pivot/wider as in R)
    #    counts_hap -> pivot -> iso_hap_noncyclo_counts, iso_hap_cyclo_counts
    counts_hap = (
        dictionary.groupby(["sample", "isoform", "condition", "haplotype"])
        .size()
        .reset_index(name="n_hap_cond")
    )

    counts_hap_wide = counts_hap.pivot_table(
        index=["sample", "isoform", "haplotype"],
        columns="condition",
        values="n_hap_cond",
        fill_value=0
    ).reset_index()

    if "cyclo" not in counts_hap_wide.columns:
        counts_hap_wide["cyclo"] = 0
    if "non-cyclo" not in counts_hap_wide.columns:
        counts_hap_wide["non-cyclo"] = 0

    counts_hap_wide["iso_hap_noncyclo_counts"] = "hn:i:" + counts_hap_wide["non-cyclo"].astype(str)
    counts_hap_wide["iso_hap_cyclo_counts"] = "hc:i:" + counts_hap_wide["cyclo"].astype(str)
    counts_hap_wide = counts_hap_wide.drop(columns=["cyclo", "non-cyclo"])

    # counts -> pivot -> iso_noncyclo_counts, iso_cyclo_counts
    counts = (
        dictionary.groupby(["sample", "isoform", "condition"])
        .size()
        .reset_index(name="n_cond")
    )

    counts_wide = counts.pivot_table(
        index=["sample", "isoform"],
        columns="condition",
        values="n_cond",
        fill_value=0
    ).reset_index()

    if "cyclo" not in counts_wide.columns:
        counts_wide["cyclo"] = 0
    if "non-cyclo" not in counts_wide.columns:
        counts_wide["non-cyclo"] = 0

    counts_wide["iso_noncyclo_counts"] = "hn:i:" + counts_wide["non-cyclo"].astype(str)
    counts_wide["iso_cyclo_counts"] = "hc:i:" + counts_wide["cyclo"].astype(str)
    counts_wide = counts_wide.drop(columns=["cyclo", "non-cyclo"])

    # 10) final merges with pivoted stats
    dictionary = pd.merge(
        dictionary,
        counts_hap_wide,
        on=["sample", "isoform", "haplotype"],
        how="left"
    )
    dictionary = pd.merge(
        dictionary,
        counts_wide,
        on=["sample", "isoform"],
        how="left"
    )

    return dictionary


# -------------------------
# 2) Main script logic
# -------------------------

def main():
    # Parse CLI arguments
    args = sys.argv[1:]
    collapsed_file = args[0]
    haplotag_dir = args[1]
    classification_file = args[2]

    print("collapsed_file:", collapsed_file)
    print("haplotag_dir:  ", haplotag_dir)
    print("classification:", classification_file)

    # Load main data
    collapsed_data = process_collapsed(collapsed_file)
    classification_data = process_classification(classification_file)

    # Create output folder if needed
    os.makedirs("tag", exist_ok=True)

    # Locate all haplotagged files
    haplotag_files = glob.glob(os.path.join(haplotag_dir, "*haplotagged.txt"))
    print(f"Found {len(haplotag_files)} haplotag files.")

    # Build small DataFrame of file_path, sample
    df_files = pd.DataFrame({"file_path": haplotag_files})

    # replicate R's mutate(sample = sub("^(.*)_(?:treated|untreated).*", "\\1"))
    df_files["sample"] = df_files["file_path"].apply(
        lambda x: re.sub(r"^(.*)_(?:treated|untreated).*", r"\1", os.path.basename(x))
    )
    df_files["sample"] = df_files["sample"].str.replace(r"whatshap//", "", regex=True)

    # group by sample
    files_by_sample = df_files.groupby("sample")["file_path"].apply(list).reset_index()

    for _, row in files_by_sample.iterrows():
        sample_name = row["sample"]
        these_files = row["file_path"]
        print("Processing sample:", sample_name)
        print("  Haplotag files for this sample:", these_files)

        hap_data_list = []
        for f in these_files:
            # Read data
            tmp = pd.read_csv(f, sep="\t", header=0)

            # Check if empty
            if tmp.empty:
                print(f"File {f} has no data (empty). Skipping.")
                continue

            # replicate R mutate calls
            tmp["sample"] = tmp["ReadName"].str.replace(
                r"^(.*?)_(?:treated|untreated)_.*", r"\1", regex=True
            )
            tmp["condition"] = tmp["ReadName"].apply(
                lambda rn: "cyclo" if "_cyclo" in rn
                else ("non-cyclo" if "_non-cyclo" in rn else None)
            )
            tmp["identifier"] = tmp["ReadName"].str.extract(r"(PS[^_]+)")

            hap_data_list.append(tmp)

        # If all files are empty => skip creating dictionary
        if not hap_data_list:
            print(f"No data for sample {sample_name}. Skipping dictionary creation.")
            continue

        hap_data = pd.concat(hap_data_list, ignore_index=True)

        # Create dictionary for this sample
        dict_one = create_dictionary_sample(collapsed_data, classification_data, hap_data)

        # Write out
        out_name = f"tag/dictionary_{sample_name}.txt"
        dict_one.to_csv(out_name, sep="\t", index=False)
        print("Wrote:", out_name)

    print("Done!")


if __name__ == "__main__":
    main()
