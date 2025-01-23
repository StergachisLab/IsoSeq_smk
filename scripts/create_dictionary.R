#!/usr/bin/env Rscript
# Adriana Sedeno
# Create dictionary files, one per haplotag.

library(dplyr)
library(data.table)
library(stringr)
library(tidyr)

# ------------------
# 1) Helper functions
# ------------------

## Process the "collapsed" data file
process_collapsed <- function(filepath) {
  collapsed_data <- read.table(filepath, header = TRUE,
                               col.names = c("read_id", "isoform"))
  return(collapsed_data)
}

## Process the "classification" file
process_classification <- function(filepath) {
  classification_data <- read.table(filepath, header = TRUE, sep = "\t") %>%
    dplyr::select(isoform, structural_category, associated_gene, associated_transcript, subcategory)
  return(classification_data)
}

## The core logic to merge collapsed + classification + a single haplotag
create_dictionary_sample <- function(
  collapsed_data, classification_data, hap_data
) {

collapsed_sub <- collapsed_data %>%
    dplyr::semi_join(hap_data, by = join_by("read_id" == "ReadName"))

classification_sub <- classification_data %>%
    dplyr::semi_join(collapsed_sub, by = "isoform")

  # 2) Join everything
  dictionary <- collapsed_sub %>%
    dplyr::left_join(classification_sub, by = "isoform") %>%
    dplyr::left_join(hap_data, by = dplyr::join_by("read_id" == "ReadName")) %>%
    mutate_all(~ replace(., is.na(.), 0)) %>%
    rename_with(tolower) %>%
    dplyr::mutate(
      isoform = gsub("PB.", "in:Z:", isoform),
      structural_category = paste0("sc:Z:", structural_category),
      associated_gene = paste0("gn:Z:", associated_gene),
      associated_transcript = paste0("tn:Z:", associated_transcript),
      subcategory = paste0("sb:Z:", subcategory)
    ) %>%
    distinct() %>%
    dplyr::select(read_id, haplotype, isoform, condition, sample,
                  structural_category, associated_gene,
                  associated_transcript, subcategory)

  # 4) Build additional stats
  counts_hap <- dictionary %>%
  group_by(sample, isoform, condition, haplotype) %>%
  summarize(n_hap_cond = n()) %>%
  ungroup() %>%
    tidyr::pivot_wider(
    id_cols = c(sample, isoform, haplotype),
    names_from = condition,
    values_from = n_hap_cond,
    values_fill = 0) %>%
  mutate(
    cyclo      = if ("cyclo"      %in% colnames(.)) cyclo      else 0,
    `non-cyclo` = if ("non-cyclo" %in% colnames(.)) `non-cyclo` else 0
  ) %>%
  mutate(
    iso_hap_noncyclo_counts = paste0("hn:i:", `non-cyclo`),
    iso_hap_cyclo_counts    = paste0("hc:i:", cyclo)
  ) %>%
  select(-cyclo, -`non-cyclo`)

counts <- dictionary %>%
  group_by(sample, isoform, condition) %>%
  summarize(n_cond = n()) %>%
  ungroup() %>%
    tidyr::pivot_wider(
    id_cols = c(sample, isoform),
    names_from = condition,
    values_from = n_cond,
    values_fill = 0) %>%
  mutate(
    cyclo      = if ("cyclo"      %in% colnames(.)) cyclo      else 0,
    `non-cyclo` = if ("non-cyclo" %in% colnames(.)) `non-cyclo` else 0
  ) %>%
  mutate(
    iso_noncyclo_counts = paste0("hn:i:", `non-cyclo`),
    iso_cyclo_counts    = paste0("hc:i:", cyclo)
  ) %>%
  select(-cyclo, -`non-cyclo`)


  dictionary <- dictionary %>%
    left_join(counts_hap) %>%
    left_join(counts)

  return(dictionary)
}

# -------------------------
# 2) Main: process arguments
# -------------------------

args <- commandArgs(trailingOnly = TRUE)
collapsed_file <- args[1]
haplotag_dir   <- args[2]
classification_file <- args[3]

cat("collapsed_file: ", collapsed_file, "\n")
cat("haplotag_dir:   ", haplotag_dir,   "\n")
cat("classification: ", classification_file, "\n")

# ------------------------
# 3) Read smaller input once
# ------------------------
collapsed_data <- process_collapsed(collapsed_file)
classification_data <- process_classification(classification_file)

# Create output folder if needed
dir.create("tag", showWarnings = FALSE)

# ----------------------------
# 4) Loop over haplotag files
# ----------------------------
haplotag_files <- list.files(haplotag_dir,
                             pattern = "haplotagged.txt",
                             full.names = TRUE)

cat("Found", length(haplotag_files), "haplotag files.\n")

df_files <- tibble(file_path = haplotag_files) %>%
  mutate(sample = sub("^(.*)_(?:treated|untreated).*", "\\1", file_path)) %>%
  mutate(sample = sub("whatshap//", "", sample))


# 3. Group by sample, then store the file paths in a list column
files_by_sample <- df_files %>%
  group_by(sample) %>%
  summarise(file_list = list(file_path))

for (i in seq_len(nrow(files_by_sample))) {
  sample_name <- files_by_sample$sample[i]
  these_files <- files_by_sample$file_list[[i]]

  cat("Processing sample:", sample_name, "\n")
  cat("  Haplotag files for this sample:", these_files, "\n")

  hap_data <- purrr::map_dfr(these_files, ~ {
    data.table::fread(.x, header = TRUE, stringsAsFactors = FALSE) %>%
      mutate(
        sample = sub("^(.*?)_(?:treated|untreated)_.*", "\\1", ReadName),
        condition = case_when(
          str_detect(ReadName, "_cyclo") ~ "cyclo",
          str_detect(ReadName, "_non-cyclo") ~ "non-cyclo",
          TRUE ~ NA_character_
        ),
        identifier = str_extract(ReadName, "PS[^_]+")
      )
  })

 dict_one <- create_dictionary_sample(
    collapsed_data,
    classification_data,
    hap_data = hap_data
  )

  # 4) Write out a single dictionary file for this sample
  out_name <- paste0("tag/dictionary_", sample_name, ".txt")
  write.table(dict_one, file = out_name,
              sep = "\t", row.names = FALSE, quote = FALSE)
}



cat("Done!\n")
