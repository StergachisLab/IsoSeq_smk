#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(dplyr)
  library(purrr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript create_dictionary.R [collapsed_file] [haplotag_dir] [classification_file] [threads] [output_dir]")
}

collapsed_file       <- args[1]
haplotag_dir         <- args[2]
classification_file  <- args[3]
n_threads            <- as.integer(args[4])
output_dir           <- args[5]

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

collapsed_data <- fread(collapsed_file, header = TRUE, col.names = c("read_id", "isoform"))
classification_data <- fread(classification_file, header = TRUE, sep = "\t",
                             select = c("isoform", "structural_category", "associated_gene", "associated_transcript", "subcategory"))

hap <- list.files(haplotag_dir, pattern = "haplotagged.txt", full.names = TRUE) %>%
  set_names() %>%  
  map_dfr(~ fread(.x, header = TRUE, stringsAsFactors = FALSE, check.names = TRUE) %>%
            mutate(Filename = .x), .id = "Filename")

isoform_merge <- left_join(collapsed_data, classification_data) 
isoform_merge <- left_join(isoform_merge, hap, by = join_by("read_id" == "ReadName"))

dictionary <- isoform_merge %>%
    mutate(
    sample = str_extract(read_id, "PS\\d+"),
    condition = str_extract(read_id, "treated|untreated")
  )

counts_hap <- dictionary %>%
  group_by(sample, isoform, condition, Haplotype) %>%
  count() %>%
  mutate(n_hap_cond = n) %>%
  pivot_wider(id_cols = c(sample, isoform, Haplotype), names_from = c(condition), values_from = n_hap_cond) %>%
 mutate(across(everything(), ~replace(., is.na(.), 0))) %>%
  mutate(
    iso_hap_noncyclo_counts = if ("untreated" %in% colnames(.)) paste0("hn:i:", untreated) else "hn:i:0",
    iso_hap_cyclo_counts = if ("treated" %in% colnames(.)) paste0("hc:i:", treated) else "hc:i:0"
  ) %>%
  select(-any_of(c("treated", "untreated")))
  #mutate(iso_hap_noncyclo_counts = paste0("hn:i:", untreated),
  #       iso_hap_cyclo_counts = paste0("hc:i:", treated)) %>%
  #select(-treated, -untreated)

counts <- dictionary %>%
  group_by(sample, isoform, condition) %>%
  count() %>%
  mutate(n_cond = n) %>%
  pivot_wider(c(sample, isoform), names_from = condition, values_from = n_cond) %>%
 mutate(across(everything(), ~replace(., is.na(.), 0))) %>%
    mutate(
    iso_noncyclo_counts = if ("untreated" %in% colnames(.)) paste0("nc:i:", untreated) else "nc:i:0",
    iso_cyclo_counts = if ("treated" %in% colnames(.)) paste0("cc:i:", treated) else "cc:i:0"
  ) %>%
 select(-any_of(c("treated", "untreated")))
  #mutate(iso_noncyclo_counts = paste0("sn:i:", untreated),
  #       iso_cyclo_counts = paste0("sc:i:", treated)) %>%
  #select(-treated, -untreated)

dictionary <- left_join(dictionary, counts_hap) %>%
  left_join(counts) %>%
  mutate(across(everything(), ~replace(., is.na(.), 0))) %>%
  mutate(isoform = gsub("PB.", "in:Z:", isoform),
         structural_category = paste0("sc:Z:", structural_category),
         associated_gene = paste0("gn:Z:", associated_gene),
         associated_transcript = paste0("tn:Z:", associated_transcript),
         subcategory = paste0("sb:Z:", subcategory)) %>%
  distinct()

out <- dictionary %>%
 select(-ReadQuality, -ReadLength, -Haplotype, -Filename, -sample, -condition)

write_gzipped_tsv <- function(df, filename) {
  gz_filename <- paste0(filename, ".gz")
  con <- gzfile(gz_filename, "wt")  # Open gzipped file for writing
  write.table(df, con, sep = "\t", quote = FALSE, row.names = FALSE)
  close(con)
  return(gz_filename)
}

write_gzipped_tsv(out, "tag/dictionary.tsv")
