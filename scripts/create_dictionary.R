#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(parallel)
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
  mutate(iso_hap_noncyclo_counts = paste0("hn:i:", `non-cyclo`),
         iso_hap_cyclo_counts = paste0("hc:i:", cyclo)) %>%
  select(-cyclo, -`non-cyclo`)

counts <- dictionary %>%
  group_by(sample, isoform, condition) %>%
  count() %>%
  mutate(n_cond = n) %>%
  pivot_wider(c(sample, isoform), names_from = condition, values_from = n_cond) %>%
  mutate(iso_noncyclo_counts = paste0("sn:i:", `non-cyclo`),
         iso_cyclo_counts = paste0("sc:i:", cyclo)) %>%
  select(-cyclo, -`non-cyclo`)

dictionary <- left_join(dictionary, counts_hap) %>%
  left_join(counts) %>%
  mutate(across(everything(), ~replace(., is.na(.), 0))) %>%
  mutate(isoform = gsub("PB.", "in:Z:", isoform),
         structural_category = paste0("sc:Z:", structural_category),
         associated_gene = paste0("gn:Z:", associated_gene),
         associated_transcript = paste0("tn:Z:", associated_transcript),
         subcategory = paste0("sb:Z:", subcategory)) %>%
  distinct()

samples_list <- unique(dictionary$sample)

mclapply(samples_list, function(sample_id) {
  dictionary_out <- dictionary %>%
    filter(sample == sample_id) %>%
    select(-ReadQuality, -ReadLength, -Haplotype, -Filename, -condition, -id)
  
  output_file <- file.path(output_dir, paste0("dictionary_", sample_id, ".txt"))
  write.table(dictionary_out, output_file, quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
}, mc.cores = n_threads)
