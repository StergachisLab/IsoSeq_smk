#!/usr/bin/env Rscript
# Adriana Sedeno
# Create dictionary file

# Load necessary libraries
library(dplyr)
library(data.table)
library(stringr)
library(tidyr)
library(purrr)
library(janitor)

# Function to read and process the collapsed read statistics
process_collapsed <- function(filepath) {
  collapsed_data <- read.table(filepath, header = TRUE, col.names = c("read_id", "isoform"))
  return(collapsed_data)
}

# Function to read and process haplotag files
process_haplotag <- function(files) {
  haplotag_data <- files %>%
    set_names() %>%
    purrr::map_dfr(
      ~ data.table::fread(.x, header = TRUE, stringsAsFactors = FALSE, check.names = TRUE) %>% 
        dplyr::mutate(Filename = .x), 
      .id = "Filename"
    )
  return(haplotag_data)
}

# Function to read and process pigeon classification
process_classification <- function(filepath) {
  classification_data <- read.table(filepath, header = TRUE, sep = "\t") %>%
    janitor::row_to_names(row_number = 1) %>%
    select(isoform, structural_category, associated_gene, associated_transcript, subcategory)
  return(classification_data)
}

# Main function to create dictionary
create_dictionary <- function(collapsed_file, haplotag_files, classification_file) {
  # Process the input files
  collapsed_data <- process_collapsed(collapsed_file)
  haplotag_data <- process_haplotag(haplotag_files)
  classification_data <- process_classification(classification_file)
  
  dictionary <- collapsed_data %>%
    dplyr::left_join(classification_data) %>%
    dplyr::left_join(haplotag_data, by = join_by("read_id" == "Readname")) 
  
counts_hap <- dictionary %>%
  group_by(isoform, condition, haplotype) %>%
  count() %>%
  mutate(n_hap_cond = n) %>%
  select(-n) %>%
  pivot_wider(id_cols = c(isoform,Haplotype), names_from = c(condition), values_from = n_hap_cond) %>%
  mutate(iso_hap_noncyclo_counts = paste0("hn:i:", `non-cyclo`),
         iso_hap_cyclo_counts = paste0("hc:i:", cyclo)) %>%
  select(-cyclo, -`non-cyclo`)

counts <- dictionary %>%
  group_by(isoform, condition) %>%
  count() %>%
  mutate(n_cond = n) %>%
  select(-n) %>%
  pivot_wider(isoform, names_from = condition, values_from = n_cond) %>%
  mutate(iso_noncyclo_counts = paste0("sn:i:", `non-cyclo`),
         iso_cyclo_counts = paste0("sc:i:", cyclo)) %>%
  select(-cyclo, -`non-cyclo`)

dictionary <- left_join(dictionary, counts_hap) %>%
  left_join(counts) %>%
  mutate(isoform = gsub("PB.", "in:Z:", isoform),
         structural_category = paste0("sc:Z:",structural_category),
         associated_gene = paste0("gn:Z:", associated_gene),
         associated_transcript = paste0("tn:Z:", associated_transcript),
         subcategory = paste0("sb:Z:", subcategory)) %>%
    distinct(.) %>%
    dplyr::select(read_id, haplotype, isoform, structural_category, associated_gene, associated_transcript, subcategory)
  
  # Write the dictionary to a file
  dir.create("tag")
  write.table(dictionary, file = "tag/dictionary_all.txt", sep = "\t", row.names = FALSE, quote = FALSE)
}

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
collapsed_file <- args[1]
haplotag_files <- args[2]
classification_file <- args[3]

# Call the main function
create_dictionary(collapsed_file, haplotag_files, classification_file)
