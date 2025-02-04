#!/usr/bin/env Rscript
# Adriana Sedeno
# Create dictionary files, one per haplotag.

##############################################################################
# create_dictionary.R
#
# Optimized to read collapsed_data and classification_data once, then
# process each haplotag file in parallel if >4 rows, creating a separate
# output dictionary per haplotag file.
#
# Usage:
#   Rscript create_dictionary.R [collapsed_file] [haplotag_dir] [classification_file] [threads] [output_dir]
#
# Example:
#   Rscript create_dictionary.R collapse/collapsed.read_stat.txt whatshap/ pigeon/pigeon_classification.txt 4 tag
#
##############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(parallel)     # or future.apply if you prefer cross-platform
})

# 1) Parse arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript create_dictionary.R [collapsed_file] [haplotag_dir] [classification_file] [threads] [output_dir]")
}

collapsed_file       <- args[1]
haplotag_dir         <- args[2]
classification_file  <- args[3]
n_threads            <- as.integer(args[4])
output_dir           <- args[5]

cat("collapsed_file:      ", collapsed_file,      "\n")
cat("haplotag_dir:        ", haplotag_dir,        "\n")
cat("classification_file: ", classification_file, "\n")
cat("threads:             ", n_threads,           "\n")
cat("output_dir:          ", output_dir,          "\n\n")

# Ensure output folder exists
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 2) Load the large files once ----------------------------------------------
cat("Reading collapsed_data...\n")
collapsed_data <- fread(
  collapsed_file,
  header = TRUE,
  col.names = c("read_id", "isoform")
)

cat("Reading classification_data...\n")
classification_data <- fread(
  classification_file,
  header = TRUE,
  sep = "\t",
  select = c("isoform", "structural_category", "associated_gene", "associated_transcript", "subcategory")
)

# 3) Define the core logic as a function ------------------------------------
create_dictionary_sample <- function(collapsed_data, classification_data, hap_data) {
  # Convert everything to data.table (if not already)
  setDT(collapsed_data)
  setDT(classification_data)
  setDT(hap_data)

  # 1) Filter
  # We only keep hap_data that merges with collapsed_data
  #   i.e. read_id in collapsed_data matches ReadName in hap_data
  setkey(collapsed_data, read_id)
  setkey(hap_data, ReadName)
  collapsed_sub <- hap_data[collapsed_data, nomatch = 0L, on = .(ReadName = read_id)]

  # Then classification_data must match isoform in collapsed_sub
  isoform_keep <- unique(collapsed_sub$isoform)
  classification_sub <- classification_data[isoform %in% isoform_keep]

  # 2) Join everything
  setnames(collapsed_sub, "ReadName", "read_id_original")
  dictionary <- merge(
    x = collapsed_sub,
    y = classification_sub,
    by = "isoform",
    all.x = TRUE
  )

  dictionary <- merge(
    x = dictionary,
    y = hap_data,
    by.x = "read_id",
    by.y = "ReadName",
    all.x = TRUE,
    suffixes = c("", "_hap")
  )

  # Replace NA with 0
  for (j in names(dictionary)) {
    set(dictionary, which(is.na(dictionary[[j]])), j, 0)
  }

  # 3) Rename columns
  dictionary[, isoform := sub("PB\\.", "in:Z:", isoform)]
  dictionary[, structural_category := paste0("sc:Z:", structural_category)]
  dictionary[, associated_gene := paste0("gn:Z:", associated_gene)]
  dictionary[, associated_transcript := paste0("tn:Z:", associated_transcript)]
  dictionary[, subcategory := paste0("sb:Z:", subcategory)]

  # Distinct rows
  setkeyv(dictionary, c("read_id", "isoform", "haplotype", "condition", "sample",
                        "structural_category", "associated_gene", "associated_transcript", "subcategory"))
  dictionary <- unique(dictionary)

  # 4) Build additional stats
  # a) counts_hap
  counts_hap <- dictionary[, .(n_hap_cond = .N), by = .(sample, isoform, condition, haplotype)]
  # Reshape wide
  counts_hap_wide <- dcast(
    counts_hap,
    sample + isoform + haplotype ~ condition,
    value.var = "n_hap_cond",
    fill = 0
  )
  if (!("cyclo" %in% names(counts_hap_wide))) counts_hap_wide[, cyclo := 0]
  if (!("non-cyclo" %in% names(counts_hap_wide))) counts_hap_wide[, `non-cyclo` := 0]

  counts_hap_wide[, iso_hap_noncyclo_counts := paste0("hn:i:", `non-cyclo`)]
  counts_hap_wide[, iso_hap_cyclo_counts := paste0("hc:i:", cyclo)]
  counts_hap_wide[, c("cyclo", "non-cyclo") := NULL]

  # b) counts
  counts_ <- dictionary[, .(n_cond = .N), by = .(sample, isoform, condition)]
  counts_wide <- dcast(
    counts_,
    sample + isoform ~ condition,
    value.var = "n_cond",
    fill = 0
  )
  if (!("cyclo" %in% names(counts_wide))) counts_wide[, cyclo := 0]
  if (!("non-cyclo" %in% names(counts_wide))) counts_wide[, `non-cyclo` := 0]

  counts_wide[, iso_noncyclo_counts := paste0("hn:i:", `non-cyclo`)]
  counts_wide[, iso_cyclo_counts := paste0("hc:i:", cyclo)]
  counts_wide[, c("cyclo", "non-cyclo") := NULL]

  # Merge back
  dictionary <- merge(dictionary, counts_hap_wide,
                      by = c("sample", "isoform", "haplotype"), all.x = TRUE)
  dictionary <- merge(dictionary, counts_wide,
                      by = c("sample", "isoform"), all.x = TRUE)

  # Final reorder
  out_cols <- c("read_id", "haplotype", "isoform", "condition", "sample",
                "structural_category", "associated_gene", "associated_transcript", "subcategory",
                "iso_hap_noncyclo_counts", "iso_hap_cyclo_counts",
                "iso_noncyclo_counts", "iso_cyclo_counts")
  out_cols <- intersect(out_cols, names(dictionary))
  setcolorder(dictionary, out_cols)

  return(dictionary)
}

# 4) List haplotag files, keep only those with >4 rows ----------------------
haplotag_files <- list.files(haplotag_dir, pattern = "haplotagged.txt", full.names = TRUE)
cat("Found", length(haplotag_files), "haplotag files in", haplotag_dir, "\n")

valid_files <- vector("list", length(haplotag_files))
valid_count <- 0

cat("Checking number of rows in each haplotag file...\n")
for (f in haplotag_files) {
  # For speed, we can quickly read the first 5 lines to see if it has >4 data lines
  # Or read the entire file if it's not huge. We'll do a minimal approach:
  n_lines <- as.integer(system2("wc", c("-l", shQuote(f)), stdout = TRUE))
  # If we trust that all lines are data lines (has a header?), adjust logic accordingly
  if (n_lines > 5) {
    valid_count <- valid_count + 1
    valid_files[[valid_count]] <- f
  }
}
valid_files <- valid_files[seq_len(valid_count)]
cat("Found", valid_count, "haplotag files with >=5 lines.\n\n")

# 5) Process in parallel per file -------------------------------------------
# If on Windows, you can't use mclapply reliably. Switch to lapply or future_lapply.
# For Mac/Linux HPC, mclapply is usually fine.

num_cores <- min(n_threads, valid_count)
cat("Using", num_cores, "cores\n")

mclapply(valid_files, function(f) {
  # Read hap_data
  hap_data <- fread(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  # If < 5 rows in the actual data (not counting the header), skip
  # (We've already done a quick check, but let's be safe.)
  if (nrow(hap_data) < 5) {
    message("Skipping file with <5 data rows: ", f)
    return(NULL)
  }

  sample_name <- sub("^(.*)_(?:treated|untreated).*", "\\1", basename(f))
  sample_name <- sub(".haplotagged.txt", "", sample_name)  # remove trailing extension
  # parse sample from ReadName:
  hap_data[, sample := sub("^(.*?)_(?:treated|untreated)_.*", "\\1", ReadName)]
  hap_data[, condition := fifelse(str_detect(ReadName, "_cyclo"), "cyclo",
                           fifelse(str_detect(ReadName, "_non-cyclo"), "non-cyclo", NA_character_))]
  hap_data[, identifier := str_extract(ReadName, "PS[^_]+")]

  # Create the dictionary data.table
  dict_dt <- create_dictionary_sample(collapsed_data, classification_data, hap_data)

  # Write output in <output_dir>/dictionary_<basename>.txt
  # Or we can do dictionary_{sample_name}.txt if you want one dictionary per haplotag file.
  out_base <- paste0("dictionary_", sample_name, ".txt")
  out_file <- file.path(output_dir, out_base)

  message("Writing dictionary for file: ", f, "\n  -> ", out_file)
  fwrite(dict_dt, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE)
  return(out_file)
}, mc.cores = num_cores)

cat("\nDone! All dictionary files created in:", output_dir, "\n")
