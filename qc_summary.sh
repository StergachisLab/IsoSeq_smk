#!/bin/bash

# Check if directory path was provided
if [ -z "$1" ]; then
  echo "Usage: $0 <path_to_qc_summary_files>"
  exit 1
fi

INPUT_DIR="$1"

# Print header
echo -e "File\tTotal_FLNC_Reads\tTotal_Unique_Isoforms\tTotal_Unique_Genes\tTotal_Unique_Transcripts"

# Loop through all *_qc_summary.tsv files in the provided directory
for f in "${INPUT_DIR}"/*_qc_summary.tsv; do
  # Extract metrics
  reads=$(awk '$1 == "Total_FLNC_Reads" {print $2}' "$f")
  isoforms=$(awk '$1 == "Total_Unique_Isoforms" {print $2}' "$f")
  genes=$(awk '$1 == "Total_Unique_Genes" {print $2}' "$f")
  transcripts=$(awk '$1 == "Total_Unique_Transcripts" {print $2}' "$f")
  
  # Print row
  echo -e "$(basename "$f")\t${reads}\t${isoforms}\t${genes}\t${transcripts}"
done
