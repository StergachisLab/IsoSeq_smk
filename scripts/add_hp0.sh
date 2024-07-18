#!/bin/bash
# Adriana Sedeno
# 5/16/2024
# Append a fixed tag to a BAM file

# Check for the correct number of arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input.bam> <output.bam>"
    exit 1
fi

input_bam=$1
output_bam=$2
tmp_prefix="./temp_$$"

# Ensure samtools is installed
if ! command -v samtools &> /dev/null; then
    echo "samtools could not be found, please install samtools."
    exit 1
fi

# Convert BAM to SAM and preserve header
samtools view -h "$input_bam" > "${tmp_prefix}_temp.sam"
if [ $? -ne 0 ]; then
    echo "Error: samtools view command failed."
    exit 1
fi

# Append the fixed tag HP:i:0 to each read in the SAM file
awk 'BEGIN { FS=OFS="\t" }
     /^@/ { print; next }
     { print $0, "HP:i:0" }
    ' "${tmp_prefix}_temp.sam" > "${tmp_prefix}_temp_with_tags.sam"
if [ $? -ne 0 ]; then
    echo "Error: awk command failed."
    exit 1
fi

# Convert the modified SAM back to BAM
samtools view -b -o "$output_bam" "${tmp_prefix}_temp_with_tags.sam"
if [ $? -ne 0 ]; then
    echo "Error: samtools view command failed."
    exit 1
fi

# Index the BAM file
samtools index "$output_bam"
if [ $? -ne 0 ]; then
    echo "Error: samtools index command failed."
    exit 1
fi

# Cleanup temporary files
rm "${tmp_prefix}_temp.sam" "${tmp_prefix}_temp_with_tags.sam"

# Confirm completion
echo "Tag HP:i:0 added and BAM file indexed successfully to $output_bam."