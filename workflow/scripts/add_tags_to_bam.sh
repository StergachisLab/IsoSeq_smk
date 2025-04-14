#!/bin/bash
# Adriana Sedeno
# Optimized script to append tags to BAM file with minimal disk usage and memory overhead

# Check for the correct number of arguments
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <input.bam> <dictionary.tsv.gz> <output.bam> <threads>"
    exit 1
fi

input_bam=$1
dict_file=$2
output_bam=$3
threads=$4

# Ensure samtools is installed
if ! command -v samtools &> /dev/null; then
    echo "Error: samtools could not be found, please install samtools."
    exit 1
fi

echo "Processing BAM file: $input_bam"

# Step 1: Extract headers from BAM (streaming output to stdout)
samtools view -@ "$threads" -H "$input_bam" > "${output_bam}.header.sam"

# Step 2: Extract read IDs from BAM and filter dictionary dynamically
echo "Extracting read IDs from BAM and filtering dictionary..."
samtools view -@ "$threads" "$input_bam" | cut -f1 | sort -u | \
    awk 'NR==FNR { read_ids[$1]=1; next } ($1 in read_ids)' - <(gzip -dc "$dict_file") | \
    sort -k1,1 > "${output_bam}.filtered_dict.tsv"

# Step 3: Stream BAM reads and merge with dictionary dynamically
echo "Processing BAM reads and appending tags..."
samtools view -@ "$threads" "$input_bam" | sort -k1,1 | \
    join -t $'\t' -1 1 -2 1 - "${output_bam}.filtered_dict.tsv" | \
    cat "${output_bam}.header.sam" - | samtools view -@ "$threads" -b -o "${output_bam}.unsorted.bam"

# Step 4: Sort BAM before indexing
echo "Sorting BAM..."
samtools sort -@ "$threads" -o "$output_bam" "${output_bam}.unsorted.bam"

# Step 5: Index BAM
echo "Indexing BAM..."
samtools index -@ "$threads" "$output_bam"

# Cleanup temporary header file and filtered dictionary
rm "${output_bam}.header.sam" "${output_bam}.filtered_dict.tsv" "${output_bam}.unsorted.bam"

# Confirm completion
if [ $? -eq 0 ]; then
    echo "✅ Tags added and BAM file indexed successfully to $output_bam."
else
    echo "❌ An error occurred during the tagging process."
    exit 1
fi
