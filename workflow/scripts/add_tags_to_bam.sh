#!/bin/bash
# Adriana Sedeno
# Optimized script to append tags to BAM file while preserving integrity

# Check for the correct number of arguments
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <input.bam> <dictionary.tsv.gz> <output.bam> <threads>"
    exit 1
fi

input_bam=$1
dict_file=$2
output_bam=$3
threads=$4
tmp_prefix="./temp_$$"

# Ensure samtools is installed
if ! command -v samtools &> /dev/null; then
    echo "Error: samtools could not be found, please install samtools."
    exit 1
fi

echo "Extracting BAM headers..."
samtools view -@ "$threads" -H "$input_bam" > "${tmp_prefix}_header.sam"

echo "Extracting read IDs from BAM..."
samtools view -@ "$threads" "$input_bam" | cut -f1 | sort -u > "${tmp_prefix}_bam_reads.list"

echo "Filtering dictionary for this BAM file..."
gzip -dc "$dict_file" | grep -F -f "${tmp_prefix}_bam_reads.list" | sort --parallel="$threads" -k1,1 > "${tmp_prefix}_filtered_dict.tsv"

echo "Sorting BAM reads..."
samtools view -@ "$threads" "$input_bam" | sort --parallel="$threads" -k1,1 > "${tmp_prefix}_sorted_reads.tsv"

echo "Merging dictionary tags with BAM reads..."
join -t $'\t' -1 1 -2 1 "${tmp_prefix}_sorted_reads.tsv" "${tmp_prefix}_filtered_dict.tsv" > "${tmp_prefix}_tagged_reads.tsv"

echo "Reassembling BAM..."
cat "${tmp_prefix}_header.sam" "${tmp_prefix}_tagged_reads.tsv" > "${tmp_prefix}_final.sam"

# Step 7: Convert back to BAM
echo "Converting to BAM..."
samtools view -@ "$threads" -b -o "${tmp_prefix}_unsorted.bam" "${tmp_prefix}_final.sam"

# Step 8: Sort BAM by coordinate order before indexing
echo "Sorting BAM by coordinate order..."
samtools sort -@ "$threads" -o "$output_bam" "${tmp_prefix}_unsorted.bam"

# Step 9: Index BAM file
echo "Indexing BAM..."
samtools index -@ "$threads" "$output_bam"

# Cleanup temporary files
rm "${tmp_prefix}_header.sam" "${tmp_prefix}_bam_reads.list" "${tmp_prefix}_filtered_dict.tsv" \
   "${tmp_prefix}_sorted_reads.tsv" "${tmp_prefix}_tagged_reads.tsv" "${tmp_prefix}_final.sam" "${tmp_prefix}_unsorted.bam"

# Confirm completion
if [ $? -eq 0 ]; then
    echo "✅ Tags added and BAM file indexed successfully to $output_bam."
else
    echo "❌ An error occurred during the tagging process."
    exit 1
fi
