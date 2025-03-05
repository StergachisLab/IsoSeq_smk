#!/bin/bash
# Adriana Sedeno
# Append tags to a BAM file using a gzipped TSV dictionary

# Check for the correct number of arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input.bam> <dictionary.tsv.gz> <output.bam>"
    exit 1
fi

input_bam=$1
dict_file=$2
output_bam=$3
tmp_prefix="./temp_$$"

# Ensure samtools is installed
if ! command -v samtools &> /dev/null; then
    echo "Error: samtools could not be found, please install samtools."
    exit 1
fi

# Convert BAM to SAM and preserve header
samtools view -@ 4 -h "$input_bam" > "${tmp_prefix}_temp.sam"
if [ $? -ne 0 ]; then
    echo "Error: samtools view command failed."
    exit 1
fi

# Prepare the dictionary for use in awk (decompress on the fly)
gzip -dc "$dict_file" | awk 'BEGIN { FS=OFS="\t" }
     { tag = "";
       for (i = 2; i <= NF; i++) tag = tag OFS $i;
       tags[$1] = tag }
     END { for (key in tags) print key, tags[key] }' > "${tmp_prefix}_tags.tmp"
if [ $? -ne 0 ]; then
    echo "Error: preparing dictionary for awk failed."
    exit 1
fi

# Handle the tab-separated dictionary and append tags to matching reads in the SAM file
awk 'BEGIN { FS=OFS="\t" }
     NR==FNR { tags[$1] = substr($0, index($0, $2)); next }
     /^@/ { print; next }
     { if ($1 in tags) print $0 tags[$1]; else print }' "${tmp_prefix}_tags.tmp" "${tmp_prefix}_temp.sam" > "${tmp_prefix}_temp_with_tags.sam"
if [ $? -ne 0 ]; then
    echo "Error: awk command failed."
    exit 1
fi

# Convert the modified SAM back to BAM
samtools view -@ 4 -b -o "$output_bam" "${tmp_prefix}_temp_with_tags.sam"
if [ $? -ne 0 ]; then
    echo "Error: samtools view command failed."
    exit 1
fi

# Index the BAM file
samtools index -@ 4 "$output_bam"
if [ $? -ne 0 ]; then
    echo "Error: samtools index command failed."
    exit 1
fi

# Cleanup temporary files
rm "${tmp_prefix}_temp.sam" "${tmp_prefix}_temp_with_tags.sam" "${tmp_prefix}_tags.tmp"

# Confirm completion
if [ $? -eq 0 ]; then
    echo "Tags added and BAM file indexed successfully to $output_bam."
else
    echo "An error occurred during the tagging process."
    exit 1
fi
