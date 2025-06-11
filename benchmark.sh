#!/bin/bash

# Check input
if [[ -z "$1" ]]; then
  echo "Usage: $0 <input_folder>"
  exit 1
fi

input_dir="$1"
output_file="${input_dir%/}/benchmark_summary.tsv"
echo -e "step\tavg_time_min\tmin_time_min\tmax_time_min\tavg_mem_mb\tmin_mem_mb\tmax_mem_mb" > "$output_file"

# Function to convert h:m:s to minutes
hms_to_min() {
    IFS=':' read -r h m s <<< "$1"
    echo "$((10#$h * 60 + 10#$m)) + 10#$s / 60" | bc -l
}

# Loop through each benchmark step folder
for folder in "$input_dir"/*/; do
    step=$(basename "$folder")
    
    runtimes=()
    memories=()

    for file in "$folder"/*.txt; do
        [[ ! -f "$file" ]] && continue

        # Get second line of file
        line=$(awk 'NR==2' "$file")
        time_hms=$(echo "$line" | awk '{print $2}')
        max_rss_kb=$(echo "$line" | awk '{print $3}')  # Adjust column as needed

        [[ -z "$time_hms" || -z "$max_rss_kb" ]] && continue

        # Convert time and memory
        time_min=$(hms_to_min "$time_hms")
        mem_mb=$(echo "$max_rss_kb / 1024" | bc -l)

        runtimes+=("$time_min")
        memories+=("$mem_mb")
    done

    [[ ${#runtimes[@]} -eq 0 ]] && continue

    avg_time=$(printf '%s\n' "${runtimes[@]}" | awk '{sum+=$1} END {printf "%.2f", sum/NR}')
    min_time=$(printf '%.2f\n' "${runtimes[@]}" | sort -n | head -1)
    max_time=$(printf '%.2f\n' "${runtimes[@]}" | sort -n | tail -1)

    avg_mem=$(printf '%s\n' "${memories[@]}" | awk '{sum+=$1} END {printf "%.2f", sum/NR}')
    min_mem=$(printf '%.2f\n' "${memories[@]}" | sort -n | head -1)
    max_mem=$(printf '%.2f\n' "${memories[@]}" | sort -n | tail -1)

    echo -e "${step}\t${avg_time}\t${min_time}\t${max_time}\t${avg_mem}\t${min_mem}\t${max_mem}" >> "$output_file"
done

echo "✅ Benchmark summary saved to: $output_file"
