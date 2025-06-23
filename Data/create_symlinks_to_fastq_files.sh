#!/bin/bash

# -------------------------------
# Overview:
# This script reads a CSV file where each row corresponds to a sample.
# It collects sample IDs from column 2 and uses them to filter FASTQ files found in specified folders (column 5).
# It creates symbolic links to matched FASTQ files in a central output directory, renaming them consistently.
# -------------------------------

# Input CSV file with sample metadata (semicolon-delimited)
csv_file="map_to_fastq_HTL2_168_samples.csv"

# Output directory where symbolic links will be created
output_dir="HTL2_fastq_files"
echo "Running on CSV: $csv_file"

# Ensure the CSV uses Unix-style line endings (in case edited in Windows)
dos2unix "$csv_file" 2>/dev/null

# Create the output directory if it doesn't already exist
mkdir -p "$output_dir"
echo "Output directory: $output_dir"

# -------------------------------
# Step 1: Build list of allowed sample IDs
# -------------------------------

# Declare an associative array to store allowed sample IDs (e.g., HTL214)
declare -A allowed_ids

echo ""
echo "Extracting allowed IDs from column 2..."

# Loop through each line (after the header), splitting on semicolon
# Only column 2 is extracted (the numeric ID), all others are ignored using `_`
while IFS=';' read -r _ id batch _ _; do
  # Skip row if column 2 contains "dont_select"
  [[ "$id" == *dont_select* ]] && continue

  # Remove whitespace from the ID, just in case
  id_trimmed=$(echo "$id" | tr -d '[:space:]')
  batch_trimmed=$(echo "$batch" | tr -d '[:space:]')

   # Combine both columns for a unique key, e.g., "HTL214_180830"
  allowed_id="HTL${id_trimmed}_${batch_trimmed}"
  

  # Store the formatted ID in the associative array
  allowed_ids["$allowed_id"]=1
done < <(tail -n +2 "$csv_file")  # Skip the first (header) line



# -------------------------------
# Step 2: Scan directories listed in column 5 and create symlinks
# -------------------------------

echo ""
echo "Scanning directories..."

# Extract one unique directory path per batch (column 5, using column 3 for uniqueness)
# This avoids scanning the same directory multiple times
awk -F';' 'NR > 1 && !seen[$3]++ { print $5 }' "$csv_file" | while read -r path; do
  echo ""
  echo "Processing directory: $path"

  # Make sure the path is a directory
  if [ -d "$path" ]; then
    # Loop over all .fastq.gz files in the folder
    for file in "$path"/*.fastq.gz; do
      # If no file matches, skip iteration
      [ -e "$file" ] || continue

      # Get just the filename (no path)
      filename=$(basename "$file")
      echo "  Found file: $filename"

      # -------------------------------
      # Step 2a: Normalize filename prefix
      # -------------------------------

      # Extract the prefix (i.e., the first token before the first underscore)
      prefix=$(echo "$filename" | cut -d'_' -f1)

      # If the prefix doesn't start with HTL, add it
      if [[ ! "$prefix" =~ ^HTL ]]; then
        prefix="HTL$prefix"
      fi

      # Normalize: remove any leading zeros after HTL (e.g., HTL00214 → HTL214)
      normalized_prefix=$(echo "$prefix" | sed -E 's/^HTL0*/HTL/')

      echo "    Original prefix: $prefix"
      echo "    Normalized prefix: $normalized_prefix"

        # -------------------------------
      # Step 2b: Check if this sample is in the allowed list
      # -------------------------------

      # Extract the batch name from the path (e.g., last folder in the path)
      batch_from_path=$(basename "$path")

      # Combine normalized prefix and batch to match allowed_ids key
      normalized_id_batch="${normalized_prefix}_${batch_from_path}"
      echo "normalized_id_batch": ${normalized_id_batch}
      

      # +_ is a way to test if an array key exists
      # It does not access the value — only checks presence of the key
      # The _ is just a placeholder and could be replaced with anything non-empty
      if [[ ${allowed_ids["$normalized_id_batch"]+_} ]]; then

        # Replace original prefix in filename with normalized prefix
        symlink_name="$normalized_prefix$(echo "$filename" | sed -E 's/^([^_]+)//')"

        echo "Match found — Creating symlink: $output_dir/$symlink_name"
        ln -s "$file" "$output_dir/$symlink_name"
      else
        echo "Skipped — ${normalized_id_batch} not in allowed list."
      fi

      # Optional: sleep for debugging or rate-limiting (no effect here)
      sleep 0
    done
  else
    echo "Directory not found: $path"
  fi
done

echo ""
echo "Done. Symlinks created in: $output_dir"
