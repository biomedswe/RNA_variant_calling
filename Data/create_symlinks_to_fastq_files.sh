#!/bin/bash

# Description:
# For each unique batch in the CSV, this script finds all files inside the corresponding directory
# (from the 'Full path in server' column) and creates symlinks to them in the subfolder HTL2_fastq_files,
# using the same filenames.

csv_file="map_to_fastq_HTL2_168_samples.csv"
output_dir="HTL2_fastq_files"
echo "Running: $csv_file"

# Ensure file has Unix line endings
dos2unix "$csv_file" 2>/dev/null

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Extract unique folder paths based on first occurrence of each batch
awk -F';' '                        # Use semicolon as field delimiter (since CSV uses ; not ,)

  NR > 1 && !seen[$3]++ {          # For all rows *except* the header (NR > 1):
                                   #   - $3 refers to the "batch" column
                                   #   - seen[$3] is an associative array that keeps track of how many times each batch value appears
                                   #   - !seen[$3]++ means: only evaluate true the *first time* we see a given batch
                                   #     (it returns true when seen[$3] == 0, then increments the count)

    print $5;                      # Print the "Full path in server" column (field 5), once per unique batch
  }

' "$csv_file" | while read -r path; do # Pipe the paths into the while-loop

  echo "Processing directory: $path"

  # Check if directory exists
  if [ -d "$path" ]; then
    # Loop over each .fastq.gz file inside the path
    for file in "$path"/*.fastq.gz; do
      # Get the base filename (without path)
      filename=$(basename "$file")
      echo "Found file: $filename"
      

      # Check if the filename starts with "HTL"
      # =~ is used for regex matching in bash
      if [[ "$filename" =~ ^HTL ]]; then
        continue  # Skip files that already start with "HTL"
      elif [[ "$filename" =~ ^[0-9] ]]; then
        filename="HTL$filename"
      else
        echo "Skipping file (invalid prefix): $filename"
        continue
      fi

      # Print what symlink is being created
      echo "Creating symlink: $output_dir/$filename -> $file"

      # Create the symbolic link in the output directory
      ln -s "$file" "$output_dir/$filename"
    done
  else
    # Print a warning if the directory doesn't exist
    echo "Directory not found: $path"
  fi

done
# End of script
echo "Symlinks created in directory: $output_dir"
