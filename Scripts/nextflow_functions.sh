# functions.sh
# This file defines reusable Bash functions for Nextflow processes.
# Source it in any process with: `source ${file('functions.sh')}`

################################################################################
# move_and_rename_logs
# --------------------
# Moves log files that start with a digit and end with a specified extension
# (e.g., .log, .err, .out) into a specified log directory, while prefixing
# each filename with a task-specific tag.
#
# Usage:
#   move_and_rename_logs log
#   move_and_rename_logs err
#
# Requires:
#   - LOG_DIR: path to the target log directory (exported before calling)
#   - TAG: a string to prepend to each file (e.g., task.tag or sample ID)
################################################################################
move_named_log() {
  # Default cluster_id to 'unknown'
  local cluster_id="unknown"

  # Find the first file starting with 6+ digits and ending in .log
  local log_file
  log_file=$(ls "$PWD"/[0-9][0-9][0-9][0-9][0-9][0-9]*.log 2>/dev/null)
  echo "[DEBUG] log_file found: ${log_file}" >> "$PWD/.command.out"
  

  if [[ -n "$log_file" ]]; then
    cluster_id=$(basename "$log_file")
    cluster_id="${cluster_id%%.log}"
  fi

  # Create symlink to the log file in the log directory
  ln -s "$log_file" "$LOG_DIR/${cluster_id}_${TAG}_condor.log"

  # Loop over all input patterns
  for pattern in "$@"; do
    for file_to_move in $pattern; do  # Expand globs like "*.out"
      [[ -f "$file_to_move" ]] || {
        echo "[WARN] File not found: $file_to_move" >> "$PWD/.command.err"
        continue
      }
      local base_name=$(basename "$file_to_move")
      local new_name="${cluster_id}_${TAG}_${base_name}"
      cp "$file_to_move" "${LOG_DIR}/${new_name}"
    done
  done
}



################################################################################
# write_summary
# -------------
# Appends a simple status message to a summary log file in the log directory.
#
# Usage:
#   write_summary sample123 DONE
#
# Arguments:
#   \$1 - Sample ID or task name
#   \$2 - Status message (e.g., DONE, FAILED)
#
# Requires:
#   - LOG_DIR: path to the log directory
################################################################################
write_summary() {
  echo "Sample: \$1, Status: \$2" >>"\${LOG_DIR}/summary.txt"
}
