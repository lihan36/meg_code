#!/usr/bin/env bash
set -euo pipefail

# =============================================================
# Filter raw PacBio reads with Filtlong
# =============================================================

# 1. Load user-defined paths from config.sh
source "$(dirname "$0")/../../config/config.sh"

# 2. Define input and output files
# You can override these by passing arguments to the script
INPUT_FASTQ="${1:-${PROJECT_DIR}/00_rawdata/pacbio_raw.fastq}"
OUTPUT_FASTQ="${2:-${PROJECT_DIR}/01_assembly/pacbio_filtered.fastq.gz}"

# 3. Run Filtlong and compress output
echo "Filtering ${INPUT_FASTQ} with Filtlong..."
filtlong --min_length 1000 --keep_percent 95 "${INPUT_FASTQ}" | gzip > "${OUTPUT_FASTQ}"
echo "Filtering completed. Output: ${OUTPUT_FASTQ}"