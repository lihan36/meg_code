#!/usr/bin/env bash
# =============================================================
# Convert BAM to FASTQ format
# =============================================================

# 1. Load user-defined paths from config.sh
source "$(dirname "$0")/../../config/config.sh"
set -euo pipefail
# 2. Define input and output files using config variables
# You can override these by passing arguments to the script
INPUT_BAM="${1:-${RAW_DATA_DIR}/P1DLD1705050-1-P2_1.bam}"
OUTPUT_FASTQ="${2:-${PROJECT_DIR}/00_rawdata/pacbio_raw.fastq}"

# 3. Create output directory if it doesn't exist
mkdir -p "$(dirname "${OUTPUT_FASTQ}")"

# 4. Run samtools fastq
echo "Converting ${INPUT_BAM} to ${OUTPUT_FASTQ}..."
samtools fastq "${INPUT_BAM}" > "${OUTPUT_FASTQ}"
echo "Conversion completed."