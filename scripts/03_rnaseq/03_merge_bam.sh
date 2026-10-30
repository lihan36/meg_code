#!/usr/bin/env bash
set -euo pipefail

# =============================================================
# Merge all sample BAM files into a single BAM
# =============================================================

# 1. Load user-defined paths from config.sh
source "$(dirname "$0")/../../config/config.sh"

# 2. Define input and output paths
# You can override these by passing arguments to the script
ALIGN_DIR="${1:-${RNASEQ_DIR}/alignment}"
OUTPUT_BAM="${2:-${RNASEQ_DIR}/alignment/all_samples_merged.bam}"
LOG_FILE="${PROJECT_DIR}/logs/merge_bam.log"

# 3. Create log directory if it doesn't exist
mkdir -p "$(dirname "${LOG_FILE}")"

# 4. Merge all BAM files
echo "Starting BAM file merge..." > "${LOG_FILE}"
samtools merge -@ "${THREADS}" \
    "${OUTPUT_BAM}" \
    "${ALIGN_DIR}"/*_Aligned.sortedByCoord.out.bam \
    >> "${LOG_FILE}" 2>&1

# 5. Index the merged BAM
echo "Indexing merged BAM..." >> "${LOG_FILE}"
samtools index "${OUTPUT_BAM}" >> "${LOG_FILE}" 2>&1

echo "BAM merge completed. Output: ${OUTPUT_BAM}"