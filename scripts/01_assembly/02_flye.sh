#!/usr/bin/env bash
set -euo pipefail

# =============================================================
# Convert BAM to FASTQ format
# =============================================================

# 1. Load user-defined paths from config.sh
source "$(dirname "$0")/../../config/config.sh"

# 2. Define input and output files using config variables
# You can override these by passing arguments to the script
INPUT_BAM="${1:-${RAW_DATA_DIR}/P1DLD1705050-1-P2_1.bam}"
OUTPUT_FASTQ="${2:-${PROJECT_DIR}/00_rawdata/pacbio_raw.fastq}"
LOG_FILE="${PROJECT_DIR}/logs/flye_assembly.log"

# 3. Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"
mkdir -p "$(dirname "${LOG_FILE}")"

# 4. Run Flye assembly in the background
echo "Starting Flye assembly..."
nohup flye --pacbio-raw "${INPUT_FASTQ}" \
    --out-dir "${OUTPUT_DIR}" \
    --threads "${THREADS}" \
    --genome-size "${GENOME_SIZE}" \
    --iterations 3 \
    --min-overlap 5000 \
    > "${LOG_FILE}" 2>&1 &

echo "Flye assembly started in background. Check log: ${LOG_FILE}"