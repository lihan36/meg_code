#!/usr/bin/env bash
set -euo pipefail

# =============================================================
# Generate STAR genome index
# =============================================================

# 1. Load user-defined paths from config.sh
source "$(dirname "$0")/../../config/config.sh"

# 2. Define input and output paths
# You can override these by passing arguments to the script
INPUT_GENOME="${1:-${PROJECT_DIR}/02_annotation/repeat_mask_combined/genome.fasta.masked}"
OUTPUT_DIR="${2:-${RNASEQ_DIR}/genome_index}"
LOG_FILE="${PROJECT_DIR}/logs/genome_index.log"

# 3. Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"
mkdir -p "$(dirname "${LOG_FILE}")"

# 4. Run STAR genomeGenerate in the background
echo "Starting STAR genome index generation..." > "${LOG_FILE}"
nohup STAR \
    --runMode genomeGenerate \
    --genomeDir "${OUTPUT_DIR}" \
    --genomeFastaFiles "${INPUT_GENOME}" \
    --runThreadN "${THREADS}" \
    --genomeSAindexNbases 13 \
    --genomeChrBinNbits 16 \
    --limitGenomeGenerateRAM 300000000000 \
    >> "${LOG_FILE}" 2>&1 &

echo "STAR index generation started in background. Check log: ${LOG_FILE}"