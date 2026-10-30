#!/usr/bin/env bash
set -euo pipefail

# =============================================================
# Polish assembly with minimap2 and racon (3 rounds)
# =============================================================

# 1. Load user-defined paths from config.sh
source "$(dirname "$0")/../../config/config.sh"

# 2. Define input and output files
# You can override these by passing arguments to the script
INPUT_GENOME="${1:-${PROJECT_DIR}/01_assembly/flye_assembly/assembly.fasta}"
INPUT_FASTQ="${2:-${PROJECT_DIR}/00_rawdata/pacbio_filtered.fastq.gz}"
OUTPUT_DIR="${3:-${PROJECT_DIR}/01_assembly/racon_polish}"
LOG_FILE="${PROJECT_DIR}/logs/racon_polish.log"

# 3. Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"
mkdir -p "$(dirname "${LOG_FILE}")"

# 4. Run polishing in 3 rounds
echo "Starting Racon polishing (3 rounds)..." > "${LOG_FILE}"
cd "${OUTPUT_DIR}"

# Round 1
echo "=== Round 1/3: $(date) ===" >> "${LOG_FILE}"
minimap2 -t "${THREADS}" -x map-pb "${INPUT_GENOME}" "${INPUT_FASTQ}" > reads_1.paf
racon -t "${THREADS}" "${INPUT_FASTQ}" reads_1.paf "${INPUT_GENOME}" > racon_1.fasta

# Round 2
echo "=== Round 2/3: $(date) ===" >> "${LOG_FILE}"
minimap2 -t "${THREADS}" -x map-pb racon_1.fasta "${INPUT_FASTQ}" > reads_2.paf
racon -t "${THREADS}" "${INPUT_FASTQ}" reads_2.paf racon_1.fasta > racon_2.fasta

# Round 3
echo "=== Round 3/3: $(date) ===" >> "${LOG_FILE}"
minimap2 -t "${THREADS}" -x map-pb racon_2.fasta "${INPUT_FASTQ}" > reads_3.paf
racon -t "${THREADS}" "${INPUT_FASTQ}" reads_3.paf racon_2.fasta > racon_polished.fasta

echo "Polishing finished: $(date)" >> "${LOG_FILE}"
echo "Final polished genome: ${OUTPUT_DIR}/racon_polished.fasta" >> "${LOG_FILE}"
echo "Polishing completed. Check log: ${LOG_FILE}"