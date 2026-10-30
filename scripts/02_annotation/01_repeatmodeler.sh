#!/usr/bin/env bash
set -euo pipefail

# =============================================================
# Build RepeatModeler database and run RepeatModeler
# =============================================================

# 1. Load user-defined paths from config.sh
source "$(dirname "$0")/../../config/config.sh"

# 2. Define input and output files
# You can override these by passing arguments to the script
INPUT_GENOME="${1:-${PROJECT_DIR}/01_assembly/racon_polish/racon_polished.fasta}"
WORK_DIR="${2:-${PROJECT_DIR}/02_annotation}"
LOG_FILE="${PROJECT_DIR}/logs/repeatmodeler_run.log"

# 3. Create output directory if it doesn't exist
mkdir -p "${WORK_DIR}"
mkdir -p "$(dirname "${LOG_FILE}")"

# 4. Run BuildDatabase and RepeatModeler in the background
cd "${WORK_DIR}"

echo "Starting RepeatModeler pipeline..." > "${LOG_FILE}"

# Build database
echo "--- Building RepeatModeler database ---" >> "${LOG_FILE}"
BuildDatabase -name chengxiong_genome -engine ncbi "${INPUT_GENOME}" >> "${LOG_FILE}" 2>&1

# Run RepeatModeler
echo "--- Running RepeatModeler ---" >> "${LOG_FILE}"
nohup RepeatModeler \
    -database chengxiong_genome \
    -threads "${THREADS}" \
    -LTRStruct \
    >> "${LOG_FILE}" 2>&1 &

echo "RepeatModeler started in background. Check log: ${LOG_FILE}"