#!/usr/bin/env bash
set -euo pipefail

# =============================================================
# Functional annotation with eggNOG-mapper
# =============================================================

# 1. Load user-defined paths from config.sh
source "$(dirname "$0")/../../config/config.sh"

# 2. Define input and output files
# You can override these by passing arguments to the script
INPUT_PROTEIN="${1:-${PROJECT_DIR}/02_annotation/braker_out/braker_final_protein.fa}"
OUTPUT_DIR="${2:-${PROJECT_DIR}/02_annotation}"
LOG_FILE="${PROJECT_DIR}/logs/eggnog_final_run.log"

# 3. Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"
mkdir -p "$(dirname "${LOG_FILE}")"

# 4. Run eggNOG-mapper in the background
echo "Starting eggNOG-mapper..." > "${LOG_FILE}"
nohup emapper.py \
    -i "${INPUT_PROTEIN}" \
    -o eggnog_output \
    --cpu "${THREADS}" \
    --tax_scope 8782 \
    --dbmem \
    --data_dir "${EGGNOG_DATA_DIR}" \
    --override \
    >> "${LOG_FILE}" 2>&1 &

echo "eggNOG-mapper started in background. Check log: ${LOG_FILE}"