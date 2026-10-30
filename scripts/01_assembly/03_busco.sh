#!/usr/bin/env bash
set -euo pipefail

# =============================================================
# Run BUSCO evaluation on the assembled genome
# =============================================================

# 1. Load user-defined paths from config.sh
source "$(dirname "$0")/../../config/config.sh"

# 2. Define input and output files
# You can override these by passing arguments to the script
INPUT_GENOME="${1:-${PROJECT_DIR}/01_assembly/flye_assembly/assembly.fasta}"
OUTPUT_DIR="${2:-${PROJECT_DIR}/01_assembly/busco_flye}"
LOG_FILE="${PROJECT_DIR}/logs/busco_flye.log"

# 3. Create output directory if it doesn't exist
mkdir -p "$(dirname "${OUTPUT_DIR}")"
mkdir -p "$(dirname "${LOG_FILE}")"

# 4. Run BUSCO in the background
echo "Starting BUSCO evaluation..."
nohup busco -i "${INPUT_GENOME}" \
    -o "$(basename "${OUTPUT_DIR}")" \
    -l aves_odb10 \
    -m genome \
    -c "${THREADS}" \
    --force \
    > "${LOG_FILE}" 2>&1 &

echo "BUSCO evaluation started in background. Check log: ${LOG_FILE}"