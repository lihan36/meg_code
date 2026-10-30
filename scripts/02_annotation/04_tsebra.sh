#!/usr/bin/env bash
set -euo pipefail

# =============================================================
# Integrate gene predictions with TSEBRA
# =============================================================

# 1. Load user-defined paths from config.sh
source "$(dirname "$0")/../../config/config.sh"

# 2. Define input and output files
# You can override these by passing arguments to the script
BRAKER_OUT_DIR="${1:-${PROJECT_DIR}/02_annotation/braker_out}"
INPUT_GENOME="${2:-${PROJECT_DIR}/01_assembly/racon_polish/racon_polished.fasta}"
TSEBRA_SCRIPT="${TSEBRA_PATH}/tsebra.py"
TSEBRA_CONFIG="${BRAKER_OUT_DIR}/tsebra_config.cfg"
LOG_FILE="${PROJECT_DIR}/logs/tsebra_medium.log"

# 3. Create output directory if it doesn't exist
mkdir -p "${BRAKER_OUT_DIR}"
mkdir -p "$(dirname "${LOG_FILE}")"

# 4. Run TSEBRA
echo "Starting TSEBRA integration..." > "${LOG_FILE}"
if [ ! -f "${TSEBRA_CONFIG}" ]; then
    echo "WARNING: tsebra_config.cfg not found in ${BRAKER_OUT_DIR}. Please make sure it exists!" >> "${LOG_FILE}"
fi

python3 "${TSEBRA_SCRIPT}" \
    -g "${BRAKER_OUT_DIR}/braker.gtf" \
    -e "${BRAKER_OUT_DIR}/hintsfile.gff" \
    -c "${TSEBRA_CONFIG}" \
    -o "${BRAKER_OUT_DIR}/braker_tsebra_medium.gtf" \
    -filter_single_exon_genes \
    -v 1 >> "${LOG_FILE}" 2>&1

# 5. Integrate final version and convert to GFF3
echo "--- Copying final files ---" >> "${LOG_FILE}"
cp "${BRAKER_OUT_DIR}/braker_tsebra_medium.gtf" "${BRAKER_OUT_DIR}/braker_final.gtf"
cp "${BRAKER_OUT_DIR}/braker_tsebra_medium_protein.fa" "${BRAKER_OUT_DIR}/braker_final_protein.fa"

echo "--- Converting final GTF to GFF3 ---" >> "${LOG_FILE}"
gffread "${BRAKER_OUT_DIR}/braker_final.gtf" \
    -g "${INPUT_GENOME}" \
    -o "${BRAKER_OUT_DIR}/braker_final.gff3" >> "${LOG_FILE}" 2>&1

echo "TSEBRA integration completed. Final GFF3: ${BRAKER_OUT_DIR}/braker_final.gff3"