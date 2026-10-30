#!/usr/bin/env bash
set -euo pipefail

# =============================================================
# Merge de novo library with RepBase and run RepeatMasker
# =============================================================

# 1. Load user-defined paths from config.sh
source "$(dirname "$0")/../../config/config.sh"

# 2. Define input and output files
# You can override these by passing arguments to the script
INPUT_GENOME="${1:-${PROJECT_DIR}/01_assembly/racon_polish/racon_polished.fasta}"
WORK_DIR="${2:-${PROJECT_DIR}/02_annotation}"

DENOVO_LIB="${WORK_DIR}/repeatmodeler_run/chengxiong_genome-families.fa"
REPBASE_EMBL="${WORK_DIR}/RMRBSeqs.embl"
COMBINED_LIB="${WORK_DIR}/combined_library_clean.fa"
OUTPUT_DIR="${WORK_DIR}/repeat_mask_combined"
LOG_FILE="${PROJECT_DIR}/logs/repeatmasker_combined.log"

# 3. Create output directory if it doesn't exist
mkdir -p "${WORK_DIR}"
mkdir -p "${OUTPUT_DIR}"
mkdir -p "$(dirname "${LOG_FILE}")"

# 4. Prepare combined library
cd "${WORK_DIR}"
echo "Starting RepeatMasker pipeline..." > "${LOG_FILE}"

# Extract FASTA from EMBL (RepBase)
echo "--- Extracting FASTA from EMBL ---" >> "${LOG_FILE}"
awk 'BEGIN { seq=""; id=""}
  /^ID/ { if (seq != "") { print ">" id; print seq; seq="" } id=$2 }
  /^SQ/ { next }
  /^\/\// { next }
  /^[^A-Z]/ { gsub(/ /,"", $0); seq = seq $0 }
  END { if (seq != "") { print ">" id; print seq } }' "${REPBASE_EMBL}" > repbase_library.fasta

# Merge de novo + RepBase
echo "--- Merging de novo and RepBase libraries ---" >> "${LOG_FILE}"
cat "${DENOVO_LIB}" repbase_library.fasta > combined_library.fa

# Remove duplicates (identical sequences)
echo "--- Removing duplicates ---" >> "${LOG_FILE}"
seqkit sort -r combined_library.fa | seqkit rmdup -s -o "${COMBINED_LIB}" >> "${LOG_FILE}" 2>&1

# 5. Run RepeatMasker in the background
echo "--- Running RepeatMasker ---" >> "${LOG_FILE}"
nohup RepeatMasker "${INPUT_GENOME}" \
    -lib "${COMBINED_LIB}" \
    -pa "${THREADS}" \
    -dir "${OUTPUT_DIR}" \
    -gff \
    -xsmall \
    -html \
    >> "${LOG_FILE}" 2>&1 &

echo "RepeatMasker started in background. Check log: ${LOG_FILE}"