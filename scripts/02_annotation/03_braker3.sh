#!/usr/bin/env bash
set -euo pipefail

# =============================================================
# Run BRAKER3 for gene prediction
# =============================================================

# 1. Load user-defined paths from config.sh
source "$(dirname "$0")/../../config/config.sh"

# 2. Define input and output files
# You can override these by passing arguments to the script
INPUT_GENOME="${1:-${PROJECT_DIR}/01_assembly/racon_polish/racon_polished.fasta}"
INPUT_BAM="${2:-${PROJECT_DIR}/03_rnaseq/alignment/all_samples_merged.bam}"
INPUT_PROTEIN="${3:-${PROJECT_DIR}/00_rawdata/protein.fa}"
WORK_DIR="${4:-${PROJECT_DIR}/02_annotation/braker_out}"
LOG_FILE="${PROJECT_DIR}/logs/braker_run.log"

# 3. Create output directory if it doesn't exist
mkdir -p "$(dirname "${WORK_DIR}")"
mkdir -p "$(dirname "${LOG_FILE}")"

# 4. Export software paths for BRAKER3 (already defined in config.sh)
export AUGUSTUS_CONFIG_PATH
export AUGUSTUS_BIN_PATH
export AUGUSTUS_SCRIPTS_PATH
export GENEMARK_PATH
export PROTHINT_PATH
export BAMTOOLS_PATH
export DIAMOND_PATH
export SAMTOOLS_PATH
export TSEBRA_PATH

# 5. Run BRAKER3 in the background
echo "Starting BRAKER3 gene prediction..." > "${LOG_FILE}"
nohup "${AUGUSTUS_BIN_PATH}/braker.pl" \
    --genome="${INPUT_GENOME}" \
    --bam="${INPUT_BAM}" \
    --prot_seq="${INPUT_PROTEIN}" \
    --threads="${THREADS}" \
    --species="chengxiong" \
    --workingdir="${WORK_DIR}" \
    --gff3 \
    --softmasking \
    --BAMTOOLS_PATH="${BAMTOOLS_PATH}" \
    --DIAMOND_PATH="${DIAMOND_PATH}" \
    --TSEBRA_PATH="${TSEBRA_PATH}" \
    --skip_fixing_broken_genes \
    >> "${LOG_FILE}" 2>&1 &

echo "BRAKER3 started in background. Check log: ${LOG_FILE}"