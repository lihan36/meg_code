#!/usr/bin/env bash
set -euo pipefail

# =============================================================
# Generate expression matrix with featureCounts
# =============================================================

# 1. Load user-defined paths from config.sh
source "$(dirname "$0")/../../config/config.sh"

# 2. Define input and output paths
# You can override these by passing arguments to the script
INPUT_GFF3="${1:-${PROJECT_DIR}/02_annotation/braker_out/braker_final.gff3}"
ALIGN_DIR="${2:-${RNASEQ_DIR}/alignment}"
OUTPUT_DIR="${3:-${RNASEQ_DIR}/counts}"
GTF_FILE="${RNASEQ_DIR}/braker_annotation.gtf"
LOG_FILE="${PROJECT_DIR}/logs/featurecounts.log"

# 3. Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"
mkdir -p "$(dirname "${LOG_FILE}")"

# 4. Convert GFF3 to GTF for featureCounts
echo "Converting GFF3 to GTF..." > "${LOG_FILE}"
gffread "${INPUT_GFF3}" -T -o "${GTF_FILE}" >> "${LOG_FILE}" 2>&1

# 5. Run featureCounts for each sample
echo "Starting featureCounts quantification..." >> "${LOG_FILE}"
for bam in "${ALIGN_DIR}"/*_Aligned.sortedByCoord.out.bam; do
    sample=$(basename "$bam" _Aligned.sortedByCoord.out.bam)
    echo "Processing ${sample}..." >> "${LOG_FILE}"
    featureCounts \
        -a "${GTF_FILE}" \
        -o "${OUTPUT_DIR}/${sample}_counts.txt" \
        -t exon \
        -g gene_id \
        -p \
        -T "${THREADS}" \
        -S 0 \
        "$bam" >> "${LOG_FILE}" 2>&1
done

echo "featureCounts quantification completed. Output: ${OUTPUT_DIR}"