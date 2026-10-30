#!/usr/bin/env bash
set -euo pipefail

# =============================================================
# Align all RNA-seq samples with STAR
# =============================================================

# 1. Load user-defined paths from config.sh
source "$(dirname "$0")/../../config/config.sh"

# 2. Define input and output paths
SAMPLE_LIST_FILE="${1:-${SAMPLE_LIST}}"
GENOME_DIR="${2:-${RNASEQ_DIR}/genome_index}"
ALIGN_DIR="${3:-${RNASEQ_DIR}/alignment}"
LOG_DIR="${PROJECT_DIR}/logs"
STAR_LOG_FILE="${LOG_DIR}/star_batch.log"

# 3. Create output directories if they don't exist
mkdir -p "${ALIGN_DIR}"
mkdir -p "${LOG_DIR}"

# 4. Start alignment pipeline
echo "=============================================" > "${STAR_LOG_FILE}"
echo "STAR Alignment Pipeline" >> "${STAR_LOG_FILE}"
echo "Start: $(date)" >> "${STAR_LOG_FILE}"
echo "=============================================" >> "${STAR_LOG_FILE}"

PROCESSED=0
FAILED=0

while IFS=$'\t' read -r sample r1 r2; do
    # Skip empty lines
    if [ -z "${sample}" ]; then continue; fi

    echo "" >> "${STAR_LOG_FILE}"
    echo "--- Processing: ${sample} ($(date)) ---" >> "${STAR_LOG_FILE}"

    # Check if already done (skip if BAM exists)
    if [ -f "${ALIGN_DIR}/${sample}_Aligned.sortedByCoord.out.bam" ]; then
        echo "   Already done: ${sample}, skipping." >> "${STAR_LOG_FILE}"
        PROCESSED=$((PROCESSED + 1))
        continue
    fi

    # Check input files
    if [ ! -f "${r1}" ] || [ ! -f "${r2}" ]; then
        echo "   Missing input files for ${sample}: ${r1} or ${r2}" >> "${STAR_LOG_FILE}"
        FAILED=$((FAILED + 1))
        continue
    fi

    # Run STAR
    STAR \
        --genomeDir "${GENOME_DIR}" \
        --readFilesIn "${r1}" "${r2}" \
        --readFilesCommand zcat \
        --outFileNamePrefix "${ALIGN_DIR}/${sample}_" \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMattributes Standard \
        --outSAMstrandField intronMotif \
        --outFilterMultimapNmax 20 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 100000 \
        --alignMatesGapMax 100000 \
        --sjdbScore 1 \
        --limitBAMsortRAM 30000000000 \
        --runThreadN "${THREADS}" \
        >> "${LOG_DIR}/${sample}_star.log" 2>&1

    # Check if successful
    if [ $? -eq 0 ] && [ -f "${ALIGN_DIR}/${sample}_Aligned.sortedByCoord.out.bam" ]; then
        echo "   Success: ${sample}" >> "${STAR_LOG_FILE}"
        samtools index "${ALIGN_DIR}/${sample}_Aligned.sortedByCoord.out.bam"
        PROCESSED=$((PROCESSED + 1))
    else
        echo "   Failed: ${sample}" >> "${STAR_LOG_FILE}"
        FAILED=$((FAILED + 1))
    fi
done < "${SAMPLE_LIST_FILE}"

echo "" >> "${STAR_LOG_FILE}"
echo "=============================================" >> "${STAR_LOG_FILE}"
echo "End: $(date)" >> "${STAR_LOG_FILE}"
echo "Processed: ${PROCESSED}" >> "${STAR_LOG_FILE}"
echo "Failed: ${FAILED}" >> "${STAR_LOG_FILE}"
echo "=============================================" >> "${STAR_LOG_FILE}"

echo "STAR batch alignment completed. Check log: ${STAR_LOG_FILE}"