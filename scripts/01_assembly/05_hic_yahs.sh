#!/usr/bin/env bash
set -euo pipefail

# =============================================================
# Hi-C scaffolding with BWA, SAMtools and YaHS
# =============================================================

# 1. Load user-defined paths from config.sh
source "$(dirname "$0")/../../config/config.sh"

# 2. Define input and output files
# You can override these by passing arguments to the script
INPUT_GENOME="${1:-${PROJECT_DIR}/01_assembly/racon_polish/racon_polished.fasta}"
HIC_R1="${2:-${HIC_DATA_DIR}/L1EJE1700081-GG002.R1.raw.fastq.gz}"
HIC_R2="${3:-${HIC_DATA_DIR}/L1EJE1700081-GG002.R2.raw.fastq.gz}"
OUTPUT_DIR="${4:-${PROJECT_DIR}/01_assembly/hic_scaffolding}"
LOG_FILE="${PROJECT_DIR}/logs/hic_yahs.log"

# 3. Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"
mkdir -p "$(dirname "${LOG_FILE}")"

# 4. Clean contig names and prepare files
echo "Starting Hi-C scaffolding..." > "${LOG_FILE}"
cd "${OUTPUT_DIR}"

echo "--- Cleaning contig names ---" >> "${LOG_FILE}"
awk '/^>/{print ">Contig_" ++i; next}{print}' "${INPUT_GENOME}" > racon_polished.clean.fasta

echo "--- Building BWA index and faidx ---" >> "${LOG_FILE}"
bwa index racon_polished.clean.fasta
samtools faidx racon_polished.clean.fasta

# 5. Align Hi-C reads
echo "--- Aligning Hi-C reads ---" >> "${LOG_FILE}"
bwa mem -t "${THREADS}" racon_polished.clean.fasta "${HIC_R1}" "${HIC_R2}" 2>> "${LOG_FILE}" | \
    samtools view -bS -@ "${THREADS}" -o hic_aligned.clean.bam -

# 6. Sort, index, and filter BAM (keep MAPQ >= 10, unique alignments)
echo "--- Sorting and filtering BAM ---" >> "${LOG_FILE}"
samtools sort -@ "${THREADS}" -m 4G -o hic_sorted.clean.bam hic_aligned.clean.bam
samtools index hic_sorted.clean.bam

samtools view -b -q 10 -F 0x100 -@ "${THREADS}" hic_sorted.clean.bam > hic_filtered.clean.bam
samtools index hic_filtered.clean.bam

samtools flagstat hic_filtered.clean.bam > "${PROJECT_DIR}/logs/flagstat_clean.txt"

# 7. Run YaHS scaffolding in background
echo "--- Running YaHS scaffolding ---" >> "${LOG_FILE}"
nohup yahs -o chengxiong_yahs racon_polished.clean.fasta hic_filtered.clean.bam \
    >> "${LOG_FILE}" 2>&1 &

echo "Hi-C scaffolding started in background. Check log: ${LOG_FILE}"
echo "Note: After YaHS finishes, run: agp_to_fasta chengxiong_yahs_scaffolds_final.agp racon_polished.clean.fasta > yahs_scaffolds.clean.fasta"