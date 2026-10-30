# Project: Genome Assembly and Transcriptomic Analysis Pipeline

## 1. Description
This repository contains a comprehensive bioinformatics pipeline for genome assembly, structural and functional annotation, and transcriptomic analysis. The workflow starts from raw sequencing data (PacBio and Hi-C), proceeds through genome assembly and polishing, gene prediction and annotation, RNA-seq quantification and differential expression analysis (DESeq2), and finally performs functional enrichment (KEGG and GSEA) to investigate the molecular mechanisms underlying phenotypic traits.

## 2. System Requirements
- **Operating System**: Linux
- **Software**: samtools, filtlong, flye, busco, minimap2, racon, bwa, yahs, RepeatModeler, RepeatMasker, BRAKER3, TSEBRA, eggNOG-mapper, STAR, featureCounts.
- **R Environment**: R (version >= 4.0 recommended). Required R packages: `DESeq2`, `ggplot2`, `ggrepel`, `readxl`, `dplyr`, `tidyr`, `stringr`, `clusterProfiler`, `enrichplot`, `fgsea`, `data.table`.
- All user-defined paths and parameters are centrally managed in `config/config.sh`.

## 3. Workflow
The pipeline is organized sequentially into the following directories and scripts:

### 3.1 Genome Assembly
- `00_rawdata/01_bam_to_fastq.sh`: Converts raw BAM files to FASTQ format.
- `01_assembly/01_filtlong.sh`: Filters raw PacBio reads to retain high-quality long reads.
- `01_assembly/02_flye.sh`: Performs de novo genome assembly using Flye.
- `01_assembly/03_busco.sh`: Evaluates assembly completeness using BUSCO (aves_odb10).
- `01_assembly/04_racon_polish.sh`: Improves assembly accuracy through 3 rounds of polishing using minimap2 and Racon.
- `01_assembly/05_hic_yahs.sh`: Scaffolds the assembly using Hi-C data via BWA and YaHS.

### 3.2 Genome Annotation
- `02_annotation/01_repeatmodeler.sh`: Builds the RepeatModeler database and identifies de novo repeats.
- `02_annotation/02_repeatmasker.sh`: Runs RepeatMasker using a combined library (de novo + RepBase) to mask repetitive elements.
- `02_annotation/03_braker3.sh`: Predicts protein-coding genes using BRAKER3 (integrating RNA-seq and protein homology).
- `02_annotation/04_tsebra.sh`: Integrates BRAKER3 predictions using TSEBRA and converts the final output to GFF3 format.
- `02_annotation/05_eggnog.sh`: Performs functional annotation of predicted proteins using eggNOG-mapper.

### 3.3 RNA-seq Analysis
- `03_rnaseq/01_star_index.sh`: Generates the STAR genome index using the soft-masked genome.
- `03_rnaseq/02_star_align_all.sh`: Aligns all RNA-seq samples to the genome and sorts the output BAM files.
- `03_rnaseq/03_merge_bam.sh`: Merges all sample BAM files into a single BAM file and indexes it.
- `03_rnaseq/04_featurecounts.sh`: Quantifies gene expression by generating count matrices using featureCounts.
- `03_rnaseq/05_merge_counts.R`: Merges individual count files into a single expression matrix for downstream analysis.

### 3.4 Differential Expression & PCA
- `04_deseq2/01_pca.R`: Performs Principal Component Analysis (PCA) on variance-stabilized transcriptomic data to assess sample clustering.
- `04_deseq2/02_deseq2_deg.R`: Identifies differentially expressed genes (DEGs) between high and low altitude groups using DESeq2, and outputs Venn diagrams and volcano plots.

### 3.5 Functional Enrichment
- `05_enrichment/01_kegg_split_deg.R`: Performs KEGG pathway enrichment analysis for split DEG lists (Female/Male, Up/Down regulated) using hypergeometric tests.
- `05_enrichment/02_gsea_custom.R`: Conducts custom Gene Set Enrichment Analysis (GSEA) based on full gene lists using predefined pathway gene sets.