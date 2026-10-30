# =============================================================
# GSEA Custom Analysis (Based on full gene Excel inputs)
# =============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(fgsea)
  library(data.table)
})

# =============================================================
# 0. Automatically locate project root
# =============================================================
find_project_root <- function(start_dir = getwd()) {
  dir <- normalizePath(start_dir)
  for (i in 1:10) {
    if (dir.exists(file.path(dir, "config")) &&
        file.exists(file.path(dir, "config", "config.sh"))) {
      return(dir)
    }
    parent <- dirname(dir)
    if (parent == dir) break
    dir <- parent
  }
  return(getwd())
}

PROJECT_ROOT <- find_project_root()
WORK_DIR <- file.path(PROJECT_ROOT, "scripts", "05_enrichment")
OUTPUT_DIR <- file.path(WORK_DIR, "GSEA_results_custom")

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
for (s in c("female", "male")) {
  d <- file.path(OUTPUT_DIR, s)
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# =============================================================
# 1. Input files
# =============================================================
FEMALE_RES <- file.path(WORK_DIR, "H-F_vs_L-F_full_results.xlsx")
MALE_RES   <- file.path(WORK_DIR, "H-M_vs_L-M_full_results.xlsx")
GENE2PATHWAY <- file.path(WORK_DIR, "gene2pathway.xlsx")

for (f in c(FEMALE_RES, MALE_RES, GENE2PATHWAY)) {
  if (!file.exists(f)) stop("ERROR: File not found: ", f)
}

# =============================================================
# 2. Read full gene results and build ranked vectors
# =============================================================
read_full_res <- function(file, label) {
  df <- read_excel(file, sheet = 1)
  colnames(df)[1] <- "gene_id"
  df <- df %>%
    select(gene_id, log2FoldChange, padj) %>%
    filter(!is.na(gene_id), !is.na(log2FoldChange), is.finite(log2FoldChange))
  cat(label, "genes loaded:", nrow(df), "\n")
  df
}

female_res <- read_full_res(FEMALE_RES, "Female")
male_res   <- read_full_res(MALE_RES, "Male")

make_genelist <- function(df) {
  gl <- df$log2FoldChange
  names(gl) <- as.character(df$gene_id)
  gl <- gl[!duplicated(names(gl))]
  sort(gl, decreasing = TRUE)
}

female_geneList <- make_genelist(female_res)
male_geneList   <- make_genelist(male_res)

cat("Female ranked vector length:", length(female_geneList), "\n")
cat("Male ranked vector length:", length(male_geneList), "\n")

# =============================================================
# 3. Read gene2pathway and build TERM2GENE
# =============================================================
g2p_raw <- read_excel(GENE2PATHWAY, sheet = 1)
colnames(g2p_raw)[1:2] <- c("gene", "pathway")

g2p_long <- g2p_raw %>%
  select(gene, pathway) %>%
  filter(!is.na(gene), !is.na(pathway), pathway != "") %>%
  separate_rows(pathway, sep = "[;,|\\s]+") %>%
  filter(pathway != "") %>%
  mutate(
    gene    = str_trim(as.character(gene)),
    pathway = str_trim(as.character(pathway))
  ) %>%
  distinct()

term2gene <- g2p_long %>% select(term = pathway, gene)

cat("TERM2GENE rows:", nrow(term2gene), "\n")
cat("Pathways:", length(unique(term2gene$term)), "\n")
cat("Genes:", length(unique(term2gene$gene)), "\n")

cat("Female rnk overlap with pathway genes:",
    length(intersect(names(female_geneList), unique(term2gene$gene))), "\n")
cat("Male rnk overlap with pathway genes:",
    length(intersect(names(male_geneList),   unique(term2gene$gene))), "\n")

# =============================================================
# 4. Define GSEA runner function
# =============================================================
run_gsea_custom <- function(geneList, term2gene, label, out_subdir,
                            pvalueCutoff = 0.05, pAdjustMethod = "BH",
                            minGSSize = 10, maxGSSize = 500,
                            nPermSimple = 10000, seed = 123) {
  set.seed(seed)
  res <- GSEA(
    geneList      = geneList,
    TERM2GENE     = term2gene,
    minGSSize     = minGSSize,
    maxGSSize     = maxGSSize,
    pvalueCutoff  = pvalueCutoff,
    pAdjustMethod = pAdjustMethod,
    verbose       = TRUE,
    nPermSimple   = nPermSimple,
    seed          = TRUE
  )
  df <- as.data.frame(res)
  cat(label, "significant pathways:", nrow(df), "\n")

  write.table(df,
              file.path(out_subdir, paste0("GSEA_", label, "_results.txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  if (nrow(df) > 0) {
    core_genes <- unique(unlist(strsplit(df$core_enrichment, "/")))
    cat(label, "core_enrichment genes:", length(core_genes),
        " Not in rnk:", length(setdiff(core_genes, names(geneList))), "\n")
  }

  if (nrow(df) > 0) {
    p_dot <- dotplot(res, showCategory = 20, split = ".sign") +
      facet_grid(. ~ .sign) +
      ggtitle(paste0(label, " KEGG-GSEA"))
    ggsave(file.path(out_subdir, paste0("GSEA_", label, "_dotplot.pdf")),
           p_dot, width = 10, height = 8)

    top_n <- min(5, nrow(df))
    p_curve <- gseaplot2(res, geneSetID = 1:top_n,
                         title = paste0(label, " top ", top_n, " pathways"))
    ggsave(file.path(out_subdir, paste0("GSEA_", label, "_top", top_n, "_curves.pdf")),
           p_curve, width = 10, height = 8)
  }
  res
}

# =============================================================
# 5. Run GSEA for Female and Male
# =============================================================
cat("\n========== Female GSEA ==========\n")
gsea_female <- run_gsea_custom(
  geneList   = female_geneList,
  term2gene  = term2gene,
  label      = "female",
  out_subdir = file.path(OUTPUT_DIR, "female")
)

cat("\n========== Male GSEA ==========\n")
gsea_male <- run_gsea_custom(
  geneList   = male_geneList,
  term2gene  = term2gene,
  label      = "male",
  out_subdir = file.path(OUTPUT_DIR, "male")
)

# =============================================================
# 6. Common significant pathways
# =============================================================
f_df <- as.data.frame(gsea_female)
m_df <- as.data.frame(gsea_male)

common_terms <- intersect(f_df$ID, m_df$ID)
cat("\nCommon significant pathways:", length(common_terms), "\n")

if (length(common_terms) > 0) {
  write.table(common_terms,
              file.path(OUTPUT_DIR, "common_significant_pathways.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

cat("\nAll done. Results saved in:", OUTPUT_DIR, "\n")