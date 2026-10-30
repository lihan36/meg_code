# =============================================================
# KEGG Enrichment for split DEGs (Female/Male, Up/Down)
# =============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
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
OUTPUT_DIR <- file.path(WORK_DIR, "KEGG_split_results")
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# =============================================================
# 1. Read gene2pathway to build background (keep only map prefix)
# =============================================================
# Note: Update this path to your actual gene2pathway.txt
pw_file <- file.path(WORK_DIR, "gene2pathway.txt")
if (!file.exists(pw_file)) stop("ERROR: File not found: ", pw_file)

pw_raw <- readLines(pw_file)
pw_raw <- pw_raw[!grepl("^#", pw_raw) & pw_raw != ""]

pw_list <- lapply(pw_raw, function(line) {
  p <- strsplit(line, "\t")[[1]]
  if (length(p) < 2) return(NULL)
  g <- gsub("^anno1\\.", "", p[1])
  g <- gsub("\\..*$", "", g) # Convert to short ID
  paths <- p[-1]
  paths <- paths[grepl("^map", paths)] # Keep only map
  paths <- paths[paths != "-" & paths != ""]
  if (length(paths) == 0) return(NULL)
  data.frame(gene = g, pathway = unique(paths), stringsAsFactors = FALSE)
})

gene2pathway <- unique(do.call(rbind, pw_list))
background_genes <- unique(gene2pathway$gene)
cat("Background genes with pathway annotation:", length(background_genes), "\n")

# =============================================================
# 2. Read Excel and loop enrichment
# =============================================================
excel_file <- file.path(WORK_DIR, "雌雄上调下调DEGs g+log.xlsx")
if (!file.exists(excel_file)) stop("ERROR: File not found: ", excel_file)

sheets <- c("female up", "female down", "male up", "male down")
sheet_names <- c("Female_Up", "Female_Down", "Male_Up", "Male_Down")

for (i in 1:length(sheets)) {
  sheet <- sheets[i]
  name <- sheet_names[i]
  cat("\n========== Analyzing:", name, "==========\n")

  df <- read_excel(excel_file, sheet = sheet)
  deps <- unique(df$id)

  deg_in_pw <- intersect(deps, background_genes)
  cat("Total genes in list:", length(deps), ", with pathway annotation:", length(deg_in_pw), "\n")

  if (length(deg_in_pw) < 3) {
    cat("Too few genes, skipping enrichment.\n")
    next
  }

  # Run enrichment
  res <- enricher(
    gene = deg_in_pw,
    universe = background_genes,
    TERM2GENE = gene2pathway[, c("pathway", "gene")],
    pvalueCutoff = 0.1,
    pAdjustMethod = "BH",
    minGSSize = 3,
    maxGSSize = 500,
    qvalueCutoff = 0.3
  )

  res_df <- as.data.frame(res)

  if (nrow(res_df) == 0) {
    cat("No significant pathways found.\n")
  } else {
    cat("Found", nrow(res_df), "significant/candidate pathways!\n")

    # Save table
    write.table(res_df, file.path(OUTPUT_DIR, paste0("KEGG_", name, ".txt")),
                sep = "\t", quote = FALSE, row.names = FALSE)

    # Dotplot
    n_show <- min(15, nrow(res_df))
    p1 <- dotplot(res, showCategory = n_show) +
      ggtitle(paste0("KEGG Pathway: ", name)) +
      theme_bw(base_size = 11)
    ggsave(file.path(OUTPUT_DIR, paste0("KEGG_", name, "_dotplot.pdf")), p1, width = 10, height = 6)

    # Barplot
    res_df <- res_df[order(res_df$Count, decreasing = FALSE), ]
    res_df$Description <- factor(res_df$Description, levels = res_df$Description)
    p2 <- ggplot(res_df, aes(x = Count, y = Description, fill = p.adjust)) +
      geom_bar(stat = "identity") +
      scale_fill_gradient(low = "red", high = "blue") +
      labs(title = paste0("KEGG Pathway: ", name), x = "Gene Count", y = "") +
      theme_bw(base_size = 11)
    ggsave(file.path(OUTPUT_DIR, paste0("KEGG_", name, "_barplot.pdf")), p2, width = 10, height = 6)

    cat("  Generated plots and tables: KEGG_", name, "\n")
  }
}
cat("\nAll done! Check", OUTPUT_DIR, "for results.\n")