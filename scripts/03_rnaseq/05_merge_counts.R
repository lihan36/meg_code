# =============================================================
# Merge featureCounts output into a single expression matrix
# =============================================================

# Note: Run this script from within the 03_rnaseq directory:
#   cd scripts/03_rnaseq
#   Rscript 05_merge_counts.R

# 1. Find all count files
files <- list.files("counts", pattern = "*_counts.txt", full.names = TRUE)
sample_names <- gsub("_counts.txt", "", basename(files))

cat("Number of samples:", length(sample_names), "\n")
cat("Sample list:\n")
print(sample_names)

# 2. Read the first file to get gene IDs
first_df <- read.table(files[1], header = TRUE, sep = "\t", skip = 1)
gene_ids <- first_df$Geneid

# 3. Initialize merged matrix
merged_counts <- data.frame(Geneid = gene_ids, stringsAsFactors = FALSE)

# 4. Loop through all files and merge
for (i in seq_along(files)) {
  df <- read.table(files[i], header = TRUE, sep = "\t", skip = 1)
  idx <- match(gene_ids, df$Geneid)
  merged_counts[[sample_names[i]]] <- df[idx, ncol(df)]
  cat("Loaded:", sample_names[i], "\n")
}

# 5. Save final expression matrix
write.table(merged_counts, "expression_matrix_24samples.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("========================================\n")
cat("Merge completed!\n")
cat("Samples:", ncol(merged_counts) - 1, "\n")
cat("Genes:", nrow(merged_counts), "\n")
cat("File: expression_matrix_24samples.txt\n")
cat("========================================\n")