# ============================================================
# Differential expression analysis with DESeq2
# ============================================================


library(DESeq2)
library(ggplot2)

# 1. Load expression matrix and sample info
counts <- read.table("../03_rnaseq/expression_matrix_24samples.txt",
                     header = TRUE, row.names = 1, check.names = FALSE)
sample_info <- read.table("../03_rnaseq/sample_info_full.txt",
                          header = TRUE, stringsAsFactors = FALSE)

# 2. Filter target samples
target_groups <- c("L-F", "H-F", "L-M", "H-M")
keep_samples <- sample_info$sample[sample_info$group %in% target_groups]
counts_filtered <- counts[, keep_samples, drop = FALSE]
coldata <- sample_info[match(colnames(counts_filtered), sample_info$sample), ]
rownames(coldata) <- coldata$sample

# Set factor levels
coldata$group <- factor(coldata$group, levels = c("L-F", "H-F", "L-M", "H-M"))
cat("Group counts:\n")
print(table(coldata$group))

# 3. Build DESeq2 object
counts_filtered <- round(counts_filtered)
dds <- DESeqDataSetFromMatrix(countData = counts_filtered,
                              colData = coldata,
                              design = ~ group)

# Pre-filtering: keep genes with counts >= 10 in at least 3 samples
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
cat("Genes after filtering:", sum(keep), "/", nrow(counts_filtered), "\n")

# 4. Run DESeq2
dds <- DESeq(dds)

# 5. Extract results
res_female <- results(dds, contrast = c("group", "H-F", "L-F"))
res_male <- results(dds, contrast = c("group", "H-M", "L-M"))

res_female <- res_female[order(res_female$padj), ]
res_male <- res_male[order(res_male$padj), ]

# 6. Filter significant DEGs (padj < 0.05 and |log2FC| > 1)
padj_cutoff <- 0.05
lfc_cutoff <- 1

filter_deg <- function(res, padj_cut, lfc_cut) {
  res_clean <- res[!is.na(res$padj) & !is.na(res$log2FoldChange), ]
  deg <- rownames(res_clean[res_clean$padj < padj_cut & abs(res_clean$log2FoldChange) > lfc_cut, ])
  return(deg)
}

female_alt_genes <- filter_deg(res_female, padj_cutoff, lfc_cutoff)
male_alt_genes <- filter_deg(res_male, padj_cutoff, lfc_cutoff)

cat("Female altitude response genes:", length(female_alt_genes), "\n")
cat("Male altitude response genes:", length(male_alt_genes), "\n")

# 7. Check up/down regulation
female_up <- rownames(res_female[res_female$padj < padj_cutoff &
                                  res_female$log2FoldChange > lfc_cutoff &
                                  !is.na(res_female$padj), ])
female_down <- rownames(res_female[res_female$padj < padj_cutoff &
                                    res_female$log2FoldChange < -lfc_cutoff &
                                    !is.na(res_female$padj), ])
male_up <- rownames(res_male[res_male$padj < padj_cutoff &
                              res_male$log2FoldChange > lfc_cutoff &
                              !is.na(res_male$padj), ])
male_down <- rownames(res_male[res_male$padj < padj_cutoff &
                                res_male$log2FoldChange < -lfc_cutoff &
                                !is.na(res_male$padj), ])

cat("Female up:", length(female_up), " down:", length(female_down), "\n")
cat("Male up:", length(male_up), " down:", length(male_down), "\n")

# 8. Compare female and male response
common <- intersect(female_alt_genes, male_alt_genes)
female_only <- setdiff(female_alt_genes, male_alt_genes)
male_only <- setdiff(male_alt_genes, female_alt_genes)

cat("Common genes:", length(common), "\n")
cat("Female-specific genes:", length(female_only), "\n")
cat("Male-specific genes:", length(male_only), "\n")

# 9. Save results
output_dir <- "DE_results"
if (!dir.exists(output_dir)) dir.create(output_dir)

write.table(female_alt_genes, file.path(output_dir, "female_altitude_genes.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(male_alt_genes, file.path(output_dir, "male_altitude_genes.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(common, file.path(output_dir, "common_altitude_genes.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(female_only, file.path(output_dir, "female_specific_genes.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(male_only, file.path(output_dir, "male_specific_genes.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)

write.csv(as.data.frame(res_female), file.path(output_dir, "H-F_vs_L-F_full_results.csv"))
write.csv(as.data.frame(res_male), file.path(output_dir, "H-M_vs_L-M_full_results.csv"))

cat("All DE results saved to:", output_dir, "\n")

# 10. Volcano plots
plot_volcano <- function(res, title, filename) {
  df <- data.frame(log2FC = res$log2FoldChange, padj = -log10(res$padj))
  df <- df[!is.na(df$padj) & !is.na(df$log2FC) & is.finite(df$padj), ]

  p <- ggplot(df, aes(x = log2FC, y = padj)) +
    geom_point(alpha = 0.5, size = 0.8) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
    labs(title = title, x = "log2 Fold Change", y = "-log10(adjusted p-value)") +
    theme_minimal()

  ggsave(filename, p, width = 8, height = 6)
}

plot_volcano(res_female, "Female: H-F vs L-F", file.path(output_dir, "volcano_female.png"))
plot_volcano(res_male, "Male: H-M vs L-M", file.path(output_dir, "volcano_male.png"))
cat("Volcano plots saved.\n")

# 11. Venn diagram (base R + plotrix)
if (length(female_alt_genes) > 0 && length(male_alt_genes) > 0) {
  if (!requireNamespace("plotrix", quietly = TRUE)) {
    install.packages("plotrix", repos = "https://cloud.r-project.org/")
  }
  library(plotrix)

  png(file.path(output_dir, "venn_diagram.png"),
      width = 800, height = 800, res = 120)

  n_female <- length(female_alt_genes)
  n_male <- length(male_alt_genes)
  n_common <- length(common)
  n_female_only <- length(female_only)
  n_male_only <- length(male_only)

  plot(1, type = "n", xlim = c(0, 10), ylim = c(0, 10),
       xlab = "", ylab = "", axes = FALSE)

  draw.circle(3.5, 5, 3, col = rgb(1, 0, 0, 0.3), border = "red", lwd = 2)
  draw.circle(6.5, 5, 3, col = rgb(0, 0, 1, 0.3), border = "blue", lwd = 2)

  text(2, 9, paste("Female", n_female), cex = 1.5, col = "red")
  text(8, 9, paste("Male", n_male), cex = 1.5, col = "blue")
  text(5, 5, paste("Common\n", n_common), cex = 1.5, font = 2)
  text(2.5, 4.5, n_female_only, cex = 1.2)
  text(7.5, 4.5, n_male_only, cex = 1.2)

  dev.off()
  cat("Venn diagram saved.\n")
}

# 12. Statistical summary
cat("\n========================================\n")
cat("Analysis Summary\n")
cat("========================================\n")
cat(sprintf("%-30s %s\n", "Item", "Count"))
cat(sprintf("%-30s %d\n", "Female altitude response genes", length(female_alt_genes)))
cat(sprintf("%-30s %d\n", "Male altitude response genes", length(male_alt_genes)))
cat(sprintf("%-30s %d\n", "Common response genes", length(common)))
cat(sprintf("%-30s %d\n", "Female-specific genes", length(female_only)))
cat(sprintf("%-30s %d\n", "Male-specific genes", length(male_only)))
cat("========================================\n")