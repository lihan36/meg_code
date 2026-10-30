# ============================================================
# PCA Analysis for Ficedula mugimaki RNA-seq
# Four groups: Low-Male, Low-Female, High-Male, High-Female
# ============================================================


library(DESeq2)
library(ggplot2)
library(ggrepel)

# ============================================================
# 1. Load expression data
# ============================================================

counts <- read.table("../03_rnaseq/expression_matrix_24samples.txt",
                     header = TRUE,
                     row.names = 1,
                     check.names = FALSE)

cat("Expression matrix:", dim(counts)[1], "genes x", dim(counts)[2], "samples\n")

# ============================================================
# 2. Define sample metadata
# ============================================================

samples <- colnames(counts)

low_male    <- c("GGZ56", "GGZ57", "GGZ67", "GGZ68", "GG001", "GG003", "GG006", "GG007", "GG009")
low_female  <- c("GG012", "GG019", "GG004", "GG005", "GG008")
high_male   <- c("GG384", "GG386", "GG390", "GG392", "GG407")
high_female <- c("GG254", "GG395", "GG398", "GG403", "GG266")

sample_info <- data.frame(
  sample = samples,
  stringsAsFactors = FALSE
)

sample_info$group <- ifelse(sample_info$sample %in% low_male, "L-M",
                     ifelse(sample_info$sample %in% low_female, "L-F",
                     ifelse(sample_info$sample %in% high_male, "H-M",
                     ifelse(sample_info$sample %in% high_female, "H-F", NA))))

sample_info$altitude <- factor(
  ifelse(grepl("^L", sample_info$group), "1800m", "3000m"),
  levels = c("1800m", "3000m")
)

sample_info$sex <- factor(
  ifelse(grepl("M$", sample_info$group), "Male", "Female"),
  levels = c("Male", "Female")
)

sample_info$group_label <- factor(sample_info$group,
                                   levels = c("L-M", "L-F", "H-M", "H-F"))

rownames(sample_info) <- sample_info$sample

cat("\nSample groups:\n")
print(table(sample_info$group))

# ============================================================
# 3. Create DESeq2 object and vst transformation
# ============================================================

counts <- counts[, rownames(sample_info)]

dds <- DESeqDataSetFromMatrix(
  countData = round(counts),
  colData = sample_info,
  design = ~ altitude + sex + altitude:sex
)

# Filter low-expressed genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
cat("\nGenes after filtering:", sum(keep), "/", nrow(counts), "\n")

# VST normalization
cat("\nPerforming VST normalization...\n")
vsd <- vst(dds, blind = TRUE)

# ============================================================
# 4. PCA - using all genes
# ============================================================

cat("\nCalculating PCA...\n")
pca_data <- plotPCA(vsd, intgroup = c("group_label", "altitude", "sex"),
                     returnData = TRUE)

pca_result <- prcomp(t(assay(vsd)))
var_explained <- round((pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100, 2)

cat("\nVariance explained:\n")
cat("  PC1:", var_explained[1], "%\n")
cat("  PC2:", var_explained[2], "%\n")

# ============================================================
# 5. PCA plot with group colors
# ============================================================

group_colors <- c(
  "L-M" = "#2E86AB",
  "L-F" = "#A23B72",
  "H-M" = "#F18F01",
  "H-F" = "#C73E1D"
)

pdf("pca_4groups.pdf", width = 10, height = 8)

p1 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = group_label, shape = altitude)) +
  geom_point(size = 5, alpha = 0.9) +
  scale_color_manual(values = group_colors, name = "Group") +
  scale_shape_manual(values = c(16, 17), name = "Altitude") +
  labs(
    title = "PCA of RNA-seq Samples",
    subtitle = "Ficedula mugimaki - Four groups",
    x = paste0("PC1: ", var_explained[1], "% variance"),
    y = paste0("PC2: ", var_explained[2], "% variance")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11)
  )

print(p1)
dev.off()

cat("PCA plot saved: pca_4groups.pdf\n")

# ============================================================
# 6. PCA plot with sample labels
# ============================================================

pdf("pca_4groups_with_labels.pdf", width = 12, height = 8)

p2 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = group_label)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(aes(label = rownames(pca_data)),
                  size = 3.5, max.overlaps = 20, box.padding = 0.5) +
  scale_color_manual(values = group_colors, name = "Group") +
  labs(
    title = "PCA of RNA-seq Samples (with sample labels)",
    x = paste0("PC1: ", var_explained[1], "% variance"),
    y = paste0("PC2: ", var_explained[2], "% variance")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11)
  )

print(p2)
dev.off()

cat("PCA plot with labels saved: pca_4groups_with_labels.pdf\n")

# ============================================================
# 7. PCA with ellipses (95% confidence)
# ============================================================

pdf("pca_4groups_ellipse.pdf", width = 10, height = 8)

p3 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = group_label, fill = group_label)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(level = 0.95, alpha = 0.2, geom = "polygon") +
  scale_color_manual(values = group_colors, name = "Group") +
  scale_fill_manual(values = group_colors, name = "Group") +
  labs(
    title = "PCA with 95% Confidence Ellipses",
    subtitle = "Ficedula mugimaki - Four groups",
    x = paste0("PC1: ", var_explained[1], "% variance"),
    y = paste0("PC2: ", var_explained[2], "% variance")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11)
  )

print(p3)
dev.off()

cat("PCA plot with ellipses saved: pca_4groups_ellipse.pdf\n")

# ============================================================
# 8. PC1 vs PC3 plot
# ============================================================

pdf("pca_pc1_vs_pc3.pdf", width = 10, height = 8)

pca_result_full <- prcomp(t(assay(vsd)))
pca_scores <- as.data.frame(pca_result_full$x)
pca_scores$group_label <- sample_info$group_label
pca_scores$sample <- rownames(pca_scores)

p4 <- ggplot(pca_scores, aes(x = PC1, y = PC3, color = group_label)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = group_colors, name = "Group") +
  labs(
    title = "PCA: PC1 vs PC3",
    x = paste0("PC1: ", var_explained[1], "% variance"),
    y = paste0("PC3: ", var_explained[3], "% variance")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11)
  )

print(p4)
dev.off()

cat("PC1 vs PC3 plot saved: pca_pc1_vs_pc3.pdf\n")

# ============================================================
# 9. Save PCA data for external use
# ============================================================

pca_coords <- pca_data
pca_coords$sample <- rownames(pca_coords)
write.table(pca_coords, "pca_coordinates.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

pca_loadings <- as.data.frame(pca_result_full$rotation)
pca_loadings$gene <- rownames(pca_loadings)
pca_loadings <- pca_loadings[, c("gene", "PC1", "PC2", "PC3")]

top_pc1 <- pca_loadings[order(abs(pca_loadings$PC1), decreasing = TRUE), ]
top_pc1 <- top_pc1[1:50, ]
write.table(top_pc1, "pca_top50_PC1_genes.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

top_pc2 <- pca_loadings[order(abs(pca_loadings$PC2), decreasing = TRUE), ]
top_pc2 <- top_pc2[1:50, ]
write.table(top_pc2, "pca_top50_PC2_genes.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nPCA data saved:\n")
cat("  pca_coordinates.txt\n")
cat("  pca_top50_PC1_genes.txt\n")
cat("  pca_top50_PC2_genes.txt\n")

# ============================================================
# 10. Session information
# ============================================================

cat("\n", rep("=", 60), "\n")
cat("Session information\n")
cat(rep("=", 60), "\n")

sessionInfo()

cat("\nPCA analysis completed:", format(Sys.time()), "\n")

cat("\nOutput files generated:\n")
cat("  pca_4groups.pdf\n")
cat("  pca_4groups_with_labels.pdf\n")
cat("  pca_4groups_ellipse.pdf\n")
cat("  pca_pc1_vs_pc3.pdf\n")
cat("  pca_coordinates.txt\n")
cat("  pca_top50_PC1_genes.txt\n")
cat("  pca_top50_PC2_genes.txt\n")