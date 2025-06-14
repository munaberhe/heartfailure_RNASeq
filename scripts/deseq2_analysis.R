# =============================================
# Heart Failure RNA-seq DESeq2 Analysis
# GSE71613 — Muna's Portfolio Project
# =============================================

library(DESeq2)
library(pheatmap)
library(ggplot2)

# 1) Load counts (skip first comment lines)
counts <- read.table("results/featurecounts/counts.txt", header=TRUE, row.names=1, comment.char="#")
# Remove extra columns to keep only counts:
counts <- counts[, 6:ncol(counts)]

# 2) Inspect:
print(dim(counts))
print(head(counts))

# 3) Build metadata
# Adjust this to match your samples:
metadata <- data.frame(
  row.names = colnames(counts),
  condition = c("Control", "Failing")
)

print(metadata)

# 4) Build DESeq2 object
dds <- DESeqDataSetFromMatrix(countData=counts, colData=metadata, design=~condition)

# Filter low counts
dds <- dds[ rowSums(counts(dds)) > 10, ]

# 5) Run DESeq2
dds <- DESeq(dds)
res <- results(dds)

# 6) Export results
write.csv(as.data.frame(res), file="results/deseq2_results.csv")

# 7) PCA
vsd <- vst(dds)
png("results/pca_plot.png", width=600, height=600)
plotPCA(vsd, intgroup="condition")
dev.off()

# 8) Volcano
res$padj[is.na(res$padj)] <- 1
png("results/volcano_plot.png", width=800, height=600)
plot(res$log2FoldChange, -log10(res$padj),
     pch=20,
     main="Volcano Plot",
     xlab="Log2 Fold Change",
     ylab="-Log10 Adjusted P-Value",
     col=ifelse(res$padj<0.05 & abs(res$log2FoldChange)>1, "red", "grey"))
dev.off()

# 9) Heatmap top 20 genes
top20 <- head(order(res$padj), 20)
mat <- assay(vsd)[top20, ]
png("results/heatmap_top20.png", width=800, height=800)
pheatmap(mat, cluster_rows=TRUE, cluster_cols=TRUE, annotation_col=metadata)
dev.off()

cat("✅ DESeq2 analysis complete! Results in 'results/' folder.\n")
