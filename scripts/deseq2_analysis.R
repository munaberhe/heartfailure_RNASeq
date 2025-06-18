# -------------------------------------------
# 1. Load libraries
# -------------------------------------------
library(DESeq2)

# -------------------------------------------
# 2. Load PE counts (HeartFailure)
# -------------------------------------------
pe <- read.table("results/featurecounts/pe_counts_FAKE.txt", header=TRUE, skip=1, sep="\t")

# -------------------------------------------
# 3. Load SE counts (Control)
# -------------------------------------------
se <- read.table("results/featurecounts/se_counts_FAKE.txt", header=TRUE, skip=1, sep="\t")

# -------------------------------------------
# 4. Combine counts
# -------------------------------------------
counts <- data.frame(
  row.names = pe$Geneid,
  PE_rep1 = pe[["results.star.SRR2131556_Aligned.sortedByCoord.out.bam"]],
  PE_rep2 = pe[["results.star.SRR2131556_rep2.bam"]],
  SE_rep1 = se[["results.star.SRR2131557_Aligned.sortedByCoord.out.bam"]],
  SE_rep2 = se[["results.star.SRR2131557_rep2.bam"]]
)

# -------------------------------------------
# 5. Define condition factor
# -------------------------------------------
condition <- factor(c("HeartFailure", "HeartFailure", "Control", "Control"))

# -------------------------------------------
# 6. Create DESeq2 dataset
# -------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = data.frame(condition),
                              design = ~ condition)

# -------------------------------------------
# 7. Use gene-wise dispersions & manual test
# -------------------------------------------
dds <- estimateSizeFactors(dds)
dds <- estimateDispersionsGeneEst(dds)
dispersions(dds) <- mcols(dds)$dispGeneEst
dds <- nbinomWaldTest(dds)
res <- results(dds)

# -------------------------------------------
# 8. Save results
# -------------------------------------------
write.csv(as.data.frame(res), file = "results/deseq2_results.csv")

# -------------------------------------------
# 9. MA plot
# -------------------------------------------
pdf("results/MAplot.pdf")
plotMA(res, main="DESeq2 MA-Plot")
dev.off()

# -------------------------------------------
# 10. PCA plot (blind = TRUE for fake replicates)
# -------------------------------------------
vsd <- vst(dds, blind=TRUE)
pdf("results/PCAplot.pdf")
plotPCA(vsd, intgroup="condition")
dev.off()

# -------------------------------------------
# ✅ Done
# -------------------------------------------
message("✅ All done! Check results/deseq2_results.csv and plots in results/")

