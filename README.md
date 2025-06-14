# HeartFailure_RNASeq

This repository demonstrates a simple, reproducible RNA-Seq analysis pipeline comparing heart failure vs. control samples.

## Project Structure


```
HeartFailure_RNASeq/
├── data/              # Raw and trimmed FASTQ files
├── ref/               # Reference genome and annotation (hg38)
├── scripts/           # Shell and R scripts
├── results/           # Alignment, counts, and output files
├── figures/           # Final plots (MA plot, PCA plot)
└── README.md          # Project description and instructions
```

```HeartFailure_RNASeq/
├── data/ # Raw and trimmed FASTQ files (local only, ignored)
├── qc/ # Quality control reports
├── ref/ # Reference genome and GTF (local only, ignored)
├── results/ # STAR alignments, featureCounts output, DESeq2 results, plots
├── scripts/ # Shell and R scripts for each step
├── HeartFailure_RNASeq.Rproj # RStudio project file
└── README.md # Project description and instructions```


## Workflow Overview

1. **Download & QC**
   - Download FASTQ files from SRA
   - Perform quality control and trimming using `fastp`

2. **Alignment**
   - Align reads to the human genome (hg38) using `STAR`
   - Sort BAM files with `samtools`

3. **Quantification**
   - Count gene-level reads with `featureCounts`

4. **Differential Expression**
   - Analyze counts using `DESeq2` in R

5. **Visualization**
   - Generate MA Plot and PCA Plot

## Example Plots

Below are the core outputs from the DESeq2 analysis:

**PCA Plot**

![PCA Plot](figures/PCAplot.pdf)

**MA Plot**

![MA Plot](figures/MAplot.pdf)

> Note: On GitHub, PDF images may not render inline. It is recommended to export PNG versions for web previews.

## How to Run

**Requirements:**

- Unix-like shell (Mac/Linux)
- `conda` environment with `STAR`, `samtools`, `R`, and `featureCounts`
- ~20 GB free disk space

**Typical commands:**

```bash
# 1) Activate environment
conda activate bio_env

# 2) Download and QC
bash scripts/download_qc.sh

# 3) Align reads
bash scripts/align_only.sh

# 4) Count reads
bash scripts/align_count.sh

# 5) Run DESeq2 analysis
Rscript scripts/deseq2_analysis.R
```

## Important

- This example uses duplicated BAM files to simulate biological replicates for demonstration.  
  For proper statistical power, always use at least three real replicates per condition.
- Large FASTQ and BAM files are not pushed to GitHub — only scripts, configs, and summary figures are versioned.

## License

This project is licensed under the MIT License.

---

**Author:** Muna Berhe  
**Date:** June 2025

