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

This repository contains a reproducible RNA-Seq analysis pipeline for a heart failure project. It demonstrates the key steps from raw FASTQ data to gene expression quantification and differential expression analysis using STAR, featureCounts, and DESeq2.

---

## Project Structure

```
HeartFailure_RNASeq/
├── data/                   # Raw and trimmed FASTQ files
├── ref/                    # Reference genome FASTA and GTF
├── results/                
│   ├── star/               # STAR alignments
│   ├── featurecounts/      # Count matrices
│   ├── MAplot.pdf          # DESeq2 MA plot
│   ├── PCAplot.pdf         # DESeq2 PCA plot
│   └── deseq2_results.csv  # DESeq2 results table
├── scripts/                
│   ├── download_qc.sh      # Download and trim FASTQ
│   ├── align_count.sh      # STAR indexing, alignment, counting
│   ├── align_only.sh       # STAR alignment only
│   └── deseq2_analysis.R   # R script for DESeq2 analysis
└── README.md               # Project documentation
```

---

## Pipeline Steps

1. **Download and Quality Control**  
   Run:
   ```bash
   bash scripts/download_qc.sh
   ```
   Downloads FASTQ files and performs trimming.

2. **Alignment and Counting**  
   Run:
   ```bash
   bash scripts/align_count.sh
   ```
   - Generates STAR genome index.  
   - Aligns reads and outputs sorted BAM files.  
   - Counts reads per gene with featureCounts.

   Or use:
   ```bash
   bash scripts/align_only.sh
   ```
   to run alignment only.

3. **Differential Expression Analysis**  
   Run:
   ```bash
   Rscript scripts/deseq2_analysis.R
   ```
   Generates `deseq2_results.csv`, `MAplot.pdf`, and `PCAplot.pdf`.

## Workflow Overview

1. **Download & QC**
   - Download FASTQ files from SRA
   - Perform quality control and trimming using `fastp`


2. **Alignment**
   - Align reads to the human genome (hg38) using `STAR`
   - Sort BAM files with `samtools`

- STAR
- Subread (featureCounts)
- R with DESeq2 and dependencies
- FASTQ data and reference genome files


3. **Quantification**
   - Count gene-level reads with `featureCounts`

4. **Differential Expression**
   - Analyze counts using `DESeq2` in R

## Notes

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

- `results/featurecounts/counts.txt`: Gene count matrix
- `results/deseq2_results.csv`: DESeq2 differential expression results
- `results/pca_plot.png`: PCA plot of samples
- `results/volcano_plot.png`: Volcano plot of significant genes
- `results/heatmap_top20.png`: Heatmap of top differentially expressed genes

---

## Reference Files

This pipeline requires the human reference genome (FASTA) and gene annotation (GTF). These files are not included in this repository due to size constraints. Please download them from [GENCODE](https://www.gencodegenes.org/human/) and place them in the `ref/` folder.

**Recommended files:**
- `GRCh38.primary_assembly.genome.fa`
- `gencode.v44.annotation.gtf`

---

## Citation

If you use or adapt this pipeline, please cite the original dataset [GSE71613](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71613) and acknowledge all tools and packages according to their respective licenses.

- This repository uses example data with simulated replicates for demonstration only.
- For real analysis, use proper biological replicates.


---

## License

This project is open for educational and research use.

---

## Citation

If you use this workflow, please cite:
- [STAR](https://github.com/alexdobin/STAR)
- [Subread / featureCounts](http://subread.sourceforge.net/)
- [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)

---

## Contact

For questions or suggestions, please use the Issues tab.

---

## How to Use

Clone this repository and adjust the scripts to match your data and reference files. All paths are relative; run commands from the project root.

---

**Author:** Muna Berhe  
**Date:** June 2025

