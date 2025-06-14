# 🫀 Heart Failure RNA-seq Analysis — GSE71613

## 📌 Overview
This project analyzes RNA-seq data comparing **failing** vs **non-failing human heart tissue**, using GEO dataset **GSE71613**.
It demonstrates a full pipeline from raw FASTQ download to DESeq2 differential expression analysis.

---

## 📂 Project Structure

HeartFailure_RNASeq/
├── data/ # Raw & trimmed FASTQ files
├── qc/ # Quality control reports
├── ref/ # Reference genome FASTA & GTF
├── results/ # Alignments, counts, plots
├── scripts/ # Shell & R scripts
├── HeartFailure_RNASeq.Rproj
└── README.md
---

## ✅ How to run

```bash
# 1) Download & QC
./scripts/download_qc.sh

# 2) Align & count
./scripts/align_count.sh

# 3) DESeq2 analysis
Rscript scripts/deseq2_analysis.R


---

### 📌 **Step 3 — Save and exit nano**

Inside nano:
- Press `CTRL + O` → then `Enter` to save
- Then `CTRL + X` to exit nano

---

### 📌 **Step 4 — Add, commit, push**

```bash
git add README.md
git commit -m "Add README.md"
git push

