#!/bin/bash

# =============================
# 1️⃣  Build STAR index (safe if it exists)
# =============================

echo -e "\n\tBuilding STAR index (if needed)...\n"
STAR --runThreadN 4 \
     --runMode genomeGenerate \
     --genomeDir ref/star_index \
     --genomeFastaFiles ref/hg38.fa \
     --sjdbGTFfile ref/hg38.gtf

# =============================
# 2️⃣  Align paired-end sample SRR2131556
# =============================

echo -e "\n\tAligning PAIRED-END SRR2131556...\n"
STAR --runThreadN 4 \
     --genomeDir ref/star_index \
     --readFilesIn data/trimmed/SRR2131556_R1.trimmed.fastq.gz data/trimmed/SRR2131556_R2.trimmed.fastq.gz \
     --readFilesCommand gunzip -c \
     --outFileNamePrefix results/star/SRR2131556_ \
     --outSAMtype BAM SortedByCoordinate

# =============================
# 3️⃣  Align single-end sample SRR2131557
# =============================

echo -e "\n\tAligning SINGLE-END SRR2131557...\n"
STAR --runThreadN 4 \
     --genomeDir ref/star_index \
     --readFilesIn data/trimmed/SRR2131557.trimmed.fastq.gz \
     --readFilesCommand gunzip -c \
     --outFileNamePrefix results/star/SRR2131557_ \
     --outSAMtype BAM SortedByCoordinate

# =============================
# 4️⃣  Run featureCounts on both BAMs
# =============================

echo -e "\n\tRunning featureCounts...\n"
featureCounts -T 4 \
  -a ref/hg38.gtf \
  -o results/featurecounts/counts.txt \
  results/star/SRR2131556_Aligned.sortedByCoord.out.bam \
  results/star/SRR2131557_Aligned.sortedByCoord.out.bam

echo -e "\n✅ DONE: Alignment + counting finished!"

