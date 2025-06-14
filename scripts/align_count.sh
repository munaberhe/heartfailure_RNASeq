#!/bin/bash

# =============================================
# Heart Failure RNA-Seq: STAR alignment + counting
# =============================================

# -------------------------------------------
# 1) Make STAR genome index (only once)
# Make sure you have ref/hg38.fa and ref/hg38.gtf in your ref/ folder
mkdir -p ref/star_index

STAR --runThreadN 4 --runMode genomeGenerate \
     --genomeDir ref/star_index \
     --genomeFastaFiles ref/hg38.fa \
     --sjdbGTFfile ref/hg38.gtf

# -------------------------------------------
# 2) Align trimmed reads
SAMPLES=(
    SRR2131556
    SRR2131557
)

for SAMPLE in "${SAMPLES[@]}"
do
    STAR --runThreadN 4 \
         --genomeDir ref/star_index \
         --readFilesIn data/trimmed/${SAMPLE}_R1.trimmed.fastq.gz data/trimmed/${SAMPLE}_R2.trimmed.fastq.gz \
         --readFilesCommand zcat \
         --outFileNamePrefix results/star/${SAMPLE}_ \
         --outSAMtype BAM SortedByCoordinate
done

# -------------------------------------------
# 3) Count genes with featureCounts
featureCounts -T 4 -a ref/hg38.gtf \
              -o results/featurecounts/counts.txt \
              results/star/*Aligned.sortedByCoord.out.bam

echo "✅ DONE: Alignment + counting finished!"
