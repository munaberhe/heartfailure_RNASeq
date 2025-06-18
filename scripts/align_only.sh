#!/usr/bin/env bash

# ===========================
# align_only.sh
# ===========================

set -e

echo "🔹 Aligning PAIRED-END SRR2131556..."
STAR --runThreadN 4 \
 --genomeDir ref/star_index \
 --readFilesIn data/trimmed/SRR2131556_R1.trimmed.fastq.gz data/trimmed/SRR2131556_R2.trimmed.fastq.gz \
 --readFilesCommand gunzip -c \
 --outFileNamePrefix results/star/SRR2131556_ \
 --outSAMtype BAM SortedByCoordinate

echo "🔹 Aligning SINGLE-END SRR2131557..."
STAR --runThreadN 4 \
 --genomeDir ref/star_index \
 --readFilesIn data/trimmed/SRR2131557.trimmed.fastq.gz \
 --readFilesCommand gunzip -c \
 --outFileNamePrefix results/star/SRR2131557_ \
 --outSAMtype BAM SortedByCoordinate

echo "✅ DONE: Alignment finished!"

