#!/bin/bash

SRR_IDS=(
    SRR2131556
    SRR2131557
)

mkdir -p data/raw_fastq data/trimmed qc/raw qc/trimmed

echo "Downloading FASTQ files..."
for SRR in "${SRR_IDS[@]}"
do
    echo "Fetching $SRR ..."
    fasterq-dump $SRR -O data/raw_fastq
done

echo "Running FastQC..."
fastqc data/raw_fastq/*.fastq -o qc/raw/

echo "Running fastp trimming..."
for SRR in "${SRR_IDS[@]}"
do
    fastp \
    -i data/raw_fastq/${SRR}_1.fastq \
    -I data/raw_fastq/${SRR}_2.fastq \
    -o data/trimmed/${SRR}_R1.trimmed.fastq.gz \
    -O data/trimmed/${SRR}_R2.trimmed.fastq.gz \
    --html qc/trimmed/fastp_${SRR}.html \
    --thread 4
done

echo "✅ Done: Download + QC finished!"
