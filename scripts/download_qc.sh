# =============================================
# Download + QC script for Heart Failure RNA-seq
# Handles mixed paired-end and single-end runs
# =============================================

# Make folders if they don't exist
mkdir -p data/raw_fastq
mkdir -p data/trimmed

# === Define your samples ===
PAIRED_SAMPLES=(
    SRR2131556
)
SINGLE_SAMPLES=(
    SRR2131557
)

# === Download & trim PAIRED-END ===
for SAMPLE in "${PAIRED_SAMPLES[@]}"
do
    echo "Downloading PAIRED-END $SAMPLE..."
    fasterq-dump --split-files $SAMPLE -O data/raw_fastq/

    echo "Running fastp for PAIRED-END $SAMPLE..."
    fastp \
      -i data/raw_fastq/${SAMPLE}_1.fastq \
      -I data/raw_fastq/${SAMPLE}_2.fastq \
      -o data/trimmed/${SAMPLE}_R1.trimmed.fastq.gz \
      -O data/trimmed/${SAMPLE}_R2.trimmed.fastq.gz \
      --detect_adapter_for_pe \
      --thread 4
done

# === Download & trim SINGLE-END ===
for SAMPLE in "${SINGLE_SAMPLES[@]}"
do
    echo "Downloading SINGLE-END $SAMPLE..."
    fasterq-dump $SAMPLE -O data/raw_fastq/

    echo "Running fastp for SINGLE-END $SAMPLE..."
    fastp \
      -i data/raw_fastq/${SAMPLE}.fastq \
      -o data/trimmed/${SAMPLE}.trimmed.fastq.gz \
      --thread 4
done

echo "✅ DONE: Download + QC finished!"
EOF
