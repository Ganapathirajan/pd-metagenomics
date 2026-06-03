#!/bin/bash
# Metagenomics Pipeline Script — PRJNA834801
# Usage: bash run_pipeline.sh <SRR_ID>

SRR=$1 
set -e
BASE=~/metagenomics_project
DB=$BASE/kraken2_db

echo ">>> Starting pipeline for $SRR"

# Step 1: Download
echo ">>> [1/6] Downloading $SRR..."
prefetch $SRR -O $BASE/raw_data/
fasterq-dump $BASE/raw_data/$SRR \
    -O $BASE/raw_data/ --threads 4 --progress

# Step 2: FastQC
echo ">>> [2/6] Running FastQC..."
fastqc $BASE/raw_data/${SRR}_1.fastq \
       $BASE/raw_data/${SRR}_2.fastq \
    -o $BASE/fastqc_reports/ --threads 4

# Step 3: Trimmomatic
echo ">>> [3/6] Trimming adapters..."
trimmomatic PE \
    $BASE/raw_data/${SRR}_1.fastq \
    $BASE/raw_data/${SRR}_2.fastq \
    $BASE/trimmed/${SRR}_1_trimmed.fastq \
    $BASE/trimmed/${SRR}_1_unpaired.fastq \
    $BASE/trimmed/${SRR}_2_trimmed.fastq \
    $BASE/trimmed/${SRR}_2_unpaired.fastq \
    ILLUMINACLIP:$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50 -threads 4

# Step 4: Kraken2
echo ">>> [4/6] Running Kraken2..."
kraken2 --db $DB \
    --report $BASE/kraken2_output/${SRR}_report.txt \
    --output $BASE/kraken2_output/${SRR}_output.txt \
    --paired --threads 4 \
    $BASE/trimmed/${SRR}_1_trimmed.fastq \
    $BASE/trimmed/${SRR}_2_trimmed.fastq

# Cleanup — delete raw FASTQ to save space
echo ">>> Cleaning raw files..."
rm -f $BASE/raw_data/${SRR}_1.fastq
rm -f $BASE/raw_data/${SRR}_2.fastq
rm -f $BASE/raw_data/${SRR}/${SRR}.sra

# Step 5: Bracken
echo ">>> [5/6] Running Bracken..."
bracken -d $DB \
    -i $BASE/kraken2_output/${SRR}_report.txt \
    -o $BASE/bracken_output/${SRR}_bracken.txt \
    -w $BASE/bracken_output/${SRR}_bracken_report.txt \
    -r 150 -l S

# Cleanup — free disk space
echo ">>> [cleanup] Removing raw and trimmed files..."
rm -f $BASE/raw_data/${SRR}_1.fastq
rm -f $BASE/raw_data/${SRR}_2.fastq
rm -f $BASE/raw_data/${SRR}/${SRR}.sra
rm -f $BASE/trimmed/${SRR}_1_trimmed.fastq
rm -f $BASE/trimmed/${SRR}_2_trimmed.fastq
rm -f $BASE/trimmed/${SRR}_1_unpaired.fastq
rm -f $BASE/trimmed/${SRR}_2_unpaired.fastq

# Step 6: Krona
echo ">>> [6/6] Generating Krona chart..."
ktImportTaxonomy \
    -t 5 -m 7 \
    $BASE/bracken_output/${SRR}_bracken.txt \
    -o $BASE/krona_output/${SRR}_krona.html

echo ">>> Pipeline complete for $SRR ✅"
echo ">>> Krona chart: $BASE/krona_output/${SRR}_krona.html"
