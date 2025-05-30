#!/bin/bash

# === Set sample name here ===
SAMPLE_NAME=CDI-2-1-M_L1
# ============================

# Define paths

STAR_DIR=/home/STAR/${SAMPLE_NAME}_Shape_set
MUTATION_DIR=/home/mutation/${SAMPLE_NAME}_Shape_set
FASTP_DIR=/home/fastp/${SAMPLE_NAME}
GENOME_DIR=/home/ref/UCSC_38
REFERENCE_GENOME=/home/ref/UCSC_38/hg38.fa
RAW_DIR=/home/rawdata/${SAMPLE_NAME}

# Create directories
mkdir -p "${STAR_DIR}"	
mkdir -p "${MUTATION_DIR}"

#fastq remove adapters
fastp -i "${RAW_DIR}/${SAMPLE_NAME}_1.fq.gz" -I "${RAW_DIR}/${SAMPLE_NAME}_2.fq.gz" \
      -o "${FASTP_DIR}/${SAMPLE_NAME}_1.fq.gz" -O "${FASTP_DIR}/${SAMPLE_NAME}_2.fq.gz" \
      -h ""${FASTP_DIR}/${SAMPLE_NAME}_1.html""

# STAR alignment
STAR --readFilesCommand zcat \
     --runThreadN 16 \
     --genomeDir "${GENOME_DIR}" \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMattrRGline ID:"${SAMPLE_NAME}_Shape_set" PL:Illumina LB:lib1 SM:"${SAMPLE_NAME}_Shape_set" PU:unit1 \
     --readFilesIn "${FASTP_DIR}/${SAMPLE_NAME}_1.fq.gz" "${FASTP_DIR}/${SAMPLE_NAME}_2.fq.gz" \
     --outFileNamePrefix "${STAR_DIR}/${SAMPLE_NAME}_Shape_set" \
     --scoreGap -1000000 \
     --scoreDelBase -1 \
     --scoreInsBase -1 \
     --outFilterMismatchNmax 999 \
     --outFilterMismatchNoverLmax 999 \
     --outMultimapperOrder Random \
     --outSAMmultNmax 1 

# Build BAM index
samtools index "${STAR_DIR}/${SAMPLE_NAME}_Shape_setAligned.sortedByCoord.out.bam"

# Generate pileup file
samtools mpileup --reference "${REFERENCE_GENOME}" \
     "${STAR_DIR}/${SAMPLE_NAME}_Shape_setAligned.sortedByCoord.out.bam" \
     -d 100000 \
    > "${MUTATION_DIR}/${SAMPLE_NAME}_Shape_set.pileup"

# Detect variants using VarScan
varscan pileup2cns "${MUTATION_DIR}/${SAMPLE_NAME}_Shape_set.pileup" \
     --min-coverage 5 \
     --min-reads2 1 \
     --min-avg-qual 0 \
     --min-var-freq 15 \
    > "${MUTATION_DIR}/${SAMPLE_NAME}_varscan_pileup_2cns.vcf"

echo "Processing completed for sample: ${SAMPLE_NAME}"