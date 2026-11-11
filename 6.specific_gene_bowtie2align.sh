#!/bin/bash

# === Set sample name here ===
SAMPLE_NAME= SRT_1
# ============================

# Define paths
BOWTIE2_DIR=/home/bowtie2/${SAMPLE_NAME}
MUTATION_DIR=/home/mutation/${SAMPLE_NAME}
FASTP_DIR=/home/fastp/${SAMPLE_NAME}
REFERENCE_INDEX=/home/ref/specific_ref
REFERENCE_GENOME=/home/ref/specific_ref/specific_ref.fa
RAW_DIR=/home/rawdata/${SAMPLE_NAME}

# Create directories
mkdir -p "${BOWTIE2_DIR}"
mkdir -p "${MUTATION_DIR}"

#fastq remove adapters
fastp -i "${RAW_DIR}/${SAMPLE_NAME}_1.fq.gz" -I "${RAW_DIR}/${SAMPLE_NAME}_2.fq.gz" \
      -o "${FASTP_DIR}/${SAMPLE_NAME}_1.fq.gz" -O "${FASTP_DIR}/${SAMPLE_NAME}_2.fq.gz" \
      -h ""${FASTP_DIR}/${SAMPLE_NAME}_1.html""

# Merge paired-end reads into interleaved FASTQ format
seqtk mergepe \
    "${FASTP_DIR}/${SAMPLE_NAME}_2.fq.gz" \
    "${FASTP_DIR}/${SAMPLE_NAME}_1.fq.gz" \
    > "${FASTP_DIR}/${SAMPLE_NAME}_PE.fq"
gzip "${FASTP_DIR}/${SAMPLE_NAME}_PE.fq"

# Remove duplicate reads by sequence
seqkit rmdup -s \
    "${FASTP_DIR}/${SAMPLE_NAME}_PE.fq.gz" \
    -o "${FASTP_DIR}/${SAMPLE_NAME}_PE_trim.fq.gz"

# Remove UMI (first 5 and last 5 bases)
cutadapt -u 5 -u -5 \
    -o "${FASTP_DIR}/${SAMPLE_NAME}_PE_trim_rmUMI.fq.gz" \
    "${FASTP_DIR}/${SAMPLE_NAME}_PE_trim.fq.gz"

# Align using Bowtie2 with local sensitive settings
bowtie2 -p 8 \
    -x "${REFERENCE_INDEX}" \
    --interleaved "${FASTP_DIR}/${SAMPLE_NAME}_PE_trim_rmUMI.fq.gz" \
    -S "${BOWTIE2_DIR}/${SAMPLE_NAME}_shape_set.sam" \
    --local \
    --sensitive-local \
    --mp 3,1 \
    --rdg 5,1 \
    --rfg 5,1 \
    --ignore-quals \
    2> "${BOWTIE2_DIR}/${SAMPLE_NAME}_shape_set_result.log"

# Sort BAM file
samtools sort -@ 3 -l 9 -m 4G \
    -o "${BOWTIE2_DIR}/${SAMPLE_NAME}_shape_set_sorted.bam" \
    "${BOWTIE2_DIR}/${SAMPLE_NAME}_shape_set.sam"

# Clean up intermediate SAM file
rm -f "${BOWTIE2_DIR}/${SAMPLE_NAME}_shape_set.sam"

# Generate pileup file
samtools mpileup --reference "${REFERENCE_GENOME}" \
    "${BOWTIE2_DIR}/${SAMPLE_NAME}_shape_set_sorted.bam" \
    -d 100000 \
    > "${MUTATION_DIR}/${SAMPLE_NAME}.pileup"

# Detect variants using VarScan
varscan pileup2cns "${MUTATION_DIR}/${SAMPLE_NAME}.pileup" \
    --min-coverage 5 \
    --min-reads2 1 \
    --min-avg-qual 0 \
    --min-var-freq 15 \
    > "${MUTATION_DIR}/${SAMPLE_NAME}_varscan_pileup_2cns.vcf"

echo "Processing completed for sample: ${SAMPLE_NAME}"
