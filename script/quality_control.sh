#!/bin/bash
# A script to perform quality control on RNA-Seq data using FastQC and FastQ Screen

# Variables for directories
DATA_DIR="../data/raw"
RESULTS_DIR="../results"
FASTQC_DIR="$RESULTS_DIR/fastqc"
MULTIQC_DIR="$RESULTS_DIR/multiqc"
LOG_DIR="$RESULTS_DIR/logs"

# Make directories
mkdir -p "$FASTQC_DIR" "$MULTIQC_DIR" "$LOG_DIR"

# Run FastQC
echo "Running FastQC..."
fastqc --outdir $FASTQC_DIR $DATA_DIR/*.fastq.gz 2>&1 | tee $LOG_DIR/fastqc.log

# Run MultiQC
echo "Running MultiQC report..."
multiqc -o $MULTIQC_DIR $FASTQC_DIR 2>&1 | tee $LOG_DIR/multiqc.log

echo "Quality control complete. Check $RESULTS_DIR for output files and $LOG_DIR for logs."
