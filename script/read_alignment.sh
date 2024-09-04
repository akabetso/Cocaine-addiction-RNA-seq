# After prepairing the reads, we proceed with aligning the 12 samples to the reference geonme. I love to use bowtie2 for paired-end reads alignment, here, for this single-end and unstranded reads, galaxy suggest the use of HISAT2 descendant of TopHat which is one of the first widely used aligners. there are many others like popular aligners like star.

# For the reference genome, the mm39/GRCm38 mouse reference genome released in june 2020 was downloaded from ncbi

# Directories
GENOME_DIR="../data/ref_genome"
RESULTS_DIR="../results"
LOG_DIR="../results/logs"
DATA_DIR="../data/raw"  # Directory containing your trimmed FASTQ files

mkdir -p $RESULTS_DIR/aligned
mkdir -p $RESULTS_DIR/un_aligned

# Build HITsAT2 Index, this will load the genome in the memory for fast mapping.
#####unzip $GENOME_DIR/mm39_GRCm38_ncbi_dataset.zip -d $GENOME_DIR #unzip the downloaded file
#####hisat2-build $GENOME_DIR/ncbi_dataset/data/GCF_000001635.27/GCF_000001635.27_GRCm39_genomic.fna $GENOME_DIR/mm39
HISAT2_INDEX=$GENOME_DIR/hisat2_index/mm39

# Mapping single-end unstranded reads
# -command finds the indexed path
# -u command spicifies the input files
# -s command speicifies the output files (sam files)
# --un outputs unmapped reads
# --rna-strandness R stands for unstranded reads.
# and finally the number of threads -p4 to speedup the alignment.
for SAMPLE in $(ls $DATA_DIR | grep ".fastq.gz" | sed 's/.fastq.gz//'); do
    echo "mapping $SAMPLE ..."
    hisat2 -q --rna-strandness R -p 4 -x $HISAT2_INDEX --un $RESULTS_DIR/un_aligned/${SAMPLE}_unmapped.fastq.gz -U $DATA_DIR/${SAMPLE}.fastq.gz | samtools sort -o $RESULTS_DIR/aligned/${SAMPLE}_aligned.bam 2>&1 | tee -a $LOG_DIR/hisat2.log            
done
echo "Alignment completed, check $RESULTS_DIR for SAM files and $LOG_DIR for log files"
