# Directories
DATA_DIR="../data/ref_genome"
RESULTS_DIR="../results"
LOG_DIR="../results/logs"
BAM_DIR="../results/aligned"

mkdir -p $RESULTS_DIR/featureCounts

# -O count reads that overlap more than one feature.
# -s for stranded reads so we set to 0 for unstranded
# -t CDS|genes|exons... count reads that overlap teses sequences, much greater portion overlapped to genes
# --ignoreDup ignore duplicates. Dups may come from highyly expressed genes, PCR amplification bias or technical artifects. 
# -T specifies the number of threads to use
# -g specifies attribute to use in the genome as gene identifier.

echo "Computing gene counts..."
featureCounts -O -t exon -T 4 -G $DATA_DIR/Mus_musculus.GRCm39.dna.primary_assembly.fa -a $DATA_DIR/Mus_musculus.GRCm39.112.gtf -o $RESULTS_DIR/featureCounts/tmp_counts.txt $BAM_DIR/*.bam > $LOG_DIR/featureCounts.log 2>&1



    
# The following commands are just making a nicer table header.
FILE_PATH=$RESULTS_DIR/aligned/
sed -i "s|$FILE_PATH||g" $RESULTS_DIR/featureCounts/tmp_counts.txt
sed "s/_trimmed_aligned.sam//g" $RESULTS_DIR/featureCounts/tmp_counts.txt > $RESULTS_DIR/featureCounts/featureCounts.txt 2>&1
rm $RESULTS_DIR/featureCounts/tmp_counts.txt
echo "End of computing genes counts."

# Sort and clean the output file
echo "cleaning counts..."
awk '{print $1"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23}' $RESULTS_DIR/featureCounts/featureCounts.txt > $RESULTS_DIR/featureCounts/cleaned_fcounts.txt 2>&1
sed -i '1d' $RESULTS_DIR/featureCounts/cleaned_fcounts.txt
echo "completed"
