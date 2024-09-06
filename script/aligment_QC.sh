# Directories
bam_dir="../results/aligned"
results_dir="../results"
log_dir="../results/logs"

# Make directory
#mkdir -p $results_dir/alignment_QC
mkdir -p $results_dir/alignment_MultiQC_report

# Loop through the aligned files and run qualimap2 on each bam file.
#for bamfile in $bam_dir/*.bam; do
#    echo "running Qualimap QC for $bamfile"
#   qualimap bamqc -bam "$bamfile" -outdir "$results_dir/alignment_QC/$(basename $bamfile .bam)_QC" > $log_dir/alignment_QC.log 2>&1
#done

#echo "Qualimap QC completed"

echo "Running MultiQC accross Qualimap files"
multiqc $results_dir -o $results_dir/alignment_MultiQC_report > $log_dir/alignment_multiQC_report.log 2>&1
echo "Process completed"