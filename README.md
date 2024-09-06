# Identification of key genes for cocaine addiction using integerated bioinformatics analysis

[Wang X et al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10352680/) proposed a study to reveal the mechanism of cocaine addiction and identification of key genes that play important roles in cocaine addiction.

In their study, they identified four key genes (JUN, FOS, EGR1 and IL6) that might serve as biomarkers of cocaine addiction. They used publicly available data from the Gene Expression Omnibus (GEO) with the experimental dataset [GSE54839](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54839) containing 10 samples of cocaine treated with three replicates each and 10 samples of contol with three replicates each. They further used other dataset as veryfication for their results as shown in the figure below.
![working dataset](data/working_dataset.png)

## Aim of this study

In this study, we aim to replicate results from Wang et al, we choosed to use a different approach. They developed their eperimental dataset to confirm with verification dataset, we developed their verification dataset and compared with their experimental dataset. our experimental dataset which is their verification dataset [GSE186981](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186981) is an RNA-Seq data in HMDP mouse strains of nucleus accumbens and prefrontal cortex brain regions. it contains 392 samples accross different strains between mail and female. Below is the dataset overview from [GEMMA](https://gemma.msl.ubc.ca/expressionExperiment/showExpressionExperiment.html?id=21038)
![gemma](data/gemma.png)

In this dataset, due to it's large size, we only consedered one strain A/J, below is the snipset overview of the selected strains. dataset contains 16 runs and 8 biosamples treated with saline and cocaine. Each sample has two replicates.














#Download files using sra run:
prefetch SRRXXXXXXX
#download assession files
prefetch --option-file /path/to/your/accession_list.txt

#convert to fastq.gz
fastq-dump --split-files SRRXXXXXXX.sra
#convert multi files
for sra_file in *.sra; do
  fastq-dump --split-files --gzip "$sra_file"
done

interleaved FASTQ format, where paired reads are stored consecutively within a single file.
