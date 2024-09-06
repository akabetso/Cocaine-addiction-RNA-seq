# Identification of key genes for cocaine addiction using integerated bioinformatics analysis

[Wang X et al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10352680/) proposed a study to reveal the mechanism of cocaine addiction and identification of key genes that play important roles in cocaine addiction.

In their study, they identified four key genes (JUN, FOS, EGR1 and IL6) that might serve as biomarkers of cocaine addiction. They used publicly available data from the Gene Expression Omnibus (GEO) with the experimental dataset [GSE54839](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54839) containing 10 samples of cocaine treated with three replicates each and 10 samples of contol with three replicates each. They further used other dataset as veryfication for their results as shown in the figure below.

![working dataset](data/working_dataset.png)

## Aim of this study

In this study, we aim to replicate results from Wang et al, we choosed to use a different approach. They developed their eperimental dataset to confirm with verification dataset, we developed their verification dataset and compared with their experimental dataset. our experimental dataset which is their verification dataset [GSE186981](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186981) is an RNA-Seq data in HMDP mouse strains of nucleus accumbens and prefrontal cortex brain regions. it contains 392 samples accross different strains between mail and female. Below is the dataset overview from [GEMMA](https://gemma.msl.ubc.ca/expressionExperiment/showExpressionExperiment.html?id=21038)
![gemma](data/gemma.png)

In this dataset, due to it's large size, we only consedered one strain A/J, below is the snipset overview of the selected strains. dataset contains 16 runs and 8 biosamples treated with saline and cocaine. Each sample has two replicates.

![dataset_overview](data/dataset_overview.png)

## Workflow

1. Pre-analysis
	Library preparation
	Quality control
2. Core-analysis
	Transcriptome profiling (alignment, quality control, Quantification)
	Differential expression
	Interpretation
	
### Preanlysis

According [Conesa A et al](https://pubmed.ncbi.nlm.nih.gov/26813401/), the first important step in RNA_seq study is to first choosing an appropriate library type, sequencing depth and number of replicates for the biological system under study. To provide you with a brief rundown of library preparation from the perspective of [Wang J et](https://pubmed.ncbi.nlm.nih.gov/30297273/) in a similar case with cancer cells. First the RNA is extracted from the cells and the subsets or RNA molecules are isolated. mRNA are mostly extracted using poly-A selection to prevent degration of mRNA. small RNAs lacking poly-A tail can be selected using size selection. the process is briefed in the diagram below.

![Library_prep](data/library_prep.png)

They further proposed a post library preparation workflow. Millions of short reads sequenced from PCR amplified of cDNA adapter ligaged fragments (reasons why it's important to check adapters and trim). Reads can be derived from single-end or paired end sequencing which can then be aligned to reference genome followed by down streams analysis. see the diagram below for clarity.

![steps](data/steps.png)

### Library preparation

As mentioned earlier, the said dataset was downloaded from GEO using the SRA_ACC_list and the preftech linux command line tool. later the data was unziped and dumped to fasq using the fastq-dump with the --split-files option for processing paired reads. 

### Quality control

As described by [Guo Y et al](https://academic.oup.com/bib/article/15/6/879/180439), Quality control is essential at three different stages: Raw data quality conrol, alignment quality conrol and variant calling quality control.

FastQC was used to study the quality control of the raw data, the quality showed less than o.1% of adapter contamination and phred score > 30. we saw minimal reason to trim the data and proceeded directly with alignment using HISAT2. After alignment, Qualimap was used to look into the alignment coverage and the overal quality of the data at this stage. multiQC was used to generate i nice report of the entire quality control. a snipset of the general statistics is shown below.

![quality_control](data/quality_control.png)

Here, we have a nice view of the general quality statistics. the alignment coverage, we observed that only a small portion of the genome is covered by our sequencing reads. but this is okay because we are not performing a whole genome sequence. 

We also see that, the QC content for Qualimap et fastQC are the same confirming minimal errors during alignment with GC alignment content > 95%.

The futureCounts of allignment that were qauntified and assigned were between 65% to 70%. FastQC showed a duplication reads rate of 40% averaged which impacts the unassigned reads.



The mm10 reference genome had been previously downloaded with and indexed with hitsat2.










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
