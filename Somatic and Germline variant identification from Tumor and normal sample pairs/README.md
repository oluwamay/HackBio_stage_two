# Identification of Somatic and Germline variant from Tumor and Normal sample pairs
## Reproduced by Oluwayowa Samuel
###### INTRODUCTION
Mutation in molecular biology is the alteration in the nucleic acid sequence of a genome
creating an abnormal protein or preventing the formation of one. It is responsible for 
normal and abnormal (most times harmful) biological processes such as; Evolution, Cancer
and development of the immune system. Mutation of tumor suppressor gene is implicated in 
most cases of cancer in humans. The Tumor Suppressor genes are responsible for the 
monitoring of the rate of cell division, repairing of mismatched DNA and the control of 
when a cell dies. So when a Tumor Suppressor gene mutates it leads to uncontrollable cell 
division. These mutations can be Somatic (acquired during a lifetime) or Germline (inherited 
from parents). Germline variants are identified by comparing a sample genome to a reference 
while Somatic variants are identified using normal and tumor tissue.
> In this project, the aim was to reproduce a workflow that identifies Germline and Somatic 
variants affected by loss of Heterozygosity using both a healthy and tumour tissue, from which 
we will report variant sites and genes affected that could likely be the cause of the disease.
Such insights can help us track tumorigenesis in patients and help in the diagnosis, prognosis,
developing and guiding therapeutic strategies.

Below is a graphical abstract summarizing the key steps taken to achieve this:
![Graphical Abstract](/20220615_201840_0000.png)

This tutorial was reproduced as a Linux pipeline.
# Go to Selection :
1. [Introduction](https://github.com/oluwamay/oluwamay/blob/main/Somatic%20and%20Germline%20variant%20identification%20from%20Tumor%20and%20normal%20sample%20pairs/README.md#introduction)
2. [Linux Pipeline](https://github.com/oluwamay/oluwamay/blob/main/Somatic%20and%20Germline%20variant%20identification%20from%20Tumor%20and%20normal%20sample%20pairs/README.md#linux-pipeline)
3. [Contributors]()

# Linux Pipeline.

## Dataset Description
The datasets used in this analysis (reads from human chromosomes 5, 12 and 17), were obtained 
from a cancer patient’s tumor and normal tissue samples. The normal tissue coudn't be the only
sample used because healthy tissue contains many variants and every individual inherits a
unique pattern of many variants from their parents. The samples (paired end) were two in number.

## Data Download 
The datasets were downloaded from Zenodo using the wget command.

#### Sample Dataset
```
echo -e "\n Downloading data... \n"
	
mkdir -p raw_data 
cd raw_data
	
wget https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r1_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r2_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r1_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r2_chr5_12_17.fastq.gz	
```
#### Reference Sequence
```
echo -e "\n Downloading reference sequence... \n"
	
wget https://zenodo.org/record/2582555/files/hg19.chr5_12_17.fa.gz

#unzip reference
unzip hg19.chr5_12_17.fa.gz
```
## Pre-Processing and Trimming
#### i. Quality Check
The reads quality were examined using fastqc and an aggregate report generated with multiqc.
###### Description
FastQC aims to provide a way to do quality control checks on sequence data. Within the > fastq file
is quality information that refers to the accuracy of each base call. This helps to determine any
irregularies or features that make affect your results such as adapter contamination.
###### Installation
```
conda install -c bioconda fastqc multiqc --yes
```
###### Command
```
echo -e "\n Data Preprocessing... \n"

mkdir -p Fastqc_Reports  #create directory for the fastqc output
```
```
#Qc on reads
for sample in `cat list.txt`
do
	fastqc raw_data/${sample}*.fastq.gz -o Fastqc_Reports
done

multiqc Fastqc_Reports -o Fastqc_Reports	
```
The multiqc report can be examined from [here](). From the report, the reads quality are great, a few 
adapters are however observed.
#### ii. Removing Low quality reads using Fastp
###### Description
> Fastp is a FASTQ data Pre-Processing tool, the algorithm has functions for quality control, trimming of 
Adapters, trimming by quality and read pruning. It also supports multi threading, it is believed to be 
faster than other FASTQ Pre-Processing tools.
###### Installation
```
sudo apt-get install -y fastp
```
###### Command
```



```
The post trimming multiqc report can be found [here]() It is evident from the report that the quality of the
reads improved having per base quality scores above 35 and no adapters observed. After trimming an average 
of 0.73% normal reads and 1.24% tumor reads were lost.
**Note:  To view the multiqc html reports download the files and view them from your browser.**
## Mapped Reads Processing.
