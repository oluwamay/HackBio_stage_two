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
'''
echo -e "\n Downloading data... \n"
	
mkdir -p raw_data 
cd raw_data
	
wget https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r1_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r2_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r1_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r2_chr5_12_17.fastq.gz	
'''
#### Reference Sequence
'''
echo -e "\n Downloading reference sequence... \n"
	
wget https://zenodo.org/record/2582555/files/hg19.chr5_12_17.fa.gz

#unzip reference
unzip hg19.chr5_12_17.fa.gz
'''
## Pre-Processing and Trimming


