# Identification of Somatic and Germline variant from Tumor and Normal sample pairs
## Reproduced by Oluwamayowa Samuel
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
gunzip hg19.chr5_12_17.fa.gz
```
## Pre-Processing and Trimming
### i. Quality Check
The reads quality were examined using fastqc and an aggregate report generated with multiqc.
##### Description
FastQC aims to provide a way to do quality control checks on sequence data. Within the `fastq` file
is quality information that refers to the accuracy of each base call. This helps to determine any
irregularies or features that make affect your results such as adapter contamination.
##### Installation
```
conda install -c bioconda fastqc multiqc --yes
```
##### Command
```
echo -e "\n Data Preprocessing... \n"

mkdir -p Fastqc_Reports  #create directory for the fastqc output
```
```
#Qc on reads
	fastqc *.fastq.gz -o Fastqc_Reports
        multiqc Fastqc_Reports -o Fastqc_Reports	
```
The multiqc report can be examined from [here](). From the report, the reads quality are great, a few 
adapters are however observed.
### ii. Removing Low quality reads using Trimmomatic
##### Description
`Trimmomatic` is a wrapper script that automate quality and adapter trimming. After analyzing data quality, the next step is to remove sequences that do not meet quality standards.
##### Installation
```
conda install -c bioconda trimmomatic --yes
```
##### Command
```
mkdir -p trimmed_reads

for sample in `cat list.txt`
do
       trimmomatic PE -threads 8 raw_data/${sample}_r1_chr5_12_17.fastq.gz raw_data/${sample}_r2_chr5_12_17.fastq.gz \
               trimmed_reads/${sample}_r1_paired.fq.gz trimmed_reads/${sample}_r1_unpaired.fq.gz \
               trimmed_reads/${sample}_r2_paired.fq.gz trimmed_reads/${sample}_r2_unpaired.fq.gz \
               ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:8:keepBothReads \
               LEADING:3 TRAILING:10 MINLEN:25
done
cd trimmed_reads

fastqc *.fq.gz -o Fastqc_results/

multiqc  Fastqc_results -o Fastgc_results
```
The post trimming multiqc report can be found [here]() It is evident from the report that the quality of the
reads improved having per base quality scores above 35 and no adapters observed. After trimming an average 
of 0.73% normal reads and 1.24% tumor reads were lost.
**Note:  To view the multiqc html reports download the files and view them from your browser.**
## Mapped Reads Processing.
### Description
Mapping of sample sequences against the reference genome is conducted with an aim of determining the most likey 
source of the observed sequencing reads. `BWA-MEM` was used for alignment. The results of mapping is a sequence 
alignment map (SAM) format. The file has a single unified format for storing read alignments to a reference genome.
#### Installation
conda install -y -c bioconda bwa
conda install -c bioconda samtools
conda install -c bioconda bamtools
#### Command
##### Read Mapping
In order to align the data, we need a reference to align against. First, a directory is created for the reference and 
then copied. The reference is indexed to be able to align the data.This is done using the command;
```
#Index reference file	
bwa index hg19.chr5_12_17.fa 
```
This produces 5 files in the reference directory that BWA uses during the alignment phase. The 5 
files have different extensions named amb,ann,bwt pac and sa. Alignment can be done using the program; 
`bwa mem`

Note that bwa is given a location , which is the path to the reference. Now, the two paired-end files 
are aligned and the alignment output (in SAM format) directed to a file. 24 threads (processors) were
used to speed up this process and a read group (i.e sample ID) information was added to the alignment:
```
mkdir Mapping
   
#Perform alignment
bwa mem -R '@RG\tID:231335\tSM:Normal' hg19.chr5_12_17.fa trimmed_reads/SLGFSK-N_231335_r1_paired.fq.gz \
      trimmed_reads/SLGFSK-N_231335_r2_paired.fq.gz > Mapping/SLGFSK-N_231335.sam

bwa mem -R '@RG\tID:231336\tSM:Tumor' hg19.chr5_12_17.fa trimmed_reads/SLGFSK-T_231336_r1_paired.fq.gz \
       trimmed_reads/SLGFSK-T_231336_r2_paired.fq.gz > Mapping/SLGFSK-T_231336.sam	


```
##### Converting SAM files to BAM files, sorting and indexing.
A Binary Alignment Map (BAM) format is an equivalent to sam but its developed for fast processing and 
indexing. It stores every read base, base quality and uses a single conventional technique for all types
 of data. The produced BAM files were sorted by read name and indexing was done for faster or rapid retrieval. 
At the end of the every BAM file, a special end of file (EOF) marker is usually written, the samtools index 
command also checks for this and produces an error message if its not found.
```
for sample in `cat list.txt`
do
        Convert SAM to BAM and sort it 
        samtools view -@ 20 -S -b Mapping/${sample}.sam | samtools sort -@ 32 > Mapping/${sample}.sorted.bam
        
        Index BAM file
        samtools index Mapping/${sample}.sorted.bam
done
```
##### Mapped Reads filtered
```
for sample in `cat list.txt`
do
	#Filter BAM files
        samtools view -q 1 -f 0x2 -F 0x8 -b Mapping/${sample}.sorted.bam > Mapping/${sample}.filtered1.bam
done
```
View the output of the result:
```
samtools flagstat <bam file>
```
##### Duplicate Removal
During library construction sometimes there's introduction of PCR (Polymerase Chain Reaction) duplicates, 
these duplicates usually can result in false SNPs (Single Nucleotide Polymorphisms), whereby the can manifest 
themselves as high read depth support. A low number of duplicates (<5%) in good libraries is considered standard.
```
#use the command markdup
for sample in `cat list.txt`
do
	samtools collate -o Mapping/${sample}.namecollate.bam Mapping/${sample}.filtered1.bam
        samtools fixmate -m Mapping/${sample}.namecollate.bam Mapping/${sample}.fixmate.bam
        samtools sort -@ 32 -o Mapping/${sample}.positionsort.bam Mapping/${sample}.fixmate.bam
        samtools markdup -@32 -r Mapping/${sample}.positionsort.bam Mapping/${sample}.clean.bam
done
	
#or rmdup
samtools rmdup SLGFSK35.sorted.bam  SLGFSK35.rdup and samtools rmdup SLGFSK36.sorted.bam  SLGFSK36.rdup
```
##### Left Align BAM
```
for sample in `cat list.txt`
do      
        cat Mapping/${sample}.clean.bam  | bamleftalign -f hg19.chr5_12_17.fa -m 5 -c > Mapping/${sample}.leftAlign.bam

#-c - compressed, -m - max-iterations
```
##### Recalibrate Read Mapping qualities
```
for sample in `cat list.txt`
do
        samtools calmd -@ 32 -b Mapping/${sample}.leftAlign.bam hg19.chr5_12_17.fa > Mapping/${sample}.recalibrate.bam
done
```
##### Refilter read mapping qualities
```
for sample in `cat list.txt`
do
        bamtools filter -in Mapping/${sample}.recalibrate.bam -mapQuality <=254 > Mapping/${sample}.refilter.bam
done
```
## Variant Calling and Classification
http://varscan.sourceforge.net/somatic-calling.html
### Description
To be able to identify variants from the mapped samples, the tool `VarScan somatic` was used. The command expects 
both a normal and tumor sample in `Samtools` `pileup` format and outputs an indel file and snp file. The command
reports germline, somatic, and LOH events at positions where both normal and tumor samples have sufficient coverage
### Installation
```
wget https://sourceforge.net/projects/varscan/files/VarScan.v2.3.9.jar		
```
### Command
##### Convert data to pileup
```
mkdir Variants

for sample in `cat list.txt`
do
        samtools mpileup -f hg19.chr5_12_17.fa Mapping/${sample}.refilter.bam --min-MQ 1 --min-BQ 28 \
                > Variants/${sample}.pileup
done
```
##### Call Variants
```
java -jar VarScan.v2.3.9.jar somatic Variants/SLGFSK-N_231335.pileup \
        Variants/SLGFSK-T_231336.pileup Variants/SLGFSK \
        --normal-purity 1  --tumor-purity 0.5 --output-vcf 1 
```
##### Merge vcf
VarScan generates 2 outputs (indel.vcf and snp.vcf), merge the two into one vcf file using `bcftools`.
```
#merge vcf
bgzip Variants/SLGFSK.snp.vcf > Variants/SLGFSK.snp.vcf.gz
bgzip Variants/SLGFSK.indel.vcf > Variants/SLGFSK.indel.vcf.gz
tabix Variants/SLGFSK.snp.vcf.gz
tabix Variants/SLGFSK.indel.vcf.gz
bcftools merge Variants/SLGFSK.snp.vcf.gz Variants/SLGFSK.indel.vcf.gz > Variants/SLGFSK.vcf
```
### Variant Annotation
##### Functional annotation using `SnpEff`
https://pcingola.github.io/SnpEff/examples/
##### Description
`SnpEff` is a variant annotator and functional effect predictor. The output is appended to the vcf file with the field `ANN`. 
A SnpEff database is required prior to performing annotation. In case the organism of interest is not present in the snpEff 
database, you can build the database using the snpEff command. If the organism is present in the database, download it using 
the SnpEff command.
##### Installation
```
#download jar file
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip

# Unzip file
unzip snpEff_latest_core.zip
		
#download snpEff database
java -jar snpEff.jar download hg19
```
##### Command
```
#annotate variants
java -Xmx8g -jar snpEff/snpEff.jar hg19 Variants/SLGFSK.vcf > Variants/SLGFSK.ann.vcf
```		
##### Clinical Annotation using `GEMINI`
https://gemini.readthedocs.io/en/latest/content/preprocessing.html
##### Description
##### Installation
```
wget https://raw.github.com/arq5x/gemini/master/gemini/scripts/gemini_install.py
python gemini_install.py /usr/local /usr/local/share/gemini
```
##### Command
```
gemini load -v Variants/SLGFSK.ann.vcf -t snpEff Annotation/gemini.db
```
Genimi annotation required additional download of annotation resources and databases that totaled to 21GB. Due to lack
of enough resources (bundles and Storage space) to download these resources, we rather reverted to upload the generated 
VCFs into Galaxy and perfrmed the rest of the Annotation with Gemini there with steps described under the Galaxy workflow 
section.
#####Conclusion

##### Contributor
