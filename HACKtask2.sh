#!/bin/bash
#Download dataset
echo -e "\n Downloading data... \n"   
mkdir -p raw_data 
cd raw_data
        
wget https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r1_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r2_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r1_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r2_chr5_12_17.fastq.gz

#Reference sequence
echo -e "\n Downloading reference sequence... \n"
        
wget https://zenodo.org/record/2582555/files/hg19.chr5_12_17.fa.gz

#unzip reference
gunzip hg19.chr5_12_17.fa.gz

#Preprocessing and trimming
echo -e "\n Data Preprocessing... \n"

mkdir -p Fastqc_Reports  #creates directory for the fastqc output
fastqc *.fastq.gz -o Fastqc_Reports
multiqc Fastqc_Reports -o Fastqc_Reports
#Removing low quality reads with trimmomatic
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

multiqc  Fastqc_results -o Fastqc_results

#Mapping
#Index reference file   
bwa index hg19.chr5_12_17.fa 
mkdir Mapping
   
#Perform alignment
bwa mem -R '@RG\tID:231335\tSM:Normal' hg19.chr5_12_17.fa trimmed_reads/SLGFSK-N_231335_r1_paired.fq.gz \
      trimmed_reads/SLGFSK-N_231335_r2_paired.fq.gz > Mapping/SLGFSK-N_231335.sam

bwa mem -R '@RG\tID:231336\tSM:Tumor' hg19.chr5_12_17.fa trimmed_reads/SLGFSK-T_231336_r1_paired.fq.gz \
       trimmed_reads/SLGFSK-T_231336_r2_paired.fq.gz > Mapping/SLGFSK-T_231336.sam      

#convert samfiles to bam files, sort it and index bam file
for sample in `cat list.txt`
do     
        samtools view -@ 20 -S -b Mapping/${sample}.sam | samtools sort -@ 32 > Mapping/${sample}.sorted.bam
        samtools index Mapping/${sample}.sorted.bam
done

#Filter mapped reads
for sample in `cat list.txt`
do
        samtools view -q 1 -f 0x2 -F 0x8 -b Mapping/${sample}.sorted.bam > Mapping/${sample}.filtered1.bam | samtools flagstat ${sample}.filtered1.bam
done

#Duplicate removal using markdup
for sample in `cat list.txt`
do
        samtools collate -o Mapping/${sample}.filtered1.bam Mapping/${sample}.namecollate.bam
        samtools fixmate -m Mapping/${sample}.namecollate.bam Mapping/${sample}.fixmate.bam
        samtools sort -@ 32 -o Mapping/${sample}.positionsort.bam Mapping/${sample}.fixmate.bam
        samtools markdup -@32 -r Mapping/${sample}.positionsort.bam Mapping/${sample}.clean.bam
 done
     #Left Align Bam   
for sample in `cat list.txt`
do      
     cat Mapping/${sample}.clean.bam | bamleftalign -f hg19.chr5_12_17.fa -m 5 -c > Mapping/${sample}.leftAlign.bam
done
#recalibrate read mapping qualities
for sample in `cat list.txt`
do
        samtools calmd -@ 32 -b Mapping/${sample}.leftAlign.bam hg19.chr5_12_17.fa > Mapping/${sample}.recalibrate.bam
done
#Refilter reads
for sample in `cat list.txt`
do
        bamtools filter -in Mapping/${sample}.recalibrate.bam -mapQuality "<=254" > Mapping/${sample}.refilter.bam
done

#Variant calling and classification

wget https://sourceforge.net/projects/varscan/files/VarScan.v2.3.9.jar

# convert data to pileup
mkdir Variants

for sample in `cat list.txt`
do    
       samtools mpileup -f hg19.chr5_12_17.fa Mapping/${sample}.refilter.bam --min-MQ 1 --min-BQ 28 \
       > Variants/${sample}.pileup
 done

# call variant
      java -jar VarScan.v2.3.9.jar somatic Variants/SLGFSK-N_231335.pileup \
      Variants/SLGFSK-T_231336.pileup Variants/SLGFSK \ 
      --normal-purity 1  --tumor-purity 0.5 --output-vcf 1 
```

#merge vcf
bgzip Variants/SLGFSK.snp.vcf > Variants/SLGFSK.snp.vcf.gz
bgzip Variants/SLGFSK.indel.vcf > Variants/SLGFSK.indel.vcf.gz
tabix Variants/SLGFSK.snp.vcf.gz
tabix Variants/SLGFSK.indel.vcf.gz
bcftools merge Variants/SLGFSK.snp.vcf.gz Variants/SLGFSK.indel.vcf.gz > Variants/SLGFSK.vcf

#Functional annotation using SnpEff
#download jar file
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip

# Unzip file
unzip snpEff_latest_core.zip
                
#download snpEff database
java -jar snpEff.jar download hg19

#annotate variants
java -Xmx8g -jar snpEff/snpEff.jar hg19 Variants/SLGFSK.vcf > Variants/SLGFSK.ann.vcf
