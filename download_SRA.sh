#!/bin/bash
module load sra-toolkit/3.0.9
module load fastp
module load fastqc/0.12.1
module load star/2.7.11b

prefetch SRR19225570 #Spleen
prefetch SRR19225571

prefetch SRR19225574 #Ovary
prefetch SRR19225575

prefetch SRR19225586 #Cerebrum
prefetch SRR19225587

prefetch SRR19225588 #Cerebellum
prefetch SRR19225589

fasterq-dump --split-files -O fastq_files SRR19225570/SRR19225570.sra
fasterq-dump --split-files -O fastq_files SRR19225571/SRR19225571.sra

fasterq-dump --split-files -O fastq_files SRR19225574/SRR19225574.sra
fasterq-dump --split-files -O fastq_files SRR19225575/SRR19225575.sra

fasterq-dump --split-files -O fastq_files SRR19225586/SRR19225586.sra
fasterq-dump --split-files -O fastq_files SRR19225587/SRR19225587.sra

fasterq-dump --split-files -O fastq_files SRR19225588/SRR19225588.sra
fasterq-dump --split-files -O fastq_files SRR19225589/SRR19225589.sra

fastp \
-i SRR19225571_1.fastq \
-I SRR19225571_2.fastq \
-o SRR19225571_1_trimmed.fastq \
-O SRR19225571_2_trimmed.fastq \
--trim_poly_g

fastp \
-i SRR19225570_1.fastq \
-I SRR19225570_2.fastq \
-o SRR19225570_1_trimmed.fastq \
-O SRR19225570_2_trimmed.fastq \
--trim_poly_g

fastp \
-i SRR19225574_1.fastq \
-I SRR19225574_2.fastq \
-o SRR19225574_1_trimmed.fastq \
-O SRR19225574_2_trimmed.fastq \
--trim_poly_g

fastp \
-i SRR19225575_1.fastq \
-I SRR19225575_2.fastq \
-o SRR19225575_1_trimmed.fastq \
-O SRR19225575_2_trimmed.fastq \
--trim_poly_g

fastp \
-i SRR19225586_1.fastq \
-I SRR19225586_2.fastq \
-o SRR19225586_1_trimmed.fastq \
-O SRR19225586_2_trimmed.fastq \
--trim_poly_g

fastp \
-i SRR19225587_1.fastq \
-I SRR19225587_2.fastq \
-o SRR19225587_1_trimmed.fastq \
-O SRR19225587_2_trimmed.fastq \
--trim_poly_g

fastp \
-i SRR19225588_1.fastq \
-I SRR19225588_2.fastq \
-o SRR19225588_1_trimmed.fastq \
-O SRR19225588_2_trimmed.fastq \
--trim_poly_g

fastp \
-i SRR19225589_1.fastq \
-I SRR19225589_2.fastq \
-o SRR19225589_1_trimmed.fastq \
-O SRR19225589_2_trimmed.fastq \
--trim_poly_g

mkdir -p trimmed_fastqc_reports

fastqc SRR19225570_1_trimmed.fastq SRR19225570_2_trimmed.fastq -o trimmed_fastqc_reports

fastqc SRR19225571_1_trimmed.fastq SRR19225571_2_trimmed.fastq -o trimmed_fastqc_reports

fastqc SRR19225574_1_trimmed.fastq SRR19225574_2_trimmed.fastq -o trimmed_fastqc_reports

fastqc SRR19225575_1_trimmed.fastq SRR19225575_2_trimmed.fastq -o trimmed_fastqc_reports

fastqc SRR19225586_1_trimmed.fastq SRR19225586_2_trimmed.fastq -o trimmed_fastqc_reports

fastqc SRR19225587_1_trimmed.fastq SRR19225587_2_trimmed.fastq -o trimmed_fastqc_reports

fastqc SRR19225588_1_trimmed.fastq SRR19225588_2_trimmed.fastq -o trimmed_fastqc_reports

fastqc SRR19225589_1_trimmed.fastq SRR19225589_2_trimmed.fastq -o trimmed_fastqc_reports

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/014/441/545/GCF_014441545.1_ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fna.gz
gunzip GCF_014441545.1_ROS_Cfam_1.0_genomic.fna.gz

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/014/441/545/GCF_014441545.1_ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.gtf.gz
gunzip GCF_014441545.1_ROS_Cfam_1.0_genomic.gtf.gz

mkdir -p STAR_index

STAR \
--runThreadN 8 \
--runMode genomeGenerate \
--genomeDir STAR_index \
--genomeFastaFiles GCF_014441545.1_ROS_Cfam_1.0_genomic.fna \
--sjdbGTFfile GCF_014441545.1_ROS_Cfam_1.0_genomic.gtf \
--sjdbOverhang 150




