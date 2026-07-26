#!/bin/bash

#_Load required packages -------------------
module load sra-toolkit/3.0.9
module load fastp
module load fastqc/0.12.1
module load star/2.7.11b
module load  samtools/1.22.1

#_Download RNA-seq data -------------------

prefetch SRR19225570 #Spleen
prefetch SRR19225571
prefetch SRR19225574 #Ovary
prefetch SRR19225575
prefetch SRR19225586 #Cerebrum
prefetch SRR19225587
prefetch SRR19225588 #Cerebellum
prefetch SRR19225589
#Download the raw RNA-seq datasets required for quality control, alignment, and TE quantification.

fasterq-dump --split-files -O fastq_files SRR19225570/SRR19225570.sra
fasterq-dump --split-files -O fastq_files SRR19225571/SRR19225571.sra
fasterq-dump --split-files -O fastq_files SRR19225574/SRR19225574.sra
fasterq-dump --split-files -O fastq_files SRR19225575/SRR19225575.sra
fasterq-dump --split-files -O fastq_files SRR19225586/SRR19225586.sra
fasterq-dump --split-files -O fastq_files SRR19225587/SRR19225587.sra
fasterq-dump --split-files -O fastq_files SRR19225588/SRR19225588.sra
fasterq-dump --split-files -O fastq_files SRR19225589/SRR19225589.sra
#Convert downloaded SRA files into paired-end FASTQ files required for downstream analyses.

#_Quality control -------------------

mkdir -p fastqc_reports
fastqc *.fastq -o fastqc_reports
#Produce quality control reports on the FASTQ files to verify read quality and determine whether trimming is required. It was identified that certain samples had a poly-G tail that required trimming

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
#The poly-G tail of all samples was trimmed for consistency.

mkdir -p trimmed_fastqc_reports
fastqc SRR19225570_1_trimmed.fastq SRR19225570_2_trimmed.fastq -o trimmed_fastqc_reports
fastqc SRR19225571_1_trimmed.fastq SRR19225571_2_trimmed.fastq -o trimmed_fastqc_reports
fastqc SRR19225574_1_trimmed.fastq SRR19225574_2_trimmed.fastq -o trimmed_fastqc_reports
fastqc SRR19225575_1_trimmed.fastq SRR19225575_2_trimmed.fastq -o trimmed_fastqc_reports
fastqc SRR19225586_1_trimmed.fastq SRR19225586_2_trimmed.fastq -o trimmed_fastqc_reports
fastqc SRR19225587_1_trimmed.fastq SRR19225587_2_trimmed.fastq -o trimmed_fastqc_reports
fastqc SRR19225588_1_trimmed.fastq SRR19225588_2_trimmed.fastq -o trimmed_fastqc_reports
fastqc SRR19225589_1_trimmed.fastq SRR19225589_2_trimmed.fastq -o trimmed_fastqc_reports
#Generate post-trimming quality control reports to verify that poly-G trimming improved read quality.

#_STAR index generation -------------------

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/014/441/545/GCF_014441545.1_ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fna.gz
gunzip GCF_014441545.1_ROS_Cfam_1.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/014/441/545/GCF_014441545.1_ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.gtf.gz
gunzip GCF_014441545.1_ROS_Cfam_1.0_genomic.gtf.gz
#Download the canine reference genome (FNA) and gene annotation (GTF) files required for STAR genome indexing.

mkdir -p STAR_index
STAR \
--runThreadN 8 \
--runMode genomeGenerate \
--genomeDir STAR_index \
--genomeFastaFiles GCF_014441545.1_ROS_Cfam_1.0_genomic.fna \
--sjdbGTFfile GCF_014441545.1_ROS_Cfam_1.0_genomic.gtf \
--sjdbOverhang 150
# Generate the STAR genome index required to align RNA-seq reads, producing the BAM files used by Telescope for TE quantification.

#_RNA-seq alignment -------------------

mkdir -p star_alignments

STAR --runThreadN 8 --genomeDir STAR_index --readFilesIn SRR19225570_1_trimmed.fastq SRR19225570_2_trimmed.fastq --outFileNamePrefix star_alignments/SRR19225570_ --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 500
STAR --runThreadN 8 --genomeDir STAR_index --readFilesIn SRR19225571_1_trimmed.fastq SRR19225571_2_trimmed.fastq --outFileNamePrefix star_alignments/SRR19225571_ --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 500
STAR --runThreadN 8 --genomeDir STAR_index --readFilesIn SRR19225574_1_trimmed.fastq SRR19225574_2_trimmed.fastq --outFileNamePrefix star_alignments/SRR19225574_ --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 500
STAR --runThreadN 8 --genomeDir STAR_index --readFilesIn SRR19225575_1_trimmed.fastq SRR19225575_2_trimmed.fastq --outFileNamePrefix star_alignments/SRR19225575_ --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 500
STAR --runThreadN 8 --genomeDir STAR_index --readFilesIn SRR19225586_1_trimmed.fastq SRR19225586_2_trimmed.fastq --outFileNamePrefix star_alignments/SRR19225586_ --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 500
STAR --runThreadN 8 --genomeDir STAR_index --readFilesIn SRR19225587_1_trimmed.fastq SRR19225587_2_trimmed.fastq --outFileNamePrefix star_alignments/SRR19225587_ --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 500
STAR --runThreadN 8 --genomeDir STAR_index --readFilesIn SRR19225588_1_trimmed.fastq SRR19225588_2_trimmed.fastq --outFileNamePrefix star_alignments/SRR19225588_ --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 500
STAR --runThreadN 8 --genomeDir STAR_index --readFilesIn SRR19225589_1_trimmed.fastq SRR19225589_2_trimmed.fastq --outFileNamePrefix star_alignments/SRR19225589_ --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 500
#Align trimmed RNA-seq reads to the canine reference genome, generating sorted BAM files for downstream TE quantification with Telescope.
#The maximum number of multiple alignments per read was set at 500 because transposable elements are repetitive genomic sequences.

#_BAM file indexing -------------------

samtools index SRR19225570_Aligned.sortedByCoord.out.bam
samtools index SRR19225571_Aligned.sortedByCoord.out.bam
samtools index SRR19225574_Aligned.sortedByCoord.out.bam
samtools index SRR19225575_Aligned.sortedByCoord.out.bam
samtools index SRR19225586_Aligned.sortedByCoord.out.bam
samtools index SRR19225587_Aligned.sortedByCoord.out.bam
samtools index SRR19225588_Aligned.sortedByCoord.out.bam
samtools index SRR19225589_Aligned.sortedByCoord.out.bam
#Generate .bai files from .bam files to enable efficient access during downstream TE quantification with Telescope.

#_Prepare Telescope TE annotation -------------------

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/014/441/545/GCF_014441545.1_ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_rm.out.gz  
#Download RepeatMasker annotation file containing transposable element coordinates for the canine reference genome. 
#Necessary to generate the TE annotation file required by Telescope for TE quantification.

awk '
NR > 3 {
    strand = ($9 == "C") ? "-" : "+";

    split($11, arr, "/");
    class = arr[1];
    family = (arr[2] == "" ? arr[1] : arr[2]);

    te_id = $10 "_" $15;

    print $5 "\tRepeatMasker\texon\t" \
          $6 "\t" $7 "\t.\t" strand "\t.\t" \
          "gene_id \"" $10 "\"; transcript_id \"" te_id "\"; class_id \"" class "\"; family_id \"" family "\";"
}
' GCF_014441545.1_ROS_Cfam_1.0_rm.out > ROS_Cfam_TE.gtf
# Convert RepeatMasker TE coordinates into a Telescope-compatible GTF file with transcript, class, and family annotations.

#_TE quantification via Telescope -------------------

source ~/miniconda3/etc/profile.d/conda.sh
conda activate telescope39

telescope assign \
  --attribute transcript_id \
  --outdir telescope_SRR19225570 \
  --exp_tag SRR19225570 \
  SRR19225570_Aligned.sortedByCoord.out.bam \
  ROS_Cfam_TE.gtf

telescope assign \
  --attribute transcript_id \
  --outdir telescope_SRR19225571 \
  --exp_tag SRR19225571 \
  SRR19225571_Aligned.sortedByCoord.out.bam \
  ROS_Cfam_TE.gtf

telescope assign \
  --attribute transcript_id \
  --outdir telescope_SRR19225574 \
  --exp_tag SRR19225574 \
  SRR19225574_Aligned.sortedByCoord.out.bam \
  ROS_Cfam_TE.gtf

  telescope assign \
  --attribute transcript_id \
  --outdir telescope_SRR19225575 \
  --exp_tag SRR19225575 \
  SRR19225575_Aligned.sortedByCoord.out.bam \
  ROS_Cfam_TE.gtf

  telescope assign \
  --attribute transcript_id \
  --outdir telescope_SRR19225586 \
  --exp_tag SRR19225586 \
  SRR19225586_Aligned.sortedByCoord.out.bam \
  ROS_Cfam_TE.gtf

  telescope assign \
  --attribute transcript_id \
  --outdir telescope_SRR19225587 \
  --exp_tag SRR19225587 \
  SRR19225587_Aligned.sortedByCoord.out.bam \
  ROS_Cfam_TE.gtf

telescope assign \
  --attribute transcript_id \
  --outdir telescope_SRR19225588 \
  --exp_tag SRR19225588 \
  SRR19225588_Aligned.sortedByCoord.out.bam \
  ROS_Cfam_TE.gtf

  telescope assign \
  --attribute transcript_id \
  --outdir telescope_SRR19225589 \
  --exp_tag SRR19225589 \
  SRR19225589_Aligned.sortedByCoord.out.bam \
  ROS_Cfam_TE.gtf
  #Produce a transcript-level Telescope report (.tsv) for each RNA-seq sample. 
  #Each report was imported into RStudio, where TE counts were merged across samples and normalized.

  
