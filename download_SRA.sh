#!/bin/bash
module load sra-toolkit/3.0.9
module load fastp

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

