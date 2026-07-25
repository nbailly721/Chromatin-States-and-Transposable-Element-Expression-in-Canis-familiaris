#!/bin/bash
module load fastp

fastp \
-i SRR19225571_1.fastq \
-I SRR19225571_2.fastq \
-o .../SRR19225571_1_trimmed.fastq \
-O .../SRR19225571_2_trimmed.fastq \
--trim_poly_g

fastp \
-i SRR19225570_1.fastq \
-I SRR19225570_2.fastq \
-o .../SRR19225570_1_trimmed.fastq \
-O .../SRR19225570_2_trimmed.fastq \
--trim_poly_g

fastp \
-i SRR19225574_1.fastq \
-I SRR19225574_2.fastq \
-o .../SRR19225574_1_trimmed.fastq \
-O .../SRR19225574_2_trimmed.fastq \
--trim_poly_g

fastp \
-i SRR19225575_1.fastq \
-I SRR19225575_2.fastq \
-o .../SRR19225575_1_trimmed.fastq \
-O .../SRR19225575_2_trimmed.fastq \
--trim_poly_g

fastp \
-i SRR19225586_1.fastq \
-I SRR19225586_2.fastq \
-o .../SRR19225586_1_trimmed.fastq \
-O .../SRR19225586_2_trimmed.fastq \
--trim_poly_g

fastp \
-i SRR19225587_1.fastq \
-I SRR19225587_2.fastq \
-o .../SRR19225587_1_trimmed.fastq \
-O .../SRR19225587_2_trimmed.fastq \
--trim_poly_g

fastp \
-i SRR19225588_1.fastq \
-I SRR19225588_2.fastq \
-o .../SRR19225588_1_trimmed.fastq \
-O .../SRR19225588_2_trimmed.fastq \
--trim_poly_g

fastp \
-i SRR19225589_1.fastq \
-I SRR19225589_2.fastq \
-o ...s/SRR19225589_1_trimmed.fastq \
-O .../SRR19225589_2_trimmed.fastq \
--trim_poly_g
