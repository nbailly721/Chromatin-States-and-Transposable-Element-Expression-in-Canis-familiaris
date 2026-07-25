#!/bin/bash
module load sra-toolkit/3.0.9

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

