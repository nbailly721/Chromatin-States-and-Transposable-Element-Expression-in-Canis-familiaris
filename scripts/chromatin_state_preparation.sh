#!/bin/bash

#_Load required packages -------------------

module load kent_tools/486

#_Download chromatin state BED files -------------

wget "https://raw.githubusercontent.com/snu-cdrc/dog-reference-epigenome/main/Data/04_chrom_states%20(canFam3.1)/CL_13_dense.bed" #cerebellum
wget "https://raw.githubusercontent.com/snu-cdrc/dog-reference-epigenome/main/Data/04_chrom_states%20(canFam3.1)/CR_13_dense.bed" #cerebrum
wget "https://raw.githubusercontent.com/snu-cdrc/dog-reference-epigenome/main/Data/04_chrom_states%20(canFam3.1)/OV_13_dense.bed" #ovary
wget "https://raw.githubusercontent.com/snu-cdrc/dog-reference-epigenome/main/Data/04_chrom_states%20(canFam3.1)/SP_13_dense.bed" #spleen
#Download tissue-specific .bed files containing chromatin state coordinate information for downstream conversion and ovelap analysis

#_Coordinate conversion with liftOver -------------

wget https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/liftOver/canFam3ToGCF_014441545.1.over.chain.gz
gunzip canFam3ToGCF_014441545.1.over.chain.gz
#Download the liftOver chain file required to convert chromatin state coordinates to the GCF_014441545.1_ROS_Cfam_1.0 reference genome. 

liftOver OV_13_dense.bed canFam3ToGCF_014441545.1.over.chain OV_ROS.bed OV_unmapped.bed
liftOver SP_13_dense.bed canFam3ToGCF_014441545.1.over.chain SP_ROS.bed SP_unmapped.bed
liftOver CR_13_dense.bed canFam3ToGCF_014441545.1.over.chain CR_ROS.bed CR_unmapped.bed
liftOver CL_13_dense.bed canFam3ToGCF_014441545.1.over.chain CL_ROS.bed CL_unmapped.bed
#The .bed files were converted to the ROS_Cfam_1.0 genome assembly to ensure coordinate compatibility with the RNA-seq reference genome.

################################################
te_chromatin_intersection.sh script follows
################################################
