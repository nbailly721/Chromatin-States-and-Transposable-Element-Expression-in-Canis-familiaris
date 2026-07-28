# Chromatin-States-and-Transposable-Element-Expression-in-Canis-familiaris

**Overview**

This repository contains the computational workflow used to investigate the association between TE expression and chromatin states across four canine tissues (ovary, spleen, cerebrum, cerebellum).

**Repository Structure**

- "scripts" folder:
  
  -te_quantification.sh: Bash script used to calculate the TE counts in every tissue and their associated genomic coordinates.
  
  -chromatin_state_preparation.sh: Bash script used to download and process the chromatin state annotation files so that their genomic coordinates can be intersected with those of the TE counts.
  
  -te_chromatin_intersection.sh: Bash script used to associate normalized TE counts with chromatin states based on overlapping genomic coordinates.
  
  -te_expression_analysis.R: R script used to (1) normalize TE counts prior to chromatin state integration and (2) perform downstream statistical analyses and visualizations of TE expression across chromatin states.

-"intersection_files" folder: Contains the tissue-specific TE–chromatin state overlap tables used for the downstream statistical analyses.

**Datasets Used**

1. Four tissue-specific .bed files containing the genomic coordinates of 13 chromatin states in canine ovary, spleen, cerebrum, and cerebellum.

2. Eight RNA-seq datasets from canine ovary, spleen, cerebrum, and cerebellum. Two biological replicates were analyzed for each tissue. Except for the ovary, each tissue included one male and one female sample.

3. 

**Workflow**
