# Chromatin-States-and-Transposable-Element-Expression-in-Canis-familiaris

**Overview**

This repository contains the computational workflow used to investigate the association between TE expression and chromatin states across four canine tissues (ovary, spleen, cerebrum, cerebellum).

**Repository Structure**

- **scripts/** folder
  
   - **te_quantification.sh: Bash script used to calculate the TE counts in every tissue and their associated genomic coordinates.
  
   - **chromatin_state_preparation.sh: Bash script used to download and process the chromatin state annotation files so that their genomic coordinates can be intersected with those of the TE counts.
  
   - **te_chromatin_intersection.sh: Bash script used to associate normalized TE counts with chromatin states based on overlapping genomic coordinates.
  
   - **te_expression_analysis.R: R script used to (1) normalize TE counts prior to chromatin state integration and (2) perform downstream statistical analyses and visualizations of TE expression across chromatin states.

- **intersection_files/** folder
   - Contains the tissue-specific TE–chromatin state overlap tables used for the downstream statistical analyses.
     
- **environment_telescope.yml** file: Contains the Conda environment used for Telescope-based TE quantification.

**Datasets Used**

1.Four tissue-specific .bed files containing the genomic coordinates of 13 chromatin states in canine ovary, spleen, cerebrum, and cerebellum (Son et al., 2025).

2.Eight RNA-seq datasets from canine ovary, spleen, cerebrum, and cerebellum. Two biological replicates were analyzed for each tissue. Except for the ovary, each tissue included one male and one female sample (Gene Expression Omnibus, 2023).

3.Canis familiaris reference genome FNA file, gene annotation GTF file, and RepeatMasker annotation file (NCBI RefSeq Genome Database).

**Workflow**

1.Run 'te_quantification.sh' to download the RNA-seq datasets, align reads to the reference genome, and quantify Transposable Element expression across the four tissues.

2.Run`chromatin_state_preparation.sh` to download the chromatin state BED files and convert their genomic coordinates to the assembly version of the reference genome.

3.Run `te_expression_analysis.R` to merge and normalize the TE counts across samples. The table containing the normalized and combined counts are exported as `Telescope_counts.tsv`.

4.Run `te_chromatin_intersection.sh' to intersect the genomic coordinates of the chromatin states with those of the normalized TE counts and produce tissue-specific overlap tables.

5.Run `te_expression_analysis.R` again to perform the statistical analyses on the tissue-specific overlap tables and produce relevant visualizations.

**References**

1.Gene Expression Omnibus (2023). Integrated mapping of the dog genome. National Center for Biotechnology Information. Available at: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203107

2.NCBI RefSeq Genome Database. Canis lupus familiaris genome assembly. ROS_Cfam_1.0 (GCF_014441545.1). National Center for Biotechnology Information. 

3.Son, K.H. et al. (2025) EpiC Dog: Epigenome Catalog of the Dog. GitHub repository. Available at: https://github.com/snu-cdrc/dog-reference-epigenome



