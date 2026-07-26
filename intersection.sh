#!/bin/bash
module load bedtools

awk -F'\t' '
{
    split($9,a,";")
    gsub(/transcript_id "/,"",a[2])
    gsub(/"/,"",a[2])
    gsub(/^ /,"",a[2])

    id=a[2]

    if(!(id in start) || $4<start[id]) start[id]=$4
    if(!(id in end) || $5>end[id]) end[id]=$5

    chr[id]=$1
}
END{
    for(i in chr)
        print chr[i]"\t"start[i]"\t"end[i]"\t"i
}
' ROS_Cfam_TE.gtf > TE_coordinates.bed
#Retrieve genomic coordinates for each transposable element from the Telescope annotation to create a BED file for downstream overlap analyses.

awk -F'\t' '
NR==FNR {
    if (FNR > 1)
        data[$1] = $0
    next
}
($4 in data) {
    print $0 "\t" data[$4]
}
' Telescope_counts.tsv TE_coordinates.bed > TE_expression.bed

awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$6,$7,$8,$9,$10,$11,$12,$13}' TE_expression.bed > TE_expression_clean.bed
# Merge normalized TE expression values with TE genomic coordinates and retain only the columns required for downstream chromatin state overlap analysis.

#_TE–chromatin state overlap analysis -------------------

bedtools intersect \
-a TE_expression_clean.bed \
-b OV_ROS.bed \
-wa -wb > OV_intersect.tsv

bedtools intersect \
-a TE_expression_clean.bed \
-b SP_ROS.bed \
-wa -wb > SP_intersect.tsv

bedtools intersect \
-a TE_expression_clean.bed \
-b CR_ROS.bed \
-wa -wb > CR_intersect.tsv

bedtools intersect \
-a TE_expression_clean.bed \
-b CL_ROS.bed \
-wa -wb > CL_intersect.tsv
#The produced .tsv files contained the normalized TE counts associated with each chromatin state based on their overlapping genomic coordinates. They were exported into Rstudio for downstream analysis.
