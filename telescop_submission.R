#########################################
#Data analyses and visualization
#########################################

#_Create environment --------------------

install.packages("FSA")
library(DESeq2)
library(dplyr)
install.packages("FSA")
library(FSA)
library(pheatmap)
library(ggplot2)
library(grid)
#Install and load required packages for downstream data analyses

#_Pre-normalization steps --------------------

s_spleen_male <- read.delim("SRR19225570-telescope_report.tsv", skip = 1)
s_spleen_female <- read.delim("SRR19225571-telescope_report.tsv", skip = 1)
s_ovary_1 <- read.delim("SRR19225574-telescope_report.tsv", skip = 1)
s_ovary_2 <- read.delim("SRR19225575-telescope_report.tsv", skip = 1)
s_cerebrum_male <- read.delim("SRR19225586-telescope_report.tsv", skip = 1)
s_cerebrum_female <- read.delim("SRR19225587-telescope_report.tsv", skip = 1)
s_cerebellum_male <- read.delim("SRR19225588-telescope_report.tsv", skip = 1)
s_cerebellum_female <- read.delim("SRR19225589-telescope_report.tsv", skip = 1)
#Load the Telescope reports imported from Bash. These reports contain the raw TE counts and their genomic locations in each tissue sample.

print(s_spleen_male )
colnames(s_spleen_male)
#Ensure that reports have the correct format for downstream analyses.

# Keep only transcript and final_count
s_spleen_male <- s_spleen_male[, c("transcript", "final_count")]
s_spleen_female <- s_spleen_female[, c("transcript", "final_count")]
s_ovary_1 <- s_ovary_1[, c("transcript", "final_count")]
s_ovary_2 <- s_ovary_2[, c("transcript", "final_count")]
s_cerebrum_male <- s_cerebrum_male[, c("transcript", "final_count")]
s_cerebrum_female <- s_cerebrum_female[, c("transcript", "final_count")]
s_cerebellum_male <- s_cerebellum_male[, c("transcript", "final_count")]
s_cerebellum_female <- s_cerebellum_female[, c("transcript", "final_count")]
#Maintain only the columns named 'transcript' and 'final_count'. This step eliminates any redundant information and keeps only the columns of interest.

colnames(s_spleen_male)[2] <- "Spleen_Male"
colnames(s_spleen_female)[2] <- "Spleen_Female"
colnames(s_ovary_1)[2] <- "Ovary_1"
colnames(s_ovary_2)[2] <- "Ovary_2"
colnames(s_cerebrum_male)[2] <- "Cerebrum_Male"
colnames(s_cerebrum_female)[2] <- "Cerebrum_Female"
colnames(s_cerebellum_male)[2] <- "Cerebellum_Male"
colnames(s_cerebellum_female)[2] <- "Cerebellum_Female"
#Assign sample-specific names to final_count columns. Necessary for merging all data sets into a single count matrix

merged_table <- merge(s_spleen_male, s_spleen_female, by = "transcript")
merged_table <- merge(merged_table, s_ovary_1, by = "transcript")
merged_table <- merge(merged_table, s_ovary_2, by = "transcript")
merged_table <- merge(merged_table, s_cerebrum_male, by = "transcript")
merged_table <- merge(merged_table, s_cerebrum_female, by = "transcript")
merged_table <- merge(merged_table, s_cerebellum_male, by = "transcript")
merged_table <- merge(merged_table, s_cerebellum_female, by = "transcript")
#Merge all data sets into a single count matrix. This step is required to create the DESeq2 input matrix for downstream normalization

dim(merged_table)
head(merged_table)
#To ensure that the merging steps worked correctly.

#_Data Normalization --------------------

count_matrix <- merged_table[, -1]
rownames(count_matrix) <- merged_table$transcript
#Extract count data and set transcript IDs as row names for DESeq2 input.

sample_info <- data.frame(
  tissue = c("spleen", "spleen",
             "ovary", "ovary",
             "cerebrum", "cerebrum",
             "cerebellum", "cerebellum")
)
#Label each sample with its respetive tissue type.

rownames(sample_info) <- colnames(count_matrix)
#Match sample metadata to the count matrix columns.

dds <- DESeqDataSetFromMatrix(
  countData = round(count_matrix),
  colData = sample_info,
  design = ~ tissue
)
#Create the DESeq2 dataset object for normalization of raw TE counts.

dds <- estimateSizeFactors(dds)
#Estimate size factors to normalize sequencing depth across samples.

normalized_counts <- counts(dds, normalized = TRUE)
tail(normalized_counts)
#Extract normalized TE counts and revise that they are correctly formatted.

normalized_table <- data.frame(
  transcript = rownames(normalized_counts),
  normalized_counts
)
print(normalized_table)
#Convert normalized counts into a table for export and downstream integration (in Bash) with chromatin state overlap data.

write.csv(
  normalized_table,
  "Telescope_counts.csv",
  row.names = FALSE
)
#Save table containing normalized TE counts.

#_Principal Component Analysis (QC) --------------------

vsd <- vst(dds)
#Transform the data to reduce differences caused by abnormally high expression values.

pca <- plotPCA(vsd, intgroup = "tissue", returnData = TRUE)
#Perform PCA to check whether samples cluster according to tissue type.

pca$sex <- c(
  "Male", "Female",
  "Female", "Female",
  "Male", "Female",
  "Male", "Female"
)
#Add sex-related information to each sample. Necessary to identify if sex may play a role in clustering patters.

percentVar <- round(100 * attr(pca, "percentVar"))
#Calculates the percentage of variance explained by the first two principal components.

ggplot(pca, aes(PC1, PC2, color = tissue, shape = sex)) +
  geom_point(size = 4) +
  scale_shape_manual(values = c(
    "Female" = 15,
    "Male" = 17
  )) +
  scale_color_manual(values = c(
    "ovary" = "#009E73",
    "spleen" = "#CC79A7",
    "cerebrum" = "#0072B2",
    "cerebellum" = "#D55E00"
  )) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of TE Expression Across Four Canine Tissues") +
  theme(
    plot.title = element_text(hjust = 0.5)
  )
#Generate a PCA plot to visualize sample clustering based on normalized TE expression.

#_Pre-analysis steps ------------------

ovary <- read.delim("OV_intersect.tsv", header = FALSE)
spleen <- read.delim("SP_intersect.tsv", header = FALSE)
cerebrum <- read.delim("CR_intersect.tsv", header = FALSE)
cerebellum <- read.delim("CL_intersect.tsv", header = FALSE)
#Import chromatin state–TE overlap data for downstream statistical analyses. 

col_names <- c(
  "Chr",
  "Start",
  "End",
  "Transcript",
  "Spleen_Male",
  "Spleen_Female",
  "Ovary_1",
  "Ovary_2",
  "Cerebrum_Male",
  "Cerebrum_Female",
  "Cerebellum_Male",
  "Cerebellum_Female",
  "Column13",
  "Column14",
  "Column15",
  "State",
  "Column17",
  "Column18",
  "Column19",
  "Column20",
  "Column21"
)
colnames(ovary) <- col_names
colnames(spleen) <- col_names
colnames(cerebrum) <- col_names
colnames(cerebellum) <- col_names
#Assign informative column names to identify tissue type, sample sex, and chromatin states.

ovary <- ovary[, c("Transcript", "Ovary_1", "Ovary_2", "State")]
spleen <- spleen[, c("Transcript", "Spleen_Male", "Spleen_Female", "State")]
cerebrum <- cerebrum[, c("Transcript", "Cerebrum_Male", "Cerebrum_Female", "State")]
cerebellum <- cerebellum[, c("Transcript", "Cerebellum_Male", "Cerebellum_Female", "State")]
#Retain only the transcript ID, sample expression values, and chromatin state columns for downstream analyses. The other columns are redundant for the analyses

#_Expression Analysis --------------

##__Descriptive statistics --------------

ovary$Mean_Expression <- rowMeans(ovary[, c("Ovary_1", "Ovary_2")])
spleen$Mean_Expression <- rowMeans(spleen[, c("Spleen_Male", "Spleen_Female")])
cerebrum$Mean_Expression <- rowMeans(cerebrum[, c("Cerebrum_Male", "Cerebrum_Female")])
cerebellum$Mean_Expression <- rowMeans(cerebellum[, c("Cerebellum_Male", "Cerebellum_Female")])
#Calculate the mean TE expression across biological replicates for each tissue.

##__Visualizations --------------

###___Boxplots --------------

state_names <- c(
  "1 Active TSS",
  "2 Weak TSS",
  "3 Flanking TSS1",
  "4 Flanking TSS2",
  "5 Strong Enhancer",
  "6 Weak Enhancer",
  "7 Poised Enhancer",
  "8 Bivalent TSS/Enh",
  "9 Repressed Polycomb",
  "10 Repressed",
  "11 ZNF genes/Repeats",
  "12 Heterochromatin",
  "13 Quiescent"
)
#Assign descriptive labels to each chromatin state in the legend. Each label was retrieved from the GitHub repository where the Bed files of the chromatin states used for analyses were located.

cols <- rainbow(13)
#Assign a distinct color to each of the chromatin states present in the legend.

par(mar = c(5, 4, 4, 10), xpd = TRUE)
#Increase the right margin of the plot so that the legend could fit.

boxplot(
  Mean_Expression ~ State,
  data = ovary,
  xlab = "Chromatin State",
  ylab = "Mean Normalized TE counts",
  main = "Distribution of TE Expression Across Chromatin States in the Ovary",
  outline = FALSE,
  col = cols
)
legend(
  "topright",
  inset = c(-0.35, 0),
  legend = state_names,
  fill = cols,
  bty = "n",
  cex = 0.8
)
#Boxplot compating the distribution of TE expression across chromatin states for each tissue in the ovary.

boxplot(
  Mean_Expression ~ State,
  data = spleen,
  xlab = "Chromatin State",
  ylab = "Mean Normalized Expression",
  main = "Normalized TE expression across chromatin states in the Spleen",
  outline = FALSE,
  col = cols
)
legend(
  "topright",
  inset = c(-0.35, 0),
  legend = state_names,
  fill = cols,
  bty = "n",
  cex = 0.8
)
#Boxplot compating the distribution of TE expression across chromatin states for each tissue in the spleen

boxplot(
  Mean_Expression ~ State,
  data = cerebrum,
  xlab = "Chromatin State",
  ylab = "Mean Normalized Expression",
  main = "Normalized TE expression across chromatin states in the Cerebrum",
  outline = FALSE,
  col = cols
)
legend(
  "topright",
  inset = c(-0.35, 0),
  legend = state_names,
  fill = cols,
  bty = "n",
  cex = 0.8
)
#Boxplot compating the distribution of TE expression across chromatin states for each tissue in the cerebrum

boxplot(
  Mean_Expression ~ State,
  data = cerebellum,
  xlab = "Chromatin State",
  ylab = "Mean Normalized Expression",
  main = "Normalized TE expression across chromatin states in the Cerebellum",
  outline = FALSE,
  col = cols
)
legend(
  "topright",
  inset = c(-0.35, 0),
  legend = state_names,
  fill = cols,
  bty = "n",
  cex = 0.8
)
#Boxplot compating the distribution of TE expression across chromatin states for each tissue in the cerebellum

par(mar = c(5.1, 4.1, 4.1, 2.1))
#Return margins to their original parameters.

###___Heatmap --------------

ovary_heat <- aggregate(Mean_Expression ~ State, data = ovary, mean)
spleen_heat <- aggregate(Mean_Expression ~ State, data = spleen, mean)
cerebrum_heat <- aggregate(Mean_Expression ~ State, data = cerebrum, mean)
cerebellum_heat <- aggregate(Mean_Expression ~ State, data = cerebellum, mean)
#Calculate the mean TE expression for each chromatin state within each tissue. Necessary for the downstream generation of the heatmap

colnames(ovary_heat)[2] <- "Ovary"
colnames(spleen_heat)[2] <- "Spleen"
colnames(cerebrum_heat)[2] <- "Cerebrum"
colnames(cerebellum_heat)[2] <- "Cerebellum"
#Rename the mean expression columns by tissue for easier identification after merging.

heatmap_data <- merge(ovary_heat, spleen_heat, by = "State")
heatmap_data <- merge(heatmap_data, cerebrum_heat, by = "State")
heatmap_data <- merge(heatmap_data, cerebellum_heat, by = "State")
#Merge tissue-specific datasets by chromatin state to create a single matrix. 

heat_matrix <- as.matrix(heatmap_data[, -1])
#Convert the merged table into a matrix for heatmap generation.

state_names <- c(
  "1 Active TSS",
  "2 Weak TSS",
  "3 Flanking TSS1",
  "4 Flanking TSS2",
  "5 Strong Enhancer",
  "6 Weak Enhancer",
  "7 Poised Enhancer",
  "8 Bivalent TSS/Enh",
  "9 Repressed Polycomb",
  "10 Repressed",
  "11 ZNF genes/Repeats",
  "12 Heterochromatin",
  "13 Quiescent"
)
#Assign descriptive labels to every chromatin state. Each label was retrieved from the GitHub repository where the Bed files of the chromatin states used for analyses were located.

rownames(heat_matrix) <- state_names
#Set chromatin state labels as the heatmap row names.

heat_matrix_log <- log2(heat_matrix + 1)
#Log-transform the data displayed in the heatmap to improve visualization by reducing the influence of highly expressed

pheatmap(
  heat_matrix_log,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = FALSE,
  color = colorRampPalette(c("white", "#9ECAE1", "#08519C"))(100),
  main = "Mean TE Expression Across Chromatin States in Four Canine Tissues"
)
#Generate a heatmap to visualize tissue-specific patterns of TE expression across chromatin states. 

##__Within-state variability ----------

ovary_summary <- ovary %>%
  group_by(State) %>%
  summarise(
    Mean = mean(Mean_Expression),
    Median = median(Mean_Expression),
    Min = min(Mean_Expression),
    Max = max(Mean_Expression),
    SD = sd(Mean_Expression),
    n = n()
  )
#Calculate summary statistics to identify within-state variability in TE expression for each chromatin state in the ovary.

spleen_summary <- spleen %>%
  group_by(State) %>%
  summarise(
    Mean = mean(Mean_Expression),
    Median = median(Mean_Expression),
    Min = min(Mean_Expression),
    Max = max(Mean_Expression),
    SD = sd(Mean_Expression),
    n = n()
  )
#Calculate summary statistics to identify within-state variability in TE expression for each chromatin state in the spleen

cerebrum_summary <- cerebrum %>%
  group_by(State) %>%
  summarise(
    Mean = mean(Mean_Expression),
    Median = median(Mean_Expression),
    Min = min(Mean_Expression),
    Max = max(Mean_Expression),
    SD = sd(Mean_Expression),
    n = n()
  )
#Calculate summary statistics to identify within-state variability in TE expression for each chromatin state in the cerebrum

cerebellum_summary <- cerebellum %>%
  group_by(State) %>%
  summarise(
    Mean = mean(Mean_Expression),
    Median = median(Mean_Expression),
    Min = min(Mean_Expression),
    Max = max(Mean_Expression),
    SD = sd(Mean_Expression),
    n = n()
  )
#Calculate summary statistics to identify within-state variability in TE expression for each chromatin state in the cerebellum

write.csv(ovary_summary, "Ovary_summary.csv", row.names = FALSE)
write.csv(spleen_summary, "Spleen_summary.csv", row.names = FALSE)
write.csv(cerebrum_summary, "Cerebrum_summary.csv", row.names = FALSE)
write.csv(cerebellum_summary, "Cerebellum_summary.csv", row.names = FALSE)
#Export the tables containing the summary statistics for a better visualization

#_Statistical tests --------------

##__Kruskal–Wallis tests ----------

kruskal.test(Mean_Expression ~ State, data = ovary)
#To test whether TE expression differs significantly among chromatin states within each tissue. Significant differences were found within the ovary.

kruskal.test(Mean_Expression ~ State, data = spleen)
#To test whether TE expression differs significantly among chromatin states within each tissue. Significant differences were found within the spleen. 

kruskal.test(Mean_Expression ~ State, data = cerebrum)
#To test whether TE expression differs significantly among chromatin states within each tissue. Significant differences were found within the cerebrum. 

kruskal.test(Mean_Expression ~ State, data = cerebellum)
#To test whether TE expression differs significantly among chromatin states within each tissue. Significant differences were found within the cerebellum. 

##__Dunn's post hoc tests ----------

dunn_ovary <- dunnTest(Mean_Expression ~ State,
                       data = ovary,
                       method = "bh")
#Pairwise comparisons to determine which specific chromatin states differ significantly in the ovary. The Benjamini–Hochberg option was chosen to adjust p-values for multiple testing.

dunn_spleen <- dunnTest(Mean_Expression ~ State,
                        data = spleen,
                        method = "bh")
#Pairwise comparisons to determine which specific chromatin states differ significantly in the spleen. The Benjamini–Hochberg option was chosen to adjust p-values for multiple testing.

dunn_cerebrum <- dunnTest(Mean_Expression ~ State,
                          data = cerebrum,
                          method = "bh")
#Pairwise comparisons to determine which specific chromatin states differ significantly in the cerebrum. The Benjamini–Hochberg option was chosen to adjust p-values for multiple testing.

dunn_cerebellum <- dunnTest(Mean_Expression ~ State,
                            data = cerebellum,
                            method = "bh")
#Pairwise comparisons to determine which specific chromatin states differ significantly in the cerebrum. The Benjamini–Hochberg option was chosen to adjust p-values for multiple testing.

dunn_cerebellum$res %>%
  filter(grepl("^1 -| - 1$", Comparison),
         P.adj < 0.05) %>%
  arrange(P.adj)
#Extract chromatin states that differ significantly (adjusted p < 0.05) from the chromatin state with the highest mean TE expression in the cerebellum.

dunn_ovary$res %>%
  filter(grepl("^2 -| - 2$", Comparison),
         P.adj < 0.05) %>%
  arrange(P.adj)
#Extract chromatin states that differ significantly (adjusted p < 0.05) from the chromatin state with the highest mean TE expression in the ovary.

dunn_spleen$res %>%
  filter(grepl("^9 -| - 9$", Comparison),
         P.adj < 0.05) %>%
  arrange(P.adj)
#Extract chromatin states that differ significantly (adjusted p < 0.05) from the chromatin state with the highest mean TE expression in the spleen.

dunn_cerebrum$res %>%
  filter(grepl("^10 -| - 10$", Comparison),
         P.adj < 0.05) %>%
  arrange(P.adj)
#Extract chromatin states that differ significantly (adjusted p < 0.05) from the chromatin state with the highest mean TE expression in the cerebrum.


