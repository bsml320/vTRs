#!/bin/bash
###
 # @author Citu
 # @email Citu.Citu@uth.tmc.edu
 # @create date 2024-12-28 23:16:53

# Define paths to various tools and files
DIAMOND_PATH="/path/to/your/folder/diamond"
HMMER_PATH="/path/to/your/folder/hmmer-3.4/bin"
MMSEQS_PATH="/path/to/your/folder/mmseqs"
PROST_PATH="/path/to/your/folder/prost.py"
HHMAKE_PATH="/path/to/your/folder/hhmake"
HHBLITS_PATH="/path/to/your/folder/hhblits"
BENCHMARK_FASTA="benchmark_dataset_SCOPe.fasta"
PREVIOUS_VTR_FASTA="previous_vtr_seq.fasta"
NCBI_UNIPROT_FASTA="human_virus_sequences"
WORKING_DIR="/path/to/your/folder/benchmark/"
FAMILY_INFO="family_info_benchmark.txt"  # File containing family information for annotation

# Define output paths for the results
DIAMOND_OUTPUT="benchmark_orthologs_blast.tsv"
JACKHMMER_OUTPUT="benchmark_orthologs_jackhmmer.txt"
MMSEQS_OUTPUT="benchmark_orthologs_mmseqs.txt"
PROST_OUTPUT="benchmark_orthologs_prost.tsv"
HHBLITS_OUTPUT="benchmark_orthologs_hhblits"

# Filtered output files
DIAMOND_FILTERED="1_diamond/benchmark_blastp_filtered.csv"
MMSEQS_FILTERED="2_mmseq/benchmark_mmseqs_filtered.csv"
HHBLITS_FILTERED="3_hhblits/benchmark_hhblits_filtered.csv"
JACKHMMER_FILTERED="4_jackhmmer/benchmark_jackhmmer_filtered.csv"
PROST_FILTERED="5_prost/benchmark_prost_filtered.tsv"

# Step 1: Create database with DIAMOND
$DIAMOND_PATH makedb --in $BENCHMARK_FASTA -d benchmark

# Step 2: Run DIAMOND blastp
$DIAMOND_PATH blastp -d benchmark -q $BENCHMARK_FASTA -o $DIAMOND_OUTPUT --evalue 1e-6 --query-cover 50 --id 50 --outfmt 6

# Step 3: Run Jackhmmer
$HMMER_PATH/jackhmmer --tblout $JACKHMMER_OUTPUT -o $JACKHMMER_OUTPUT $BENCHMARK_FASTA $BENCHMARK_FASTA

# Step 4: Run MMseqs2
$MMSEQS_PATH easy-rbh $BENCHMARK_FASTA $BENCHMARK_FASTA $MMSEQS_OUTPUT tmp

# Step 5: Run HHblits
$HHBLITS_PATH -i $BENCHMARK_FASTA -d benchmark -o $HHBLITS_OUTPUT

# Step 6: Run PROST
$PROST_PATH makedb $BENCHMARK_FASTA benchmark
$PROST_PATH search --thr 0.05 --jobs 60 $BENCHMARK_FASTA benchmark $PROST_OUTPUT

# Function to add family information to the results
add_family_info() {
    local input_file=$1
    local output_file=$2
    local family_info=$3

    # Add family information to results
    awk 'BEGIN {FS=OFS="\t"} NR==FNR {fam[$1]=$2; next} {print $0, fam[$1], fam[$2]}' $family_info $input_file > $output_file
}

# Step 7: Add family information to each tool's output
add_family_info $DIAMOND_OUTPUT $DIAMOND_FILTERED $FAMILY_INFO
add_family_info $JACKHMMER_OUTPUT $JACKHMMER_FILTERED $FAMILY_INFO
add_family_info $MMSEQS_OUTPUT $MMSEQS_FILTERED $FAMILY_INFO
add_family_info $PROST_OUTPUT $PROST_FILTERED $FAMILY_INFO
add_family_info $HHBLITS_OUTPUT $HHBLITS_FILTERED $FAMILY_INFO

# Step 8: Run R script to generate AUC values and ROC curves

Rscript -e "
library(pROC)
library(dplyr)
library(ggplot2)

# DIAMOND AUC and ROC
data <- read.csv('1_diamond/benchmark_blastp_filtered.csv', check.names = FALSE)
data\$True_Labels <- ifelse(data\$Query_Fam == data\$Target_Fam, 1, 0)
data <- data %>% filter(Query != Target)
roc(data\$True_Labels, data\$bitScore, plot = TRUE, legacy.axis = TRUE, col='orange', lwd=2, print.auc=TRUE, print.auc.x = 0.2, print.auc.y = 0.2)

auc_diamond <- roc(data$True_Labels, data$bitScore)$auc

# MMseqs2 AUC and ROC
data <- read.csv('2_mmseq/benchmark_mmseqs_filtered.csv', check.names = FALSE)
data\$True_Labels <- ifelse(data\$Query_Fam == data\$Target_Fam, 1, 0)
data <- data %>% filter(Query != Target & eValue < 1e-6)
roc(data\$True_Labels, data\$bitScore, plot = TRUE, legacy.axis = TRUE, col='darkgreen', lwd=2, print.auc=TRUE, add=TRUE, print.auc.x = 0.2, print.auc.y = 0.1)
auc_mmseqs <- roc(data$True_Labels, data$bitScore)$auc

# Jackhmmer AUC and ROC
data <- read.csv('4_jackhmmer/benchmark_jackhmmer_filtered.csv', check.names = FALSE)
data\$True_Labels <- ifelse(data\$Query_Fam == data\$Target_Fam, 1, 0)
data <- data %>% filter(Query != Target & E-value < 1e-6)
roc(data\$True_Labels, data\$score, plot = TRUE, legacy.axis = TRUE, col='blue', lwd=2, print.auc=TRUE, add=TRUE, print.auc.x = 0.2, print.auc.y = 0.4)

auc_jackhmmer <- roc(data$True_Labels, data$score)$auc

# PROST AUC and ROC
data <- read.csv('5_prost/benchmark_prost_filtered.tsv', check.names = FALSE, sep='\t')
data\$True_Labels <- ifelse(data\$Query_Fam == data\$Target_Fam, 1, 0)
data <- data %>% filter(Query != Target & evalue < 1e-6)
roc(data\$True_Labels, data\$distance, plot = TRUE, legacy.axis = TRUE, col='violet', lwd=2, print.auc=TRUE, add=TRUE, print.auc.x = 0.2, print.auc.y = 0.5)
auc_prost <- roc(data$True_Labels, data$distance)$auc

# HHblits AUC and ROC
data <- read.csv('3_hhblits/benchmark_hhblits_filtered.csv', check.names = FALSE, sep=',')
data\$True_Labels <- ifelse(data\$Query_Fam == data\$Target_Fam, 1, 0)
data <- data %>% filter(Query != Target & E-value < 1e-6)
roc(data\$True_Labels, data\$Probablility, plot = TRUE, legacy.axis = TRUE, col='brown', lwd=2, print.auc=TRUE, add=TRUE, print.auc.x = 0.2, print.auc.y = 0.3)

auc_hhblits <- roc(data$True_Labels, data$Probablility)$auc

# Aggregated score ROC
result\$True_Labels <- ifelse(result\$Query_Fam == result\$Target_Fam, 1, 0)
result <- result %>% filter(Query != Target)
roc(result\$True_Labels, result\$sum_aggregate_score, plot = TRUE, legacy.axis = TRUE, col='red', lwd=2, print.auc=TRUE, add=TRUE, print.auc.x = 0.2, print.auc.y = 0.6)

dev.off()

# Step 9: Calculate weights based on AUC values
weights <- c(
  diamond = 1 / (1 - auc_diamond),
  mmseqs = 1 / (1 - auc_mmseqs),
  jackhmmer = 1 / (1 - auc_jackhmmer),
  prost = 1 / (1 - auc_prost),
  hhblits = 1 / (1 - auc_hhblits)
)
weights <- weights / sum(weights)

# Assign weights to each tool
weight_diamond <- unname(weights["diamond"])
weight_mmseqs <- unname(weights["mmseqs"])
weight_jackhmmer <- unname(weights["jackhmmer"])
weight_prost <- unname(weights["prost"])
weight_hhblits <- unname(weights["hhblits"])

# Step 10: Calculate aggregated scores for each tool

# Diamond data
data_diamond <- read.csv('1_diamond/benchmark_blastp_filtered.csv', check.names = FALSE)
data_diamond\$True_Labels <- ifelse(data_diamond\$Query_Fam == data_diamond\$Target_Fam, 1, 0)
data_diamond\$min_max_score <- (data_diamond\$bitScore - min(data_diamond\$bitScore)) / (max(data_diamond\$bitScore) - min(data_diamond\$bitScore))
data_diamond\$aggregated_score <- data_diamond\$min_max_score * weight_diamond

# MMseqs data
data_mmseqs <- read.csv('2_mmseq/benchmark_mmseqs_filtered.csv', check.names = FALSE)
data_mmseqs\$True_Labels <- ifelse(data_mmseqs\$Query_Fam == data_mmseqs\$Target_Fam, 1, 0)
data_mmseqs\$min_max_score <- (data_mmseqs\$bitScore - min(data_mmseqs\$bitScore)) / (max(data_mmseqs\$bitScore) - min(data_mmseqs\$bitScore))
data_mmseqs\$aggregated_score <- data_mmseqs\$min_max_score * weight_mmseqs

# Jackhmmer data
data_jackhmmer <- read.csv('4_jackhmmer/benchmark_jackhmmer_filtered.csv', check.names = FALSE)
data_jackhmmer\$True_Labels <- ifelse(data_jackhmmer\$Query_Fam == data_jackhmmer\$Target_Fam, 1, 0)
data_jackhmmer\$min_max_score <- (data_jackhmmer\$score - min(data_jackhmmer\$score)) / (max(data_jackhmmer\$score) - min(data_jackhmmer\$score))
data_jackhmmer\$aggregated_score <- data_jackhmmer\$min_max_score * weight_jackhmmer

# Prost data
data_prost <- read.csv('5_prost/benchmark_prost_filtered.tsv', check.names = FALSE, sep='\t')
data_prost\$True_Labels <- ifelse(data_prost\$Query_Fam == data_prost\$Target_Fam, 1, 0)
data_prost\$min_max_score <- (data_prost\$distance - min(data_prost\$distance)) / (max(data_prost\$distance) - min(data_prost\$distance))
data_prost\$aggregated_score <- data_prost\$min_max_score * weight_prost

# HHblits data
data_hhblits <- read.csv('3_hhblits/benchmark_hhblits_filtered.csv', check.names = FALSE, sep=',')
data_hhblits\$True_Labels <- ifelse(data_hhblits\$Query_Fam == data_hhblits\$Target_Fam, 1, 0)
data_hhblits\$min_max_score <- (data_hhblits\$Probablility - min(data_hhblits\$Probablility)) / (max(data_hhblits\$Probablility) - min(data_hhblits\$Probablility))
data_hhblits\$aggregated_score <- data_hhblits\$min_max_score * weight_hhblits

#Step 11:Subsetting of data and combining
column_numbers <- c(1,2,3,4,5,6,16,17,18)
data_diamond1 <- data_diamond[, column_numbers]
data_diamond1$Tool <- "Diamond"

column_numbers_mmseq <- c(1,2,3,4,5,6,16,17,18)
data_mmseq1 <- data_mmseq[, column_numbers_mmseq]
data_mmseq1$Tool <- "mmseq"

column_numbers_hmmer <- c(1,2,3,4,5,6,13,14,15)
data_hmmer1 <- data_hmmer[, column_numbers_hmmer]
data_hmmer1$Tool <- "hmmer"

column_numbers_prost <- c(1,2,3,4,5,6,9,10,11)
data_prost1 <- data_prost[, column_numbers_prost]
data_prost1$Tool <- "prost"

column_numbers_hhblits <- c(1,2,3,4,5,6,11,12,13)
data_hhblits1 <- data_hhblits[, column_numbers_hhblits]
data_hhblits1$Tool <- "hhblits"
## binding
combined_df <- rbind(data_diamond1, data_mmseq1 , data_hmmer1, data_hhblits1, data_prost1)
cols <- c("Query","Query_Fam","Target", "Target_Fam")
result <- combined_df %>%
  group_by(across(all_of(cols))) %>%
  summarize(sum_aggregate_score = sum(aggregated_score))

#Step 12: Determing the cut-off using Youden J statistic
roc_curve <- roc(result$True_Labels, result$sum_aggregate_score)
coords(roc_curve, "best", "threshold")#, best.method = "closest.topleft")

#Step 13: Filtering based on optimal threshold
cutoff <- 0.076
result$significant <- result$sum_aggregate_score > cutoff
subset_combined_df <- subset(result, significant == TRUE)

# Step 12: Generate final output for comparison
write.csv(subset_combined_df, 'final_aggregated_scores.csv')
"
#Repeat Step 1-6,Step 10,11,13,14 for vTRs homologs identification
#Optional cd-hit -i vtr_homologs_refined.fasta -c 0.8 -o clustered_homologs.txt
