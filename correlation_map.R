# Filename: correlation_map.R
# Part of: 1KGP_sweeps_enformer_project
# Contains: Script to create a heatmap based on the normalized predictions at
# the LCT trancscription start site and the LCT causal SNP in samples containing
# the reference allele.
# Debbie Lilly, July 2023

library(readr)
library("gplots")
library(tidyverse)

# Load files
SNP_mat <- as.matrix(readr::read_tsv("/vast/palmer/home.mccleary/dal83/project/1KGP_sweeps_enformer_project/scripts/enformer_scripts/LCT_middle_window/norm_LCT_SNP_window_sums.tsv"))
LCT_mat <- as.matrix(readr::read_tsv("/vast/palmer/home.mccleary/dal83/project/1KGP_sweeps_enformer_project/scripts/enformer_scripts/norm_transcription_window_sums_labeled.tsv"))

# Clean up tables
rownames(SNP_mat) <- SNP_mat[,1]
SNP_mat <- SNP_mat[,-1]
rownames(LCT_mat) <- LCT_mat[,1]
LCT_mat <- LCT_mat[,-1]

# Filter matrices to exclude non-DNASE tracks
filtered_SNP_data <- as.data.frame(SNP_mat) %>% filter(grepl("DNASE", rownames(SNP_mat)))
filtered_LCT_data <- as.data.frame(LCT_mat) %>% filter(grepl("CAGE", rownames(LCT_mat)))

# Filter matrices to exclude alt allele samples
gens <- t(as.data.frame(readr::read_tsv("/vast/palmer/home.mccleary/dal83/project/1KGP_sweeps_enformer_project/scripts/heatmaps/Book3.tsv")))
colnames(gens) <- "val"
filtered_gens <- filter(as.data.frame(gens), val == "0")
common_row_names <- c(rownames(filtered_gens))
reverse_SNP <- t(filtered_SNP_data)
combined_SNP <- reverse_SNP[rownames(reverse_SNP) %in% common_row_names, ]
reverse_LCT <- t(filtered_LCT_data)
combined_LCT <- reverse_LCT[rownames(reverse_LCT) %in% common_row_names, ]

# Build correlation matrix
cor_mat <- matrix(nrow = ncol(combined_SNP), ncol = ncol(combined_LCT))
rownames(cor_mat) <- c(colnames(combined_SNP))
colnames(cor_mat) <- c(colnames(combined_LCT))
for(i in 1:ncol(combined_SNP)) {
  for(j in 1:ncol(combined_LCT)) {
    cor_mat[i, j] <- cor(as.double(combined_SNP[,i]), as.double(combined_LCT[,j]))
  }
}

# Create heatmap. Red values are negative, blue values are positive.
cols <- RColorBrewer:::brewer.pal(11,"RdBu")
heatmap(cor_mat[1:20, 1:20], col = cols)
