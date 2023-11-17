# Filename: heatmap.R
# Part of: 1KGP_sweeps_enformer_project
# Contains: Script to create a heatmap basd on the normalized Enformer predictions.
# Debbie Lilly, July 2023

library(readr)
library("gplots")

tsv_file <- readr::read_tsv("/Users/deborahlilly/Documents/R_prac/norm_middle_window_sums_labeled2.tsv")
matrix_data <- as.matrix(tsv_file)
rownames(matrix_data) <- matrix_data[,1]
num_matrix <- matrix(as.numeric(matrix_data), nrow=nrow(matrix_data), ncol=ncol(matrix_data))
rownames(num_matrix) <- rownames(matrix_data)
colnames(num_matrix) <- colnames(matrix_data)

# Creates heatmap. X-axis: 1KGP samples, y-axis: Enformer tracks.
heatmap(num_matrix)
