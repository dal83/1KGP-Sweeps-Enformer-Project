# Filename: normalize_preds.R
# Part of: 1KGP_sweeps_enformer_project
# Contains: Script to retrieve sums of Enformer predictions (across all 896 prediction windows) and normalize data across Enformer tracks.
# Debbie Lilly, July 2023

library(reticulate)
np <- import("numpy")

options(digits = 12)
options(max.print = 1000)

# Replace with location of all the directories containing the predictions (as npz files)
npz_dir_path = "enformer_scripts/LCT_middle_window_preds"
npz_dirs <- list.dirs(path = npz_dir_path, full.names = TRUE, recursive = FALSE)

# Every npz file should have the same name, just within a unique directory.
file_name <- "chr2:135653892-136047107.npz"

# Function sums all 896 predictions given for each Enformer track.
# Result is one prediction for each track.
# Function is applied to all npz files.
full_matrix <- sapply(npz_dirs, function(dir) {
		file <- np$load(file.path(dir, file_name))
		array <- file$f[["arr_0"]]
		colSums(array)
})

# Writing summed prediction values to a tabbed file.
# This isn't necessary if you just need normalized values.
write.table(full_matrix, "enformer_scripts/middle_window_sums2.tsv", sep="\t")


# Normalizing values:

# Creating empty matrix with same dimensions as full_matrix.
new_matrix <- matrix(NA, nrow = nrow(full_matrix), ncol = ncol(full_matrix))

# Normalizing values within each track.
for (r in 1:nrow(full_matrix)) {
		row <- full_matrix[r, ]
		row_mean <- mean(row)
		row_sd <- sd(row)
		new_matrix[r, ] <- (row - row_mean) / row_sd
}

# Writing normalized prediction values to a tabbed file.
file_out <- "enformer_scripts/norm_middle_window_sums2.tsv" # Replace with your file out
write.table(new_matrix, file_out, sep="\t")
