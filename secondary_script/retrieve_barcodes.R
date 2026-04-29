args <- commandArgs(trailingOnly = TRUE)

lane_dir <- if (length(args) >= 1) args[1] else "."
file_input <- file.path(lane_dir, "filtered_feature_bc_matrix.h5")
file_output <- file.path(lane_dir, "barcodes_R.tsv")

h5_in <- Seurat::Read10X_h5(file_input)
bc <- h5_in@Dimnames[[2]]

write.table(
	x = data.frame(val = bc),
	file = file_output,
	quote = FALSE,
	sep = "\n",
	col.names = FALSE,
	row.names = FALSE
)
