#args <- commandArgs(trailingOnly = TRUE)

#setwd(paste0("/media/MPEdge16/MM137/sc/10x_cloud_dl/10872-MM-",args[1],"_standard"))
file_input <- "filtered_feature_bc_matrix.h5"
h5_in <- Seurat::Read10X_h5(file_input)
bc <- h5_in@Dimnames[[2]]
write.table(x = data.frame(val = bc), file = "barcodes_R.tsv", quote = FALSE, sep = "\n", col.names = FALSE, row.names = FALSE)
