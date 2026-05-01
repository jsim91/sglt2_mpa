# MM148 scds doublet finding with scds
# Produces per-lane {lane}_scds_results.tsv files in intermediate/pbmc/.
# Must run before 6c_build_doublet_scores.ipynb (which merges all algorithm results).
rm(list = ls()); gc()

library(Matrix)
library(SingleCellExperiment)
library(scds)
library(here)

outdir <- here::here("intermediate/pbmc")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

lanes <- paste0("MM-", 1:8)

for (i in 1:length(lanes)) {
  print(paste0('starting scds on: ', lanes[i]))
  ct <- Matrix::readMM(file = here::here(paste0("primary_dependents/counts_from_10x_h5/10872-", lanes[i], "_counts.mtx")))
  obs <- read.csv(file = here::here(paste0("primary_dependents/counts_from_10x_h5/10872-", lanes[i], "_obs_names.csv")))
  var <- read.csv(file = here::here(paste0("primary_dependents/counts_from_10x_h5/10872-", lanes[i], "_var_names.csv")))
  tmp_sce <- SingleCellExperiment::SingleCellExperiment(list(counts = ct))
  rownames(tmp_sce) <- var[, 1]
  colnames(tmp_sce) <- obs[, 1]

  soc <- read.csv(file = here::here(paste0("primary_dependents/souporcell/10872-", lanes[i], "_soc/clusters.tsv")), sep = '\t')
  singlet_bc <- soc$barcode[soc$status == "singlet"]
  singlet_bc <- gsub(pattern = "\\-[0-9]+$", replacement = gsub(pattern = "MM", replacement = "", x = lanes[i]), x = singlet_bc)
  tmp_sce <- tmp_sce[, colnames(tmp_sce) %in% singlet_bc]

  exprs_matrix <- assay(tmp_sce)
  cell_counts <- rowSums(exprs_matrix > 0)
  cell_proportions <- cell_counts / ncol(exprs_matrix)
  threshold <- 0.01
  genes_to_keep <- cell_proportions >= threshold
  tmp_sce <- tmp_sce[genes_to_keep, ]

  tmp_sce <- cxds(tmp_sce, retRes = TRUE, estNdbl = TRUE) # just using cxds score; bcds does not seem to work well on this dataset
  Doublets <- as.data.frame(cbind(rownames(colData(tmp_sce)),
                                  colData(tmp_sce)$cxds_score,
                                  colData(tmp_sce)$cxds_call))
  colnames(Doublets) <- c("Barcode", "cxds_score", "cxds_DropletType")
  Doublets$cxds_DropletType <- gsub("FALSE", "singlet", Doublets$cxds_DropletType)
  Doublets$cxds_DropletType <- gsub("TRUE", "doublet", Doublets$cxds_DropletType)
  print(paste0(sum(Doublets$cxds_DropletType == "doublet"), " doublets found for ", lanes[i], " of ", nrow(Doublets), " cells ..."))
  write.table(x = Doublets, file = file.path(outdir, paste0(lanes[i], "_scds_results.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

for (i in 1:length(lanes)) {
  res <- read.csv(file = file.path(outdir, paste0(lanes[i], "_scds_results.tsv")), sep = '\t')
  write.table(x = res, file = file.path(outdir, paste0(lanes[i], "_scds_results.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
}
