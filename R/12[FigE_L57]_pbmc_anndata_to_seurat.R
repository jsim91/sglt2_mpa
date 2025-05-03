# /media/MPEdge16/MM137/sc/py/py_out/pbmc/anndata_elements/pbmc_anndata_to_seurat.R
rm(list = ls()); gc()

library(Seurat)
library(seutools)
library(Matrix)
source(file.path(system.file(package = "seutools"),"source_scripts/mast_src.R"))

setwd("/media/MPEdge16/MM137/sc/py/py_out/pbmc/anndata_elements")

mtx <- Matrix::readMM(file = "adata_pbmc_counts.mtx")
meta <- read.csv(file = "adata_pbmc_obs.csv", row.names = 1); colnames(meta)[which(colnames(meta)=="barcode_2.1")] <- "barcode_lane"; head(meta,3)
meta$barcode <- meta$barcode_lane
genes <- read.csv(file = "adata_pbmc_var.csv")[,1]; head(genes)
umap_coord <- read.csv(file = "adata_pbmc_umap_coordinates.csv", header = FALSE); colnames(umap_coord) <- c("UMAP1","UMAP2"); row.names(umap_coord) <- meta$barcode_lane; head(umap_coord)
latent <- read.csv(file = "adata_pbmc_latent_coordinates.csv", header = FALSE); colnames(latent) <- paste0("latent",1:ncol(latent)); row.names(latent) <- meta$barcode_lane; head(latent)

seu <- Seurat::CreateSeuratObject(counts = Matrix::t(mtx), assay = "RNA", meta.data = meta)
colnames(seu) <- seu@meta.data$barcode_lane
row.names(seu) <- genes
umap_reduc <- SeuratObject::CreateDimReducObject(embeddings = as.matrix(umap_coord), assay = "RNA", key = "umap_")
seu[['umap']] <- umap_reduc
latent_reduc <- SeuratObject::CreateDimReducObject(embeddings = as.matrix(latent), assay = "RNA", key = "latent_")
seu[['scvi_latent']] <- latent_reduc

DefaultAssay(seu)
seu <- NormalizeData(object = seu, assay = "RNA", normalization.method = "LogNormalize")
Idents(seu) <- seu@meta.data$subset_cluster

# saveRDS(object = seu, file = "/media/MPEdge16/MM137/sc/py/py_out/pbmc/anndata_elements/pbmc_seurat.rds")
seu <- readRDS(file = "/media/MPEdge16/MM137/sc/py/py_out/pbmc/anndata_elements/pbmc_seurat.rds")

seu <- AddMetaData(object = seu, metadata = paste0(seu$study_id,"_",seu$study_day), col.name = "full_id")
uid <- unique(seu$full_id); uid <- uid[order(uid)]
uclus <- unique(seu$subset_cluster); uclus <- uclus[order(uclus)]

fmat <- matrix(data = NA, nrow = length(uid), ncol = length(uclus))
row.names(fmat) <- uid; colnames(fmat) <- uclus
for(i in 1:nrow(fmat)) {
  clnums <- seu$subset_cluster[seu$full_id==row.names(fmat)[i]]
  for(j in 1:ncol(fmat)) {
    # fmat[i,j] <- mean(clnums==colnames(fmat)[j]) * 100
    fmat[i,j] <- sum(clnums==colnames(fmat)[j])
  }
}
write.csv(x = fmat, file = "/media/MPEdge16/MM137/sc/py/py_out/pbmc/anndata_elements/mm_pbmc_subset_cluster_cell_counts.csv")

uct <- unique(seu$cell_type); uct <- uct[order(uct)]
fmat2 <- matrix(data = NA, nrow = length(uid), ncol = length(uct))
row.names(fmat2) <- uid; colnames(fmat2) <- uct
for(i in 1:nrow(fmat2)) {
  clnames <- seu$cell_type[seu$full_id==row.names(fmat2)[i]]
  for(j in 1:ncol(fmat2)) {
    fmat2[i,j] <- sum(clnames==colnames(fmat2)[j])
  }
}
# PANEL E
write.csv(x = fmat2, file = "/media/MPEdge16/MM137/sc/py/py_out/pbmc/anndata_elements/mm_pbmc_named_cluster_cell_counts.csv")
