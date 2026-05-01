rm(list = ls()); gc()

library(Seurat)
library(Matrix)
library(here)

# Read myeloid-refined annotations and matrices from script 5
mtx <- Matrix::readMM(file = here::here("intermediate/pbmc/anndata_elements/adata_pbmc_myeloid_counts.mtx"))
meta <- read.csv(file = here::here("intermediate/pbmc/anndata_elements/adata_pbmc_myeloid_obs.csv"), row.names = 1); colnames(meta)[which(colnames(meta)=="barcode_2.1")] <- "barcode_lane"; head(meta,3)
meta$barcode <- meta$barcode_lane
genes <- read.csv(file = here::here("intermediate/pbmc/anndata_elements/adata_pbmc_myeloid_var.csv"))[,1]; head(genes)
umap_coord <- read.csv(file = here::here("intermediate/pbmc/anndata_elements/adata_pbmc_myeloid_umap_coordinates.csv"), header = FALSE); colnames(umap_coord) <- c("UMAP1","UMAP2"); row.names(umap_coord) <- meta$barcode_lane; head(umap_coord)
latent <- read.csv(file = here::here("intermediate/pbmc/anndata_elements/adata_pbmc_myeloid_latent_coordinates.csv"), header = FALSE); colnames(latent) <- paste0("latent",1:ncol(latent)); row.names(latent) <- meta$barcode_lane; head(latent)

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

saveRDS(object = seu, file = here::here("intermediate/pbmc/anndata_elements/pbmc_seurat.rds"))

seu <- AddMetaData(object = seu, metadata = paste0(seu$study_id,"_",seu$study_day), col.name = "full_id")
uid <- unique(seu$full_id); uid <- uid[order(uid)]

uct <- unique(seu$cell_type); uct <- uct[order(uct)]
fmat <- matrix(data = NA, nrow = length(uid), ncol = length(uct))
row.names(fmat) <- uid; colnames(fmat) <- uct
for(i in 1:nrow(fmat)) {
  clnames <- seu$cell_type[seu$full_id==row.names(fmat)[i]]
  for(j in 1:ncol(fmat)) {
    fmat[i,j] <- sum(clnames==colnames(fmat)[j])
  }
}
write.csv(x = fmat, file = here::here("intermediate/pbmc/anndata_elements/mm_myeloid_named_cluster_cell_counts.csv"))
