rm(list = ls()); gc()

library(here)
source(here::here("primary_dependents/seurat_dge_local.R"))
library(Seurat)
library(SeuratWrappers)
library(Matrix)
dir.create(here::here("intermediate/mast/mpa_real_sim"), showWarnings = FALSE, recursive = TRUE)

ct <- Matrix::t(Matrix::readMM(file = here::here("intermediate/pbmc/anndata_elements/adata_pbmc_counts_mpa_sim.mtx")))
obs <- read.csv(file = here::here("intermediate/pbmc/anndata_elements/adata_pbmc_obs_mpa_sim.csv"))
var <- read.csv(file = here::here("intermediate/pbmc/anndata_elements/adata_pbmc_var_mpa_sim.csv"))

colnames(ct) <- obs$barcode_2
row.names(ct) <- var[,1]

seu <- Seurat::CreateSeuratObject(counts = ct, assay = "RNA", meta.data = obs)
head(seu@assays$RNA@layers$counts@x,20) # should be ints
seu <- Seurat::NormalizeData(object = seu, assay = "RNA", normalization.method = "LogNormalize")

Idents(seu) <- seu$cell_type

mast_platelet_gene <- read.csv(file = here::here("intermediate/mast/platelet_genes/mast_platelet_gene_result.csv"))
drop_genes <- mast_platelet_gene$primerid[mast_platelet_gene$fdr < 0.05 & mast_platelet_gene$logFC > 1]

seu_test_mpa <- subset(x = seu, subset = droplet_type %in% c("MPA_real", "MPA_sim"))
seu_test_mpa$condition_custom <- "media"
seu_test_mpa$cluster_custom <- "MPA"

start_mast <- Sys.time()
mpa_mast <- seurat_dge(seurat_object = seu_test_mpa,
                                 dge_method = "mast",
                                 assay = "RNA",
                                 freq_expressed = 0.1,
                                 fc_threshold = log2(1.25),
                                 test_clusters = "MPA",
                                 cluster_column = "cluster_custom",
                                 category_column = "droplet_type",
                                 test_categories = c("MPA_sim","MPA_real"),
                                 test_condition = "media",
                                 condition_column = "condition_custom", # pseudo condition to prevent errors and customize output list name(s); cells were not perturbed, therefore going with 'media'
                                 pid_column = "study_id",
                                 pseudobulk_test_mode = "cluster_by_category",
                                 filter_genes = "outer", 
                                 gene_set_blacklist = drop_genes)
difftime(Sys.time(), start_mast, units = "mins")

saveRDS(object = mpa_mast[["media"]][["MPA"]][["raw_res"]], file = here::here("intermediate/mast/mpa_real_sim/pbmc_mpa_sim_real_no_platelet_genes_deg_6FEB2025.rds"))
