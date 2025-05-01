rm(list = ls()); gc()

library(seutools)
library(Seurat)
library(SeuratWrappers)
library(Matrix)
dir.create(path = "/media/WD24/sadie10x/rna/mast", showWarnings = FALSE)
outpath <- "/media/WD24/sadie10x/rna/mast/mpa_real_sim_drop_platelet_genes"
dir.create(path = outpath, showWarnings = FALSE)
setwd(outpath)

seu <- readRDS(file = "/media/WD24/sadie10x/rna/object/MM148_mpa_sim_seurat.rds")

mast_platelet_gene <- read.csv(file = "/media/WD24/sadie10x/rna/mast/platelet_genes/mast_platelet_gene_result.csv")
drop_genes <- mast_platelet_gene$primerid[intersect(mast_platelet_gene$fdr<0.05, mast_platelet_gene$logFC>1)]

seu_test_mpa <- subset(x = seu, subset = droplet_type %in% c("MPA_real", "MPA_sim"))
seu_test_mpa$condition_custom <- "media"
seu_test_mpa$cluster_custom <- "MPA"

start_mast <- Sys.time()
mpa_mast <- seutools::seurat_dge(seurat_object = seu_test_mpa,
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

saveRDS(object = mpa_mast[["media"]][["MPA"]][["raw_res"]], file = "/media/WD24/sadie10x/rna/mast/mpa_real_sim/pbmc_mpa_sim_real_no_platelet_genes_deg_6FEB2025.rds")
