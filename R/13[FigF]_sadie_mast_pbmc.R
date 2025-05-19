rm(list = ls()); gc()

# library(seutools)
library(Seurat)
library(SeuratWrappers)
library(Matrix)
# dir.create(path = "/media/WD24/sadie10x/rna/mast", showWarnings = FALSE)
dir.create(path = "/media/WD24/sadie10x/rna/mast/pbmc_mpa_mono_sd", showWarnings = FALSE)
setwd("/media/WD24/sadie10x/rna/mast/pbmc_mpa_mono_sd")

seu <- readRDS(file = "/media/WD24/sadie10x/rna/object/MM148_pbmc_seurat.rds")
seu$cell_type <- ifelse(seu$cell_type=="M-platelet", "MPA", seu$cell_type)
seu$cell_type <- ifelse(seu$cell_type=="CD14 Mono", "cMono", seu$cell_type)
seu$cell_type <- ifelse(seu$cell_type=="CD16 Mono", "nMono", seu$cell_type)

# PANEL F
if(F) {
  test_clusters <- c("MPA","cMono","nMono")
  test_cats <- c("SD1","SD3")
  cat_col <- "study_day"
  
  seu_test <- subset(x = seu, subset = cell_type %in% test_clusters)
  seu_test@meta.data$barcode <- seu_test$barcode_2
  
  mast_res <- vector("list", length = length(test_clusters)); names(mast_res) <- test_clusters; mast_gsea_res <- mast_res
  
  for(i in 1:length(mast_res)) {
    start_mast_i <- Sys.time()
    print(paste0("starting: ",names(mast_res)[i]," [",i,"] of [",length(mast_res),"] clusters"))
    target_cluster <- names(mast_res)[i]
    seu_test$condition_custom <- "media"
    
    mast_seu <- subset(x = seu_test, subset = cell_type == target_cluster) # subset for cluster
    mast_seu <- subset(x = mast_seu, cells = which(mast_seu@meta.data[,cat_col] %in% test_cats)) # subset for categories
    
    try(seu_mast <- seutools::seurat_dge(seurat_object = mast_seu,
                                         dge_method = "mast",
                                         assay = "RNA",
                                         freq_expressed = 0.1,
                                         fc_threshold = log2(1.25),
                                         test_clusters = target_cluster,
                                         cluster_column = "cell_type",
                                         category_column = cat_col,
                                         test_categories = test_cats,
                                         test_condition = "all",
                                         condition_column = "condition_custom", # pseudo condition to prevent errors and customize output list name(s); cells were not perturbed, therefore going with 'media'
                                         pid_column = "study_id",
                                         pseudobulk_test_mode = "cluster_by_category",
                                         filter_genes = "outer"), silent = TRUE)
    mast_res[[i]] <- seu_mast
  }
  
  mast_res_df <- do.call(rbind, lapply(X = mast_res, FUN = function(arg1) return(arg1[[1]][[1]][["raw_res"]]))) # if full object does not save out.. Error: C stack usage  7972852 is too close to the limit
  
  saveRDS(object = mast_res_df, file = "/media/WD24/sadie10x/rna/mast/pbmc_mpa_mono_sd/pbmc_mpa_mono_mast_deg_sd1sd3_15JAN2025.rds")
}
