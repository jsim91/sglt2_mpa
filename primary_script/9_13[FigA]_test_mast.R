rm(list = ls()); gc()

library(here)
source(here::here("primary_dependents/seurat_dge_local.R"))
library(Seurat)
library(SeuratWrappers)
library(Matrix)
library(ggplot2)
dir.create(here::here("intermediate/mast/pbmc_mpa_mono_sd"), showWarnings = FALSE, recursive = TRUE)
dir.create(here::here("output/figures"), showWarnings = FALSE, recursive = TRUE)

if(F) { # create seurat object from anndata elements once
  ct <- Matrix::t(Matrix::readMM(file = here::here("intermediate/pbmc/anndata_elements/adata_pbmc_counts.mtx")))
  obs <- read.csv(file = here::here("intermediate/pbmc/anndata_elements/adata_pbmc_obs.csv"))
  var <- read.csv(file = here::here("intermediate/pbmc/anndata_elements/adata_pbmc_var.csv"))
  umap_coord <- read.csv(file = here::here("intermediate/pbmc/anndata_elements/adata_pbmc_umap_coordinates.csv"), header = FALSE)
  colnames(umap_coord) <- c("UMAP1","UMAP2"); row.names(umap_coord) <- obs$barcode_2
  latent_coord <- read.csv(file = here::here("intermediate/pbmc/anndata_elements/adata_pbmc_latent_coordinates.csv"), header = FALSE)
  colnames(latent_coord) <- paste0("latent",1:ncol(latent_coord)); row.names(latent_coord) <- obs$barcode_2
  
  colnames(ct) <- obs$barcode_2
  row.names(ct) <- var[,1]
  
  seu <- Seurat::CreateSeuratObject(counts = ct, assay = "RNA", meta.data = obs)
  head(seu@assays$RNA@layers$counts@x,20) # should be ints
  seu <- Seurat::NormalizeData(object = seu, assay = "RNA", normalization.method = "LogNormalize")
  
  latent_dr <- Seurat::CreateDimReducObject(embeddings = as.matrix(latent_coord), assay = "RNA", key = "latent_")
  seu[['latent']] <- latent_dr
  
  umap_dr <- Seurat::CreateDimReducObject(embeddings = as.matrix(umap_coord), assay = "RNA", key = "umap_")
  seu[['umap']] <- umap_dr
  
  Seurat::Idents(seu) <- seu$cell_type
  saveRDS(object = seu, file = here::here("intermediate/seurat/MM148_pbmc_seurat.rds"))
} else { # then read in the seurat object on re-run
  seu <- readRDS(file = here::here("intermediate/seurat/MM148_pbmc_seurat.rds"))
  # format names
  seu$cell_type <- ifelse(seu$cell_type=="M-platelet", "MPA", seu$cell_type)
  seu$cell_type <- ifelse(seu$cell_type=="CD14 Mono", "cMono", seu$cell_type)
  seu$cell_type <- ifelse(seu$cell_type=="CD16 Mono", "nMono", seu$cell_type)
}

if(T) { # myeloid singlet umap
  myeloid_obs <- read.csv(file = here::here("intermediate/pbmc/anndata_elements/adata_pbmc_myeloid_obs.csv"), check.names = FALSE)
  myeloid_umap <- read.csv(file = here::here("intermediate/pbmc/anndata_elements/adata_pbmc_myeloid_umap_coordinates.csv"), check.names = FALSE, header = FALSE)
  colnames(myeloid_umap) <- c("UMAP1","UMAP2")
  
  myeloid_obs$merged_type <- ifelse(myeloid_obs$merged_type=="M-platelet", "MPA", myeloid_obs$merged_type)
  myeloid_obs$merged_type <- ifelse(myeloid_obs$merged_type=="CD14 Mono", "cMono", myeloid_obs$merged_type)
  myeloid_obs$merged_type <- ifelse(myeloid_obs$merged_type=="CD16 Mono", "nMono", myeloid_obs$merged_type)
  uclus <- unique(myeloid_obs$merged_type)
  
  myeloid_umap$cluster <- factor(myeloid_obs$merged_type)
  set.seed(123); myeloid_umap <- myeloid_umap[sample(1:nrow(myeloid_umap),nrow(myeloid_umap),replace=F),]
  
  clusx <- rep(NA, length=length(uclus)); names(clusx) <- uclus; clusy <- clusx
  
  for(i in 1:length(uclus)) {
    clusx[i] <- median(myeloid_umap$UMAP1[myeloid_umap$cluster==names(clusx)[i]])
    clusy[i] <- median(myeloid_umap$UMAP2[myeloid_umap$cluster==names(clusy)[i]])
  }
  anno_df <- data.frame(xval = clusx, yval = clusy, lab = names(clusx)); anno_df$lab <- factor(anno_df$lab)
  
  custom_col <- c('cMono' = '#FFB266', 
                  'nMono' = '#FF6666', 
                  'MPA' = '#66FF66', 
                  'cDC2' = '#66B2FF', 
                  'pDC' = 'pink', 
                  'cDC1' = '#B266FF', 
                  'Platelet' = '#66FFFF', 
                  'doublet' = '#A9A9A9')
  
  myl_umap <- ggplot() + 
    ggrastr::geom_point_rast(data = myeloid_umap, aes(x = UMAP1, y = UMAP2, color = cluster), alpha = 0.7) +
    annotate("text", x = anno_df$xval, y = anno_df$yval, label = anno_df$lab,
             hjust = 0.5, color = "black", fontface = "bold", size = 9) + 
    scale_color_manual(values = custom_col) +
    theme_bw() + 
    theme(legend.position = "none", 
          axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          axis.title = element_text(size = 24))
  ggsave(filename = "pbmc_myeloid_singlet_umap_16APR2025.pdf", plot = myl_umap, device = "pdf", path = here::here("output/figures"), 
         width = 8, height = 8, units = "in", dpi = 300, limitsize = F, bg = "white")
}

if(T) { # myeloid doublet umap
  myeloid_obs <- read.csv(file = here::here("intermediate/pbmc_myeloid_platelet_int_dbl_obs_to_r.csv"), check.names = FALSE)
  myeloid_umap <- read.csv(file = here::here("intermediate/pbmc_myeloid_platelet_int_dbl_umap_coordinates_to_r.csv"), check.names = FALSE, header = FALSE)
  colnames(myeloid_umap) <- c("UMAP1","UMAP2")
  
  myeloid_obs$ann_types <- ifelse(myeloid_obs$ann_types=="cMo", "cMono", myeloid_obs$ann_types)
  myeloid_obs$ann_types <- ifelse(myeloid_obs$ann_types=="nMo", "nMono", myeloid_obs$ann_types)
  uclus <- unique(myeloid_obs$ann_types)
  
  myeloid_umap$cluster <- factor(myeloid_obs$ann_types)
  set.seed(123); myeloid_umap <- myeloid_umap[sample(1:nrow(myeloid_umap),nrow(myeloid_umap),replace=F),]
  
  clusx <- rep(NA, length=length(uclus)); names(clusx) <- uclus; clusy <- clusx
  
  for(i in 1:length(uclus)) {
    clusx[i] <- median(myeloid_umap$UMAP1[myeloid_umap$cluster==names(clusx)[i]])
    clusy[i] <- median(myeloid_umap$UMAP2[myeloid_umap$cluster==names(clusy)[i]])
  }
  anno_df <- data.frame(xval = clusx, yval = clusy, lab = names(clusx)); anno_df$lab <- factor(anno_df$lab)
  anno_df$yval[anno_df$lab=="MPA"] <- 4.25
  
  custom_col <- c('cMono' = '#FFB266', 
                  'nMono' = '#FF6666', 
                  'MPA' = '#66FF66', 
                  'cDC2' = '#66B2FF', 
                  'pDC' = 'pink', 
                  'cDC1' = '#B266FF', 
                  'Platelet' = '#66FFFF', 
                  'doublet' = '#A9A9A9')
  
  myl_umap <- ggplot() + 
    ggrastr::geom_point_rast(data = myeloid_umap, aes(x = UMAP1, y = UMAP2, color = cluster), alpha = 0.7) +
    annotate("text", x = anno_df$xval, y = anno_df$yval, label = anno_df$lab,
             hjust = 0.5, color = "black", fontface = "bold", size = 9) + 
    scale_color_manual(values = custom_col) +
    theme_bw() + 
    theme(legend.position = "none", 
          axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          axis.title = element_text(size = 24))
  ggsave(filename = "pbmc_myeloid_doublet_umap_16APR2025.pdf", plot = myl_umap, device = "pdf", path = here::here("output/figures"), 
         width = 8, height = 8, units = "in", dpi = 300, limitsize = F, bg = "white")
}

if(T) { # myeloid sim, mpa
  myeloid_obs <- read.csv(file = here::here("intermediate/pbmc/anndata_elements/adata_pbmc_obs_mpa_sim.csv"), check.names = FALSE)
  myeloid_umap <- read.csv(file = here::here("intermediate/pbmc/anndata_elements/pbmc_mpa_sim_umap_coordinates.csv"), 
                           check.names = FALSE, header = FALSE, row.names = NULL)
  colnames(myeloid_umap) <- c("UMAP1","UMAP2")
  myeloid_umap$barcode <- myeloid_obs$barcode_2
  
  myeloid_mapping <- read.csv(file = '/media/MPEdge16/MM137/sc/py/py_out/pbmc/anndata_elements/adata_pbmc_obs_mpa_sim_2.csv')
  myeloid_mapping$sim_cell_type <- ifelse(myeloid_mapping$sim_cell_type=="M-platelet", "MPA", myeloid_mapping$sim_cell_type)
  myeloid_mapping$sim_cell_type <- ifelse(myeloid_mapping$sim_cell_type=="CD14 Mono", "cMono", myeloid_mapping$sim_cell_type)
  myeloid_mapping$sim_cell_type <- ifelse(myeloid_mapping$sim_cell_type=="CD16 Mono", "nMono", myeloid_mapping$sim_cell_type)
  myeloid_map <- myeloid_mapping$sim_cell_type; names(myeloid_map) <- myeloid_mapping$barcode_2
  
  myeloid_umap$cluster <- myeloid_map[myeloid_umap$barcode]
  myeloid_umap$cluster[grep(pattern = "\\-dbl$", x = myeloid_umap$barcode)] <- "SimMPA"
  myeloid_umap$cluster <- factor(myeloid_umap$cluster)
  
  myl_front <- myeloid_umap
  
  uclus <- unique(myl_front$cluster)

  set.seed(123); myl_front <- myl_front[sample(1:nrow(myl_front),nrow(myl_front),replace=F),]

  clusx <- rep(NA, length=length(uclus)); names(clusx) <- uclus; clusy <- clusx

  for(i in 1:length(uclus)) {
    clusx[i] <- median(myl_front$UMAP1[myl_front$cluster==names(clusx)[i]])
    clusy[i] <- median(myl_front$UMAP2[myl_front$cluster==names(clusy)[i]])
  }
  anno_df <- data.frame(xval = clusx, yval = clusy, lab = names(clusx)); anno_df$lab <- factor(anno_df$lab)
  
  custom_col <- c('cMono' = '#FFB266', 
                  'nMono' = '#FF6666', 
                  'MPA' = '#66FF66', 
                  'cDC2' = '#66B2FF', 
                  'pDC' = 'pink', 
                  'cDC1' = '#B266FF', 
                  'Platelet' = '#66FFFF', 
                  'doublet' = '#A9A9A9',
                  'SimMPA' = '#d50000')
  
  myl_umap <- ggplot() + 
    ggrastr::geom_point_rast(data = myl_front, aes(x = UMAP1, y = UMAP2, color = cluster, fill = cluster), alpha = 0.4) +
    annotate("text", x = anno_df$xval, y = anno_df$yval, label = anno_df$lab,
             hjust = 0.5, color = "black", fontface = "bold", size = 9) +
    scale_color_manual(values = custom_col) + 
    scale_fill_manual(values = custom_col) + 
    theme_bw() + 
    theme(legend.text = element_text(size = 18, face = 'bold'), 
          # legend.position = "bottom", 
          legend.position = 'none', 
          legend.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          axis.title = element_text(size = 24))
  ggsave(filename = "pbmc_myeloid_sim_umap_16APR2025.pdf", plot = myl_umap, device = "pdf", path = here::here("output/figures"), 
         width = 8, height = 8, units = "in", dpi = 300, limitsize = F, bg = "white")
}

if(T) {
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
    
    try(seu_mast <- seurat_dge(seurat_object = mast_seu,
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
    
    mast_gsea_res[[i]] <- seu_mast_gsea(mast_dge_result = seu_mast[['media']][[target_cluster]],
                                                   seu_mast_sets = c("C2_CP_Reactome","C5_GO_BP",
                                                                     "C5_GO_MF","H_Hallmark",
                                                                     "C2_CP_KEGG_LEGACY","C2_CP_KEGG_MEDICUS", 
                                                                     "C7_ImmuneSigDB","C7_VAX"),
                                                   num_boots = 50, gs_min = 3, prepare_plot_n = 20,
                                                   gs_max = Inf, gs_regex = NULL, nthread = 2)
    
    print(paste0("finished: ",names(mast_res)[i]," in ",round(as.numeric(difftime(Sys.time(), start_mast_i, units = "mins")),1)," mins"))
  }
  
  mast_res_df <- do.call(rbind, lapply(X = mast_res, FUN = function(arg1) return(arg1[[1]][[1]][["raw_res"]]))) # if full object does not save out.. Error: C stack usage  7972852 is too close to the limit
  
  saveRDS(object = mast_res_df, file = here::here("intermediate/mast/pbmc_mpa_mono_sd/pbmc_mpa_mono_mast_deg_sd1sd3_15JAN2025.rds"))
  
  saveRDS(object = mast_gsea_res, file = here::here("intermediate/mast/pbmc_mpa_mono_sd/pbmc_mpa_mono_mast_gsea_sd1sd3_15JAN2025.rds"))
}

if(T) { # find genes associated with platelets to drop when doing MPA vs CD14 Mono
  dir.create(here::here("intermediate/mast/platelet_genes"), showWarnings = FALSE, recursive = TRUE)
  
  seu_test_plt <- subset(x = seu, subset = cell_type %in% c("CD14 Mono", "Platelet"))
  seu_test_plt$condition_custom <- "media"
  seu_test_plt$cluster_custom <- "cMo_Platelet"
  
  start_mast <- Sys.time()
  seu_mast_plt <- seurat_dge(seurat_object = seu_test_plt,
                                       dge_method = "mast",
                                       assay = "RNA",
                                       freq_expressed = 0.1,
                                       fc_threshold = log2(1.25),
                                       test_clusters = "cMo_Platelet",
                                       cluster_column = "cluster_custom",
                                       category_column = "cell_type",
                                       test_categories = c("CD14 Mono","Platelet"),
                                       test_condition = "all",
                                       condition_column = "condition_custom", # pseudo condition to prevent errors and customize output list name(s); cells were not perturbed, therefore going with 'media'
                                       pid_column = "study_id",
                                       pseudobulk_test_mode = "cluster_by_category",
                                       filter_genes = "outer")
  difftime(Sys.time(), start_mast, units = "mins")
  
  save(seu_mast_plt, file = here::here("intermediate/mast/platelet_genes/cMo_vs_Platelet_mast.RData"))
  write.csv(x = seu_mast_plt[["media"]][["cMo_Platelet"]][["filtered_res"]][order(seu_mast_plt[["media"]][["cMo_Platelet"]][["filtered_res"]]$logFC,decreasing=T),], 
            file = here::here("intermediate/mast/platelet_genes/mast_platelet_gene_result.csv"), row.names = FALSE)
}

if(T) { # SD1 and SD3 (separately) MPA vs cMo
  platelet_mast_genes <- read.csv(file = here::here("intermediate/mast/platelet_genes/mast_platelet_gene_result.csv"), check.names = FALSE)
  blacklist_platelet_genes <- platelet_mast_genes$primerid[intersect(which(platelet_mast_genes$fdr<=0.05), which(platelet_mast_genes$logFC>0))]
  
  dir.create(here::here("intermediate/mast/mpa_sd_compare"), showWarnings = FALSE, recursive = TRUE)
  
  seu_test_m <- subset(x = seu, subset = cell_type %in% c("cMono", "MPA"))
  seu_test_sd1 <- subset(x = seu_test_m, subset = study_day == "SD1")
  seu_test_sd1$condition_custom <- "media"
  seu_test_sd1$cluster_custom <- "cMo_MPA"
  
  start_mast <- Sys.time()
  seu_mast_sd1 <- seurat_dge(seurat_object = seu_test_sd1,
                                       dge_method = "mast",
                                       assay = "RNA",
                                       freq_expressed = 0.1,
                                       fc_threshold = log2(1.25),
                                       test_clusters = "cMo_MPA",
                                       cluster_column = "cluster_custom",
                                       category_column = "cell_type",
                                       test_categories = c("cMono","MPA"),
                                       test_condition = "all",
                                       condition_column = "condition_custom", # pseudo condition to prevent errors and customize output list name(s); cells were not perturbed, therefore going with 'media'
                                       pid_column = "study_id",
                                       pseudobulk_test_mode = "cluster_by_category",
                                       filter_genes = "outer", 
                                       gene_set_blacklist = blacklist_platelet_genes)
  difftime(Sys.time(), start_mast, units = "mins")
  
  write.csv(x = seu_mast_sd1$media$cMo_MPA$raw_res, file = here::here("intermediate/mast/mpa_sd_compare/mast_dge_mpa_vs_cmono_sd1_29Jan2024.csv"), row.names = F)

  seu_test_sd3 <- subset(x = seu_test_m, subset = study_day == "SD3")
  seu_test_sd3$condition_custom <- "media"
  seu_test_sd3$cluster_custom <- "cMo_MPA"
  
  start_mast <- Sys.time()
  seu_mast_sd3 <- seurat_dge(seurat_object = seu_test_sd3,
                                       dge_method = "mast",
                                       assay = "RNA",
                                       freq_expressed = 0.1,
                                       fc_threshold = log2(1.25),
                                       test_clusters = "cMo_MPA",
                                       cluster_column = "cluster_custom",
                                       category_column = "cell_type",
                                       test_categories = c("cMono","MPA"),
                                       test_condition = "all",
                                       condition_column = "condition_custom", # pseudo condition to prevent errors and customize output list name(s); cells were not perturbed, therefore going with 'media'
                                       pid_column = "study_id",
                                       pseudobulk_test_mode = "cluster_by_category",
                                       filter_genes = "outer", 
                                       gene_set_blacklist = blacklist_platelet_genes)
  difftime(Sys.time(), start_mast, units = "mins")
  
  write.csv(x = seu_mast_sd3$media$cMo_MPA$raw_res, file = here::here("intermediate/mast/mpa_sd_compare/mast_dge_mpa_vs_cmono_sd3_29Jan2024.csv"), row.names = F)
  
  start_mast_gsea <- Sys.time()
  mast_gsea_sd1 <- seu_mast_gsea(mast_dge_result = seu_mast_sd1[["media"]][["cMo_MPA"]],
                                            seu_mast_sets = c("C2_CP_Reactome","C5_GO_BP",
                                                              "C5_GO_MF","H_Hallmark",
                                                              "C2_CP_KEGG_LEGACY","C2_CP_KEGG_MEDICUS", 
                                                              "C7_ImmuneSigDB","C7_VAX"),
                                            num_boots = 50, gs_min = 3, prepare_plot_n = 20,
                                            gs_max = Inf, gs_regex = NULL, nthread = 2)
  difftime(Sys.time(), start_mast_gsea)
  saveRDS(object = mast_gsea_sd1, file = here::here("intermediate/mast/mpa_sd_compare/mast_gsea_mpa_vs_cmono_sd1_29Jan2024.rds"))
  
  start_mast_gsea <- Sys.time()
  mast_gsea_sd3 <- seu_mast_gsea(mast_dge_result = seu_mast_sd3[["media"]][["cMo_MPA"]],
                                            seu_mast_sets = c("C2_CP_Reactome","C5_GO_BP",
                                                              "C5_GO_MF","H_Hallmark",
                                                              "C2_CP_KEGG_LEGACY","C2_CP_KEGG_MEDICUS", 
                                                              "C7_ImmuneSigDB","C7_VAX"),
                                            num_boots = 50, gs_min = 3, prepare_plot_n = 20,
                                            gs_max = Inf, gs_regex = NULL, nthread = 2)
  difftime(Sys.time(), start_mast_gsea)
  saveRDS(object = mast_gsea_sd3, file = here::here("intermediate/mast/mpa_sd_compare/mast_gsea_mpa_vs_cmono_sd3_29Jan2024.rds"))
}