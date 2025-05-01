rm(list = ls()); gc()

# library(seutools)
library(Seurat)
library(SeuratWrappers)
library(Matrix)
dir.create(path = "/media/WD24/sadie10x/rna/mast/pbmc_mpa_mono_sd", showWarnings = FALSE)
setwd("/media/WD24/sadie10x/rna/mast/pbmc_mpa_mono_sd")

seu <- readRDS(file = "/media/WD24/sadie10x/rna/object/MM148_pbmc_seurat.rds")
seu$cell_type <- ifelse(seu$cell_type=="M-platelet", "MPA", seu$cell_type)
seu$cell_type <- ifelse(seu$cell_type=="CD14 Mono", "cMono", seu$cell_type)
seu$cell_type <- ifelse(seu$cell_type=="CD16 Mono", "nMono", seu$cell_type)

# PANEL A (1,2,3) if(F)-wrapped
if(F) { # myeloid singlet umap
  library(ggplot2)
  library(ggrepel)
  library(shadowtext)
  
  myeloid_obs <- read.csv(file = "/media/MPEdge16/MM137/sc/py/py_out/pbmc/anndata_elements/adata_pbmc_myeloid_obs.csv", check.names = FALSE)
  myeloid_umap <- read.csv(file = "/media/MPEdge16/MM137/sc/py/py_out/pbmc/anndata_elements/adata_pbmc_myeloid_umap_coordinates.csv", check.names = FALSE, header = FALSE)
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
  ggsave(filename = "pbmc_myeloid_singlet_umap_16APR2025.pdf", plot = myl_umap, device = "pdf", path = "/media/WD24/sadie10x/rna/figure", 
         width = 8, height = 8, units = "in", dpi = 300, limitsize = F, bg = "white")
}

if(F) { # myeloid doublet umap
  library(ggplot2)
  library(ggrepel)
  library(shadowtext)
  
  myeloid_obs <- read.csv(file = '/media/MPEdge16/MM137/sc/py/py_out/pbmc_myeloid_platelet_int_dbl_obs_to_r.csv', check.names = FALSE)
  myeloid_umap <- read.csv(file = '/media/MPEdge16/MM137/sc/py/py_out/pbmc_myeloid_platelet_int_dbl_umap_coordinates_to_r.csv', check.names = FALSE, header = FALSE)
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
  ggsave(filename = "pbmc_myeloid_doublet_umap_16APR2025.pdf", plot = myl_umap, device = "pdf", path = "/media/WD24/sadie10x/rna/figure", 
         width = 8, height = 8, units = "in", dpi = 300, limitsize = F, bg = "white")
}

if(F) { # myeloid sim, mpa
  library(ggplot2)
  library(ggrepel)
  library(shadowtext)
  
  myeloid_obs <- read.csv(file = '/media/MPEdge16/MM137/sc/py/py_out/pbmc/anndata_elements/adata_pbmc_obs_mpa_sim_2.csv', check.names = FALSE)
  myeloid_umap <- read.csv(file = '/media/MPEdge16/MM137/sc/py/py_out/pbmc/anndata_elements/pbmc_mpa_sim_umap_coordinates_2.csv', 
                           check.names = FALSE, header = FALSE, row.names = NULL)
  colnames(myeloid_umap) <- c("UMAP1","UMAP2")
  
  # sim_obs <- read.csv('/media/MPEdge16/MM137/sc/py/py_out/pbmc/check_bc_sim.csv')
  # mean(myeloid_obs$barcode_2 %in% sim_obs$barcode_2)
  # myeloid_obs_singlet <- read.csv(file = "/media/MPEdge16/MM137/sc/py/py_out/pbmc/anndata_elements/adata_pbmc_myeloid_obs.csv", check.names = FALSE)
  # mean(myeloid_obs_singlet$barcode_2 %in% sim_obs$barcode_2)
  
  # myeloid_obs$droplet_type <- ifelse(myeloid_obs$droplet_type=="MPA_real", "MPA", myeloid_obs$droplet_type)
  # myeloid_obs$droplet_type <- ifelse(myeloid_obs$droplet_type=="MPA_sim", "SimMPA", myeloid_obs$droplet_type)
  # myeloid_obs$droplet_type <- ifelse(myeloid_obs$droplet_type=="singlet", "other", myeloid_obs$droplet_type)
  # myeloid_obs$droplet_type <- ifelse(myeloid_obs$barcode_2 %in% cmo_bc, 'cMono', myeloid_obs$droplet_type)
  # myeloid_obs$droplet_type <- ifelse(myeloid_obs$cell_type=='Platelet', 'Platelet', myeloid_obs$droplet_type)
  myeloid_obs$sim_cell_type <- ifelse(myeloid_obs$sim_cell_type=="M-platelet", "MPA", myeloid_obs$sim_cell_type)
  myeloid_obs$sim_cell_type <- ifelse(myeloid_obs$sim_cell_type=="CD14 Mono", "cMono", myeloid_obs$sim_cell_type)
  myeloid_obs$sim_cell_type <- ifelse(myeloid_obs$sim_cell_type=="CD16 Mono", "nMono", myeloid_obs$sim_cell_type)
  
  myeloid_umap$cluster <- factor(myeloid_obs$sim_cell_type)
  # myl_back <- myeloid_umap[myeloid_umap$cluster=='other',]
  # myl_front <- myeloid_umap[myeloid_umap$cluster!='other',]
  
  # uclus <- unique(myl_front$cluster)
  # 
  # set.seed(123); myl_front <- myl_front[sample(1:nrow(myl_front),nrow(myl_front),replace=F),]
  # 
  # clusx <- rep(NA, length=length(uclus)); names(clusx) <- uclus; clusy <- clusx
  # 
  # for(i in 1:length(uclus)) {
  #   clusx[i] <- median(myl_front$UMAP1[myl_front$cluster==names(clusx)[i]])
  #   clusy[i] <- median(myl_front$UMAP2[myl_front$cluster==names(clusy)[i]])
  # }
  # anno_df <- data.frame(xval = clusx, yval = clusy, lab = names(clusx)); anno_df$lab <- factor(anno_df$lab)
  uclus <- unique(myeloid_umap$cluster)
  set.seed(123); myl_plot_df <- myeloid_umap[sample(1:nrow(myeloid_umap),nrow(myeloid_umap),replace=F),]
  clusx <- rep(NA, length=length(uclus)); names(clusx) <- uclus; clusy <- clusx
  for(i in 1:length(uclus)) {
    clusx[i] <- median(myl_plot_df$UMAP1[myl_plot_df$cluster==names(clusx)[i]])
    clusy[i] <- median(myl_plot_df$UMAP2[myl_plot_df$cluster==names(clusy)[i]])
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
                  'SimMPA' = 'navy')
  
  myl_umap <- ggplot() + 
    ggrastr::geom_point_rast(data = myl_plot_df, aes(x = UMAP1, y = UMAP2, color = cluster, fill = cluster), alpha = 0.4) + 
    annotate("text", x = anno_df$xval, y = anno_df$yval, label = anno_df$lab,
             hjust = 0.5, color = "black", fontface = "bold", size = 9) +
    scale_color_manual(values = custom_col) + 
    scale_fill_manual(values = custom_col) + 
    theme_bw() + 
    theme(legend.text = element_text(size = 18, face = 'bold'), 
          legend.position = 'none', 
          legend.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          axis.title = element_text(size = 24))
  ggsave(filename = "pbmc_myeloid_mpa_simmpa_umap_16APR2025_alt.pdf", plot = myl_umap, device = "pdf", path = "/media/WD24/sadie10x/rna/figure", 
         width = 8, height = 8, units = "in", dpi = 300, limitsize = F, bg = "white")
}

if(F) { # myeloid sim, mpa
  library(ggplot2)
  library(ggrepel)
  library(shadowtext)
  
  myeloid_obs <- read.csv(file = '/media/MPEdge16/MM137/sc/py/py_out/pbmc/anndata_elements/adata_pbmc_obs_mpa_sim.csv', check.names = FALSE)
  myeloid_umap <- read.csv(file = '/media/MPEdge16/MM137/sc/py/py_out/pbmc/anndata_elements/pbmc_mpa_sim_umap_coordinates.csv', 
                           check.names = FALSE, header = FALSE, row.names = NULL)
  colnames(myeloid_umap) <- c("UMAP1","UMAP2")
  myeloid_umap$barcode <- myeloid_obs$barcode_2
  
  myeloid_mapping <- read.csv(file = '/media/MPEdge16/MM137/sc/py/py_out/pbmc/anndata_elements/adata_pbmc_obs_mpa_sim_2.csv')
  myeloid_mapping$sim_cell_type <- ifelse(myeloid_mapping$sim_cell_type=="M-platelet", "MPA", myeloid_mapping$sim_cell_type)
  myeloid_mapping$sim_cell_type <- ifelse(myeloid_mapping$sim_cell_type=="CD14 Mono", "cMono", myeloid_mapping$sim_cell_type)
  myeloid_mapping$sim_cell_type <- ifelse(myeloid_mapping$sim_cell_type=="CD16 Mono", "nMono", myeloid_mapping$sim_cell_type)
  myeloid_map <- myeloid_mapping$sim_cell_type; names(myeloid_map) <- myeloid_mapping$barcode_2
  
  # sim_obs <- read.csv('/media/MPEdge16/MM137/sc/py/py_out/pbmc/check_bc_sim.csv')
  # mean(myeloid_obs$barcode_2 %in% sim_obs$barcode_2)
  # myeloid_obs_singlet <- read.csv(file = "/media/MPEdge16/MM137/sc/py/py_out/pbmc/anndata_elements/adata_pbmc_myeloid_obs.csv", check.names = FALSE)
  # mean(myeloid_obs_singlet$barcode_2 %in% sim_obs$barcode_2)
  
  # myeloid_obs$droplet_type <- ifelse(myeloid_obs$droplet_type=="MPA_real", "MPA", myeloid_obs$droplet_type)
  # myeloid_obs$droplet_type <- ifelse(myeloid_obs$droplet_type=="MPA_sim", "SimMPA", myeloid_obs$droplet_type)
  # myeloid_obs$droplet_type <- ifelse(myeloid_obs$droplet_type=="singlet", "other", myeloid_obs$droplet_type)
  # myeloid_obs$droplet_type <- ifelse(myeloid_obs$barcode_2 %in% cmo_bc, 'cMono', myeloid_obs$droplet_type)
  # myeloid_obs$droplet_type <- ifelse(myeloid_obs$cell_type=='Platelet', 'Platelet', myeloid_obs$droplet_type)
  # myeloid_obs$sim_cell_type <- ifelse(myeloid_obs$sim_cell_type=="M-platelet", "MPA", myeloid_obs$sim_cell_type)
  # myeloid_obs$sim_cell_type <- ifelse(myeloid_obs$sim_cell_type=="CD14 Mono", "cMono", myeloid_obs$sim_cell_type)
  # myeloid_obs$sim_cell_type <- ifelse(myeloid_obs$sim_cell_type=="CD16 Mono", "nMono", myeloid_obs$sim_cell_type)
  
  # myeloid_umap$cluster <- factor(myeloid_obs$sim_cell_type)
  myeloid_umap$cluster <- myeloid_map[myeloid_umap$barcode]
  myeloid_umap$cluster[grep(pattern = "\\-dbl$", x = myeloid_umap$barcode)] <- "SimMPA"
  myeloid_umap$cluster <- factor(myeloid_umap$cluster)
  
  myl_back <- myeloid_umap#[myeloid_umap$cluster=='other',]
  myl_front <- myeloid_umap#[myeloid_umap$cluster!='other',]
  
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
  # saveRDS(object = custom_col, file = file.path("/media/WD24/sadie10x/rna/figure","MPA_project_color_palette.rds"))
  
  myl_umap <- ggplot() +  
    ggrastr::geom_point_rast(data = myl_front, aes(x = UMAP1, y = UMAP2, color = cluster, fill = cluster), alpha = 0.4) +
    annotate("text", x = anno_df$xval, y = anno_df$yval, label = anno_df$lab,
             hjust = 0.5, color = "black", fontface = "bold", size = 9) +
    scale_color_manual(values = custom_col) + 
    scale_fill_manual(values = custom_col) + 
    theme_bw() + 
    theme(legend.text = element_text(size = 18, face = 'bold'), 
          legend.position = 'none', 
          legend.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          axis.title = element_text(size = 24))
  ggsave(filename = "pbmc_myeloid_sim_umap_16APR2025.pdf", plot = myl_umap, device = "pdf", path = "/media/WD24/sadie10x/rna/figure", 
         width = 8, height = 8, units = "in", dpi = 300, limitsize = F, bg = "white")
}

# PANEL G
# MAST SD1 and SD3 (separately) MPA vs cMo
platelet_mast_genes <- read.csv(file = "/media/WD24/sadie10x/rna/mast/platelet_genes/mast_platelet_gene_result.csv", check.names = FALSE)
blacklist_platelet_genes <- platelet_mast_genes$primerid[intersect(which(platelet_mast_genes$fdr<=0.05), which(platelet_mast_genes$logFC>0))]

dir.create("/media/WD24/sadie10x/rna/mast/mpa_sd_compare", showWarnings = FALSE)

seu_test_m <- subset(x = seu, subset = cell_type %in% c("cMono", "MPA"))
seu_test_sd1 <- subset(x = seu_test_m, subset = study_day == "SD1")
seu_test_sd1$condition_custom <- "media"
seu_test_sd1$cluster_custom <- "cMo_MPA"

start_mast <- Sys.time()
seu_mast_sd1 <- seutools::seurat_dge(seurat_object = seu_test_sd1,
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

# saveRDS(object = seu_mast_sd1, file = "/media/WD24/sadie10x/rna/mast/mpa_sd_compare/mast_dge_mpa_vs_cmono_sd1_29Jan2024.rds")
write.csv(x = seu_mast_sd1$media$cMo_MPA$raw_res, file = '/media/WD24/sadie10x/rna/mast/mpa_sd_compare/mast_dge_mpa_vs_cmono_sd1_29Jan2024.csv', row.names = F)

seu_test_sd3 <- subset(x = seu_test_m, subset = study_day == "SD3")
seu_test_sd3$condition_custom <- "media"
seu_test_sd3$cluster_custom <- "cMo_MPA"

start_mast <- Sys.time()
seu_mast_sd3 <- seutools::seurat_dge(seurat_object = seu_test_sd3,
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

saveRDS(object = seu_mast_sd3, file = "/media/WD24/sadie10x/rna/mast/mpa_sd_compare/mast_dge_mpa_vs_cmono_sd3_29Jan2024.rds")
write.csv(x = seu_mast_sd3$media$cMo_MPA$raw_res, file = '/media/WD24/sadie10x/rna/mast/mpa_sd_compare/mast_dge_mpa_vs_cmono_sd3_29Jan2024.csv', row.names = F)
