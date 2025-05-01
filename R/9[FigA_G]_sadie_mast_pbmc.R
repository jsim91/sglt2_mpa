# Mona slack Mon 9 DEC 2024: So I think sticking with SD1-3 for now is great, and for 
# clusters let's start with the MPA, CD14+ monos, CD16+ monos to start. Tiger likes 
# the Hallmark gene sets and it was sort of interesting in the past, so we can include 
# that one. KEGG and Reactome are always good. I don't know much about C7 immunologic 
# signature gene sets, but sounds like it could be relevant!
rm(list = ls()); gc()

# library(seutools)
library(Seurat)
library(SeuratWrappers)
library(Matrix)
# dir.create(path = "/media/WD24/sadie10x/rna/mast", showWarnings = FALSE)
dir.create(path = "/media/WD24/sadie10x/rna/mast/pbmc_mpa_mono_sd", showWarnings = FALSE)
setwd("/media/WD24/sadie10x/rna/mast/pbmc_mpa_mono_sd")

if(F) {
  ct <- Matrix::t(Matrix::readMM(file = '/media/MPEdge16/MM137/sc/py/py_out/pbmc/anndata_elements/adata_pbmc_counts.mtx'))
  obs <- read.csv(file = "/media/MPEdge16/MM137/sc/py/py_out/pbmc/anndata_elements/adata_pbmc_obs.csv")
  var <- read.csv(file = "/media/MPEdge16/MM137/sc/py/py_out/pbmc/anndata_elements/adata_pbmc_var.csv")
  umap_coord <- read.csv(file = "/media/MPEdge16/MM137/sc/py/py_out/pbmc/anndata_elements/adata_pbmc_umap_coordinates.csv", header = FALSE)
  colnames(umap_coord) <- c("UMAP1","UMAP2"); row.names(umap_coord) <- obs$barcode_2
  latent_coord <- read.csv(file = "/media/MPEdge16/MM137/sc/py/py_out/pbmc/anndata_elements/adata_pbmc_latent_coordinates.csv", header = FALSE)
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
  saveRDS(object = seu, file = "/media/WD24/sadie10x/rna/object/MM148_pbmc_seurat.rds")
} else {
  seu <- readRDS(file = "/media/WD24/sadie10x/rna/object/MM148_pbmc_seurat.rds")
  seu$cell_type <- ifelse(seu$cell_type=="M-platelet", "MPA", seu$cell_type)
  seu$cell_type <- ifelse(seu$cell_type=="CD14 Mono", "cMono", seu$cell_type)
  seu$cell_type <- ifelse(seu$cell_type=="CD16 Mono", "nMono", seu$cell_type)
}

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
  
  # gg_color_hue <- function(n) {
  #   hues = seq(15, 375, length = n + 1)
  #   hcl(h = hues, l = 65, c = 100)[1:n]
  # }
  # custom_col <- gg_color_hue(6); names(custom_col) <- uclus
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
    # geom_point(data = myeloid_umap, aes(x = UMAP1, y = UMAP2, fill = cluster), alpha = 0.9, pch = 21, color = "black", stroke = 0.1) + 
    # shadowtext::geom_shadowtext(data = anno_df, aes(x = xval, y = yval, label = lab, color = lab), 
    #                             size = 8, bg.color = "black", bg.r = 0.1) + 
    annotate("text", x = anno_df$xval, y = anno_df$yval, label = anno_df$lab,
             hjust = 0.5, color = "black", fontface = "bold", size = 9) + 
    scale_color_manual(values = custom_col) +
    # scale_fill_manual(values = custom_col) +
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
    # geom_point(data = myeloid_umap, aes(x = UMAP1, y = UMAP2, fill = cluster), alpha = 0.9, pch = 21, color = "black", stroke = 0.1) + 
    # shadowtext::geom_shadowtext(data = anno_df, aes(x = xval, y = yval, label = lab, color = lab), 
    #                             size = 8, bg.color = "black", bg.r = 0.1) + 
    annotate("text", x = anno_df$xval, y = anno_df$yval, label = anno_df$lab,
             hjust = 0.5, color = "black", fontface = "bold", size = 9) + 
    scale_color_manual(values = custom_col) +
    # scale_fill_manual(values = custom_col) +
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
    # ggrastr::geom_point_rast(data = myl_front, aes(x = UMAP1, y = UMAP2, color = cluster, fill = cluster), alpha = 0.4) + 
    # ggrastr::geom_point_rast(data = myl_back, aes(x = UMAP1, y = UMAP2), color = "black", alpha = 0.2) +
    # geom_point(data = myl_front, aes(x = UMAP1, y = UMAP2, fill = cluster), alpha = 0.9, pch = 21, color = "black", stroke = 0.1) + 
    # shadowtext::geom_shadowtext(data = anno_df, aes(x = xval, y = yval, label = lab, color = lab), 
    #                             size = 8, bg.color = "black", bg.r = 0.1) + 
    annotate("text", x = anno_df$xval, y = anno_df$yval, label = anno_df$lab,
             hjust = 0.5, color = "black", fontface = "bold", size = 9) +
    scale_color_manual(values = custom_col) + 
    scale_fill_manual(values = custom_col) + 
    # guides(color = guide_legend(override.aes = list(size = 7.5, pch = 21, stroke = 0.2, color = "black", alpha = 1))) + 
    theme_bw() + 
    theme(legend.text = element_text(size = 18, face = 'bold'), 
          # legend.position = "bottom", 
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
    # ggrastr::geom_point_rast(data = myl_plot_df, aes(x = UMAP1, y = UMAP2, color = cluster, fill = cluster), alpha = 0.4) + 
    ggrastr::geom_point_rast(data = myl_front, aes(x = UMAP1, y = UMAP2, color = cluster, fill = cluster), alpha = 0.4) +
    # ggrastr::geom_point_rast(data = myl_back, aes(x = UMAP1, y = UMAP2), color = "black", alpha = 0.2) +
    # geom_point(data = myl_front, aes(x = UMAP1, y = UMAP2, fill = cluster), alpha = 0.9, pch = 21, color = "black", stroke = 0.1) + 
    # shadowtext::geom_shadowtext(data = anno_df, aes(x = xval, y = yval, label = lab, color = lab), 
    #                             size = 8, bg.color = "black", bg.r = 0.1) + 
    annotate("text", x = anno_df$xval, y = anno_df$yval, label = anno_df$lab,
             hjust = 0.5, color = "black", fontface = "bold", size = 9) +
    scale_color_manual(values = custom_col) + 
    scale_fill_manual(values = custom_col) + 
    # guides(color = guide_legend(override.aes = list(size = 7.5, pch = 21, stroke = 0.2, color = "black", alpha = 1))) + 
    theme_bw() + 
    theme(legend.text = element_text(size = 18, face = 'bold'), 
          # legend.position = "bottom", 
          legend.position = 'none', 
          legend.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          axis.title = element_text(size = 24))
  ggsave(filename = "pbmc_myeloid_sim_umap_16APR2025.pdf", plot = myl_umap, device = "pdf", path = "/media/WD24/sadie10x/rna/figure", 
         width = 8, height = 8, units = "in", dpi = 300, limitsize = F, bg = "white")
}

if(F) {
  dir.create("/media/WD24/sadie10x/rna/mast/pbmc_mpa_cmono_by_sd", showWarnings = FALSE)
  
  test_day <- c("SD1","SD3")
  
  seu@meta.data$barcode <- seu$barcode_2
  
  mast_res <- vector("list", length = length(test_day)); names(mast_res) <- test_day; mast_gsea_res <- mast_res
  
  for(i in 1:length(mast_res)) {
    start_mast_i <- Sys.time()
    print(paste0("starting: ",names(mast_res)[i]," [",i,"] of [",length(mast_res),"] clusters"))
    
    seu_test <- subset(x = seu, subset = cell_type %in% c("cMono","MPA"))
    
    target_cluster <- names(mast_res)[i]
    seu_test$condition_custom <- "media"
    seu_test$category_custom <- seu_test$cell_type
    seu_test$cluster_custom <- "cell_group"
    
    mast_seu <- subset(x = seu_test, subset = study_day == test_day[i]) # subset for categories
    
    try(seu_mast <- seutools::seurat_dge(seurat_object = mast_seu,
                                         dge_method = "mast",
                                         assay = "RNA",
                                         freq_expressed = 0.1,
                                         fc_threshold = log2(1.25),
                                         test_clusters = "cell_group",
                                         cluster_column = "cluster_custom",
                                         category_column = "category_custom",
                                         test_categories = c("cMono","MPA"),
                                         test_condition = "all",
                                         condition_column = "condition_custom", # pseudo condition to prevent errors and customize output list name(s); cells were not perturbed, therefore going with 'media'
                                         pid_column = "study_id",
                                         pseudobulk_test_mode = "cluster_by_category",
                                         filter_genes = "outer"), silent = TRUE)
    mast_res[[i]] <- seu_mast
    
    mast_gsea_res[[i]] <- seutools:::seu_mast_gsea(mast_dge_result = seu_mast[['media']][['cell_group']],
                                                   seu_mast_sets = c("C2_CP_Reactome","C5_GO_BP",
                                                                     "C5_GO_MF","H_Hallmark",
                                                                     "C2_CP_KEGG_LEGACY","C2_CP_KEGG_MEDICUS", 
                                                                     "C7_ImmuneSigDB","C7_VAX"),
                                                   num_boots = 50, gs_min = 3, prepare_plot_n = 20,
                                                   gs_max = Inf, gs_regex = NULL, nthread = 2)
    
    print(paste0("finished: ",names(mast_res)[i]," in ",round(as.numeric(difftime(Sys.time(), start_mast_i, units = "mins")),1)," mins"))
  }
  
  mast_res_df <- do.call(rbind, lapply(X = mast_res, FUN = function(arg1) return(arg1[[1]][[1]][["raw_res"]]))) # if full object does not save out.. Error: C stack usage  7972852 is too close to the limit
  
  saveRDS(object = mast_res_df, file = "/media/WD24/sadie10x/rna/mast/pbmc_mpa_cmono_by_sd/pbmc_mpa_cmono_mast_deg_sd1_sd3_31JAN2025.rds")
  
  saveRDS(object = mast_gsea_res, file = "/media/WD24/sadie10x/rna/mast/pbmc_mpa_cmono_by_sd/pbmc_mpa_cmono_mast_gsea_sd1_sd3_31JAN2025.rds")
}

install.packages(c("dplyr", "igraph", "ggplot2", "future", "future.apply", "pbapply", "irlba", "NMF", "ggalluvial", "stringr", "svglite", "Matrix", "ggrepel", "circlize", "RColorBrewer", "cowplot", "methods", "ComplexHeatmap", "RSpectra", "Rcpp", "reticulate", "scales", "sna", "reshape2", "FNN", "shape", "BiocGenerics", "magrittr", "patchwork", "colorspace", "plyr", "ggpubr", "ggnetwork", "BiocNeighbors", "plotly", "shiny", "bslib"))


if(F) {
  dir.create("/media/WD24/sadie10x/rna/figure/feature_violins", showWarnings = FALSE)
  
  feature_str <- c('CCR2, Slc9a1, Slc5a2, GP9, ITGA2B, ITGB3, F2RL3, PPBP, IL1B, IL1R1, CD40LG, CD40, SELP, SELPLG, CCL5, CX3CR1, CX3CL1, PF4, F2R, F3, ITGAM, ITGB2')
  feature_set <- toupper(strsplit(x = feature_str, split = ", ")[[1]])
  if(!"seu" %in% ls()) {
    seu <- readRDS(file = "/media/WD24/sadie10x/rna/object/MM148_pbmc_seurat.rds")
  }
  feature_set <- feature_set[which(feature_set %in% row.names(seu))]
  
  # Seurat::VlnPlot(object = seu, features = feature_set, pt.size = 1, alpha = 0, idents = "cell_type", assay = "RNA", )
  seu_adj <- seu
  target_names <- c('cMono', 'nMono', 'MPA', 'platelets')
  seu_adj$cell_type <- ifelse(seu_adj$cell_type=="CD14 Mono", "cMono", seu_adj$cell_type)
  seu_adj$cell_type <- ifelse(seu_adj$cell_type=="CD16 Mono", "nMono", seu_adj$cell_type)
  seu_adj$cell_type <- ifelse(seu_adj$cell_type=="M-platelet", "MPA", seu_adj$cell_type)
  seu_adj$cell_type <- ifelse(seu_adj$cell_type=="Platelet", "platelets", seu_adj$cell_type)
  seu_adj$cell_type <- ifelse(seu_adj$cell_type %in% target_names, seu_adj$cell_type, "other")
  seu_adj$cell_type <- factor(seu_adj$cell_type, levels = c('cMono', 'nMono', 'MPA', 'platelets', 'other'))
  Idents(seu_adj) <- seu_adj$cell_type
  seu_vln <- Seurat::VlnPlot(object = seu_adj, features = feature_set, alpha = 0, assay = "RNA", layer = "data")
  seu_vln <- lapply(X = seu_vln, FUN = function(arg1) return(arg1 + geom_boxplot(color = "black", fill = scales::alpha("black",0), width = 0.05, outlier.shape = NA) + theme(axis.title.x = element_blank())))
  seu_vln[[15]]# + geom_boxplot(color = "black", fill = scales::alpha("black",0), width = 0.05, outlier.shape = NA)
  
  pdf(file = "/media/WD24/sadie10x/rna/figure/feature_violins/drug_target_aggregation_violins_20JAN2025.pdf", width = 4, height = 4)
  lapply(X = seu_vln, FUN = function(x) x)
  dev.off()
  
  seu$long_id <- gsub("\\-","",paste0(seu$study_id,"_",seu$study_day))
  seu_aggr <- Seurat::AggregateExpression(object = seu, assays = "RNA", features = feature_set, 
                                          group.by = c("long_id","cell_type"), return.seurat = TRUE)
  
  # avg_cluster_expression <- Seurat::AggregateExpression(object = seu, assays = "RNA", features = feature_set, 
  #                                                       group.by = c("long_id","cell_type"), return.seurat = TRUE)
  rown <- avg_cluster_expression$RNA@Dimnames[[1]]
  coln <- avg_cluster_expression$RNA@Dimnames[[2]]
  ace <- as.data.frame(avg_cluster_expression)
  colnames(ace) <- gsub("RNA\\.","",coln)
  row.names(ace) <- rown
  ace_melt <- reshape2::melt(ace)
  colnames(ace_melt) <- c("variable","pseudocount")
  ace_melt$study_id <- stringr::str_extract(string = ace_melt$variable, pattern = "Sadie[0-9]+")
  ace_melt$study_day <- stringr::str_extract(string = ace_melt$variable, pattern = "SD[0-9]")
  ace_melt$cluster <- gsub("_","",stringr::str_extract(ace_melt$variable,"_.+$"))
}

if(F) {
  # check expression of SLC9A1 in platelets, MPA, CD14 Mono, CD16 Mono split by SD1 SD2 SD3
  # SLC9A1 may be a SADIE drug target
  library(Seurat); library(ggplot2)
  
  # Define the gene of interest
  gene_of_interest <- "SLC9A1"
  target_cell_type <- "CD16 Mono" # one of: "Platelet", "M-platelet", "CD14 Mono", "CD16 Mono"
  
  # subset for cell group of interest
  seu_subset <- subset(x = seu, subset = cell_type == target_cell_type)
  
  # Extract meta data and expression data
  meta_data <- seu_subset@meta.data
  expression_data <- GetAssayData(seu_subset, assay = "RNA", layer = "data")
  
  # Initialize data frames to store results
  percent_expressed <- data.frame(patient = character(), time_point = character(), 
                                  percent_expressed = numeric(), stringsAsFactors = FALSE)
  average_expression <- data.frame(patient = character(), time_point = character(), 
                                   average_expression = numeric(), stringsAsFactors = FALSE)
  
  # Loop through each patient and time point to calculate % expressed and average expression
  patients <- unique(meta_data$study_id); patients <- patients[order(as.numeric(gsub("Sadie-","",patients)))]
  time_points <- unique(meta_data$study_day); time_points <- time_points[order(as.numeric(gsub("SD","",time_points)))]
  
  for (patient in patients) {
    for (time_point in time_points) {
      # Subset data for the specific patient and time point
      cells <- rownames(meta_data)[meta_data$study_id == patient & meta_data$study_day == time_point]
      if (length(cells) > 0) {
        expr_values <- expression_data[gene_of_interest, cells]
        
        # Calculate % expressed
        pct_expressed <- sum(expr_values > 0) / length(expr_values) * 100
        
        # Calculate average expression
        avg_expr <- mean(expr_values)
        
        # Append results to data frames
        percent_expressed <- rbind(percent_expressed, data.frame(patient = patient, time_point = time_point, 
                                                                 percent_expressed = pct_expressed, stringsAsFactors = FALSE))
        average_expression <- rbind(average_expression, data.frame(patient = patient, time_point = time_point, 
                                                                   average_expression = avg_expr, stringsAsFactors = FALSE))
      }
    }
  }
  
  # Combine results into one data frame
  results <- merge(percent_expressed, average_expression, by = c("patient", "time_point"))
  
  # Convert time_point to factor with levels in chronological order
  results$time_point <- factor(results$time_point, levels = c("SD1", "SD2", "SD3"))
  
  # Plot % expressed
  ggplot(results, aes(x = time_point, y = percent_expressed)) +
    geom_boxplot() +
    geom_line(aes(group = patient), alpha = 0.3) +
    labs(title = "Percentage of Cells Expressing SLC9A1", x = "Time Point", y = "% Expressed")
  
  # Plot average expression
  ggplot(results, aes(x = time_point, y = average_expression)) +
    geom_boxplot() +
    geom_line(aes(group = patient), alpha = 0.3) +
    labs(title = "Average Expression Level of SLC9A1", x = "Time Point", y = "Average Expression")
  
  # Statistical testing
  # Perform ANOVA for % expressed
  percent_anova <- aov(percent_expressed ~ time_point + Error(patient/time_point), data = results)
  summary(percent_anova)
  
  # Perform ANOVA for average expression
  average_anova <- aov(average_expression ~ time_point + Error(patient/time_point), data = results)
  summary(average_anova)
}

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
    
    mast_gsea_res[[i]] <- seutools:::seu_mast_gsea(mast_dge_result = seu_mast[['media']][[target_cluster]],
                                                   seu_mast_sets = c("C2_CP_Reactome","C5_GO_BP",
                                                                     "C5_GO_MF","H_Hallmark",
                                                                     "C2_CP_KEGG_LEGACY","C2_CP_KEGG_MEDICUS", 
                                                                     "C7_ImmuneSigDB","C7_VAX"),
                                                   num_boots = 50, gs_min = 3, prepare_plot_n = 20,
                                                   gs_max = Inf, gs_regex = NULL, nthread = 2)
    
    print(paste0("finished: ",names(mast_res)[i]," in ",round(as.numeric(difftime(Sys.time(), start_mast_i, units = "mins")),1)," mins"))
  }
  
  mast_res_df <- do.call(rbind, lapply(X = mast_res, FUN = function(arg1) return(arg1[[1]][[1]][["raw_res"]]))) # if full object does not save out.. Error: C stack usage  7972852 is too close to the limit
  
  saveRDS(object = mast_res_df, file = "/media/WD24/sadie10x/rna/mast/pbmc_mpa_mono_sd/pbmc_mpa_mono_mast_deg_sd1sd3_15JAN2025.rds")
  
  saveRDS(object = mast_gsea_res, file = "/media/WD24/sadie10x/rna/mast/pbmc_mpa_mono_sd/pbmc_mpa_mono_mast_gsea_sd1sd3_15JAN2025.rds")
}

if(F) { # find genes associated with platelets to drop when doing MPA vs CD14 Mono
  dir.create("/media/WD24/sadie10x/rna/mast/platelet_genes", showWarnings = FALSE)
  
  seu_test_plt <- subset(x = seu, subset = cell_type %in% c("CD14 Mono", "Platelet"))
  seu_test_plt$condition_custom <- "media"
  seu_test_plt$cluster_custom <- "cMo_Platelet"
  
  start_mast <- Sys.time()
  seu_mast_plt <- seutools::seurat_dge(seurat_object = seu_test_plt,
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
  
  save(seu_mast_plt, file = "/media/WD24/sadie10x/rna/mast/platelet_genes/cMo_vs_Platelet_mast.RData")
  write.csv(x = seu_mast_plt[["media"]][["cMo_Platelet"]][["filtered_res"]][order(seu_mast_plt[["media"]][["cMo_Platelet"]][["filtered_res"]]$logFC,decreasing=T),], 
            file = "/media/WD24/sadie10x/rna/mast/platelet_genes/mast_platelet_gene_result.csv", row.names = FALSE)
}

if(T) { # SD1 and SD3 (separately) MPA vs cMo
  # platelet_mast_genes <- read.csv(file = "/media/WD24/sadie10x/rna/mast/platelet_genes/mast_platelet_gene_result.csv", check.names = FALSE)
  # blacklist_platelet_genes <- platelet_mast_genes$primerid[intersect(which(platelet_mast_genes$fdr<=0.05), 
  #                                                                    which(platelet_mast_genes$logFC>0))]
  
  dir.create("/media/WD24/sadie10x/rna/mast/mpa_sd_compare", showWarnings = FALSE)
  
  seu_test_m <- subset(x = seu, subset = cell_type %in% c("CD14 Mono", "M-platelet"))
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
                                       test_categories = c("CD14 Mono","M-platelet"),
                                       test_condition = "all",
                                       condition_column = "condition_custom", # pseudo condition to prevent errors and customize output list name(s); cells were not perturbed, therefore going with 'media'
                                       pid_column = "study_id",
                                       pseudobulk_test_mode = "cluster_by_category",
                                       filter_genes = "outer")#, 
                                       # gene_set_blacklist = blacklist_platelet_genes)
  difftime(Sys.time(), start_mast, units = "mins")
  
  saveRDS(object = seu_mast_sd1, file = "/media/WD24/sadie10x/rna/mast/mpa_sd_compare/mast_dge_mpa_vs_cmono_sd1_29Jan2024.rds")

  seu_test_m <- subset(x = seu, subset = cell_type %in% c("CD14 Mono", "M-platelet"))
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
                                       test_categories = c("CD14 Mono","M-platelet"),
                                       test_condition = "all",
                                       condition_column = "condition_custom", # pseudo condition to prevent errors and customize output list name(s); cells were not perturbed, therefore going with 'media'
                                       pid_column = "study_id",
                                       pseudobulk_test_mode = "cluster_by_category",
                                       filter_genes = "outer")#, 
                                       # gene_set_blacklist = blacklist_platelet_genes)
  difftime(Sys.time(), start_mast, units = "mins")
  
  saveRDS(object = seu_mast_sd3, file = "/media/WD24/sadie10x/rna/mast/mpa_sd_compare/mast_dge_mpa_vs_cmono_sd3_29Jan2024.rds")
            
  start_mast_gsea <- Sys.time()
  mast_gsea_sd1 <- seutools:::seu_mast_gsea(mast_dge_result = seu_mast_sd1[["media"]][["cMo_MPA"]],
                                            seu_mast_sets = c("C2_CP_Reactome","C5_GO_BP",
                                                              "C5_GO_MF","H_Hallmark",
                                                              "C2_CP_KEGG_LEGACY","C2_CP_KEGG_MEDICUS", 
                                                              "C7_ImmuneSigDB","C7_VAX"),
                                            num_boots = 50, gs_min = 3, prepare_plot_n = 20,
                                            gs_max = Inf, gs_regex = NULL, nthread = 2)
  difftime(Sys.time(), start_mast_gsea)
  saveRDS(object = mast_gsea_sd1, file = "/media/WD24/sadie10x/rna/mast/mpa_sd_compare/mast_gsea_mpa_vs_cmono_sd1_29Jan2024.rds")
  
  start_mast_gsea <- Sys.time()
  mast_gsea_sd3 <- seutools:::seu_mast_gsea(mast_dge_result = seu_mast_sd3[["media"]][["cMo_MPA"]],
                                            seu_mast_sets = c("C2_CP_Reactome","C5_GO_BP",
                                                              "C5_GO_MF","H_Hallmark",
                                                              "C2_CP_KEGG_LEGACY","C2_CP_KEGG_MEDICUS", 
                                                              "C7_ImmuneSigDB","C7_VAX"),
                                            num_boots = 50, gs_min = 3, prepare_plot_n = 20,
                                            gs_max = Inf, gs_regex = NULL, nthread = 2)
  difftime(Sys.time(), start_mast_gsea)
  saveRDS(object = mast_gsea_sd3, file = "/media/WD24/sadie10x/rna/mast/mpa_sd_compare/mast_gsea_mpa_vs_cmono_sd3_29Jan2024.rds")
}