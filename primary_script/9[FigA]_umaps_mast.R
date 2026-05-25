rm(list = ls()); gc()

library(here)
source(here::here("primary_dependents/seurat_dge_local.R"))
library(Seurat)
library(SeuratWrappers)
library(Matrix)
library(ggplot2)
dir.create(here::here("intermediate/mast/pbmc_mpa_mono_sd"), showWarnings = FALSE, recursive = TRUE)
dir.create(here::here("output/figures"), showWarnings = FALSE, recursive = TRUE)

seu_cache_path <- here::here("intermediate/seurat/MM148_pbmc_seurat.rds")

if(!file.exists(seu_cache_path)) { # build cache on first run; reuse it on subsequent runs
  myeloid_req <- c(
    "intermediate/pbmc/anndata_elements/adata_pbmc_myeloid_counts.mtx",
    "intermediate/pbmc/anndata_elements/adata_pbmc_myeloid_obs.csv",
    "intermediate/pbmc/anndata_elements/adata_pbmc_myeloid_var.csv"
  )
  platelet_req <- c(
    "intermediate/pbmc/anndata_elements/adata_pbmc_platelet_counts.mtx",
    "intermediate/pbmc/anndata_elements/adata_pbmc_platelet_obs.csv",
    "intermediate/pbmc/anndata_elements/adata_pbmc_platelet_var.csv"
  )
  req_paths <- c(myeloid_req, platelet_req)
  req_full <- vapply(req_paths, here::here, character(1))
  if(!all(file.exists(req_full))) {
    missing_paths <- req_paths[!file.exists(req_full)]
    stop(
      paste0(
        "Missing myeloid/platelet anndata elements needed to rebuild Seurat object:\n",
        paste0(" - ", missing_paths, collapse = "\n"),
        "\n\nWrite the required myeloid and platelet anndata elements first, then rerun this script."
      )
    )
  }

  build_seu_from_elements <- function(counts_path, obs_path, var_path, prefix_name = "obj") {
    ct <- Matrix::t(Matrix::readMM(file = here::here(counts_path)))
    obs <- read.csv(file = here::here(obs_path), check.names = FALSE)
    var <- read.csv(file = here::here(var_path), check.names = FALSE)

    if(!"barcode_2" %in% colnames(obs)) {
      if("barcode_2.1" %in% colnames(obs)) {
        obs$barcode_2 <- obs$barcode_2.1
      } else {
        stop(paste0("Could not find barcode_2 in ", prefix_name, " obs file."))
      }
    }

    if(!"cell_type" %in% colnames(obs)) {
      if("merged_type" %in% colnames(obs)) {
        obs$cell_type <- as.character(obs$merged_type)
      } else {
        stop(paste0("Could not find cell_type or merged_type in ", prefix_name, " obs file."))
      }
    }
    if(!"merged_type" %in% colnames(obs)) {
      obs$merged_type <- as.character(obs$cell_type)
    }

    required_meta <- c("barcode_2", "study_id", "study_day", "cell_type")
    missing_meta <- required_meta[!required_meta %in% colnames(obs)]
    if(length(missing_meta) > 0) {
      stop(paste0("Missing required metadata columns in ", prefix_name, " obs file: ", paste(missing_meta, collapse = ", ")))
    }

    obs$barcode_2 <- as.character(obs$barcode_2)
    obs$study_id <- as.character(obs$study_id)
    obs$study_day <- as.character(obs$study_day)
    obs$cell_type <- as.character(obs$cell_type)
    obs$merged_type <- as.character(obs$merged_type)

    genes <- as.character(var[[1]])
    if(any(is.na(genes)) || any(genes == "")) {
      stop(paste0("Missing or blank gene names in ", prefix_name, " var file."))
    }
    if(anyDuplicated(genes) > 0) {
      stop(paste0("Duplicate gene names in ", prefix_name, " var file."))
    }
    if(any(is.na(obs$barcode_2)) || any(obs$barcode_2 == "")) {
      stop(paste0("Missing or blank barcodes in ", prefix_name, " obs file."))
    }
    if(any(is.na(obs$study_id)) || any(obs$study_id == "")) {
      stop(paste0("Missing or blank study_id values in ", prefix_name, " obs file."))
    }
    if(any(is.na(obs$study_day)) || any(obs$study_day == "")) {
      stop(paste0("Missing or blank study_day values in ", prefix_name, " obs file."))
    }
    if(any(is.na(obs$cell_type)) || any(obs$cell_type == "")) {
      stop(paste0("Missing or blank cell_type values in ", prefix_name, " obs file."))
    }
    if(anyDuplicated(obs$barcode_2) > 0) {
      stop(paste0("Duplicate barcodes in ", prefix_name, " obs file."))
    }
    if(nrow(ct) != length(genes)) {
      stop(paste0("Gene count mismatch in ", prefix_name, ": nrow(counts)=", nrow(ct), " vs length(var)=", length(genes)))
    }
    if(ncol(ct) != nrow(obs)) {
      stop(paste0("Cell count mismatch in ", prefix_name, ": ncol(counts)=", ncol(ct), " vs nrow(obs)=", nrow(obs)))
    }

    row.names(ct) <- genes
    colnames(ct) <- obs$barcode_2
    row.names(obs) <- obs$barcode_2
    obs <- obs[colnames(ct), , drop = FALSE]

    seu_obj <- Seurat::CreateSeuratObject(counts = ct, assay = "RNA", meta.data = obs)
    return(seu_obj)
  }

  seu_myeloid <- build_seu_from_elements(
    counts_path = "intermediate/pbmc/anndata_elements/adata_pbmc_myeloid_counts.mtx",
    obs_path = "intermediate/pbmc/anndata_elements/adata_pbmc_myeloid_obs.csv",
    var_path = "intermediate/pbmc/anndata_elements/adata_pbmc_myeloid_var.csv",
    prefix_name = "myeloid"
  )

  seu_platelet <- build_seu_from_elements(
    counts_path = "intermediate/pbmc/anndata_elements/adata_pbmc_platelet_counts.mtx",
    obs_path = "intermediate/pbmc/anndata_elements/adata_pbmc_platelet_obs.csv",
    var_path = "intermediate/pbmc/anndata_elements/adata_pbmc_platelet_var.csv",
    prefix_name = "platelet"
  )

  common_genes <- intersect(rownames(seu_myeloid), rownames(seu_platelet))
  if(length(common_genes) == 0) {
    stop("No common genes across myeloid and platelet anndata elements.")
  }
  common_genes <- sort(common_genes)

  seu_myeloid <- subset(seu_myeloid, features = common_genes)
  seu_platelet <- subset(seu_platelet, features = common_genes)

  if(!identical(rownames(seu_myeloid), rownames(seu_platelet))) {
    stop("Gene order mismatch between myeloid and platelet Seurat objects after feature intersection.")
  }

  seu <- merge(seu_myeloid, y = seu_platelet, merge.data = FALSE)
  if(ncol(seu) != (ncol(seu_myeloid) + ncol(seu_platelet))) {
    stop("Merged Seurat object has unexpected cell count.")
  }

  head(seu@assays$RNA@layers$counts@x,20) # should be ints
  seu <- Seurat::NormalizeData(object = seu, assay = "RNA", normalization.method = "LogNormalize")
  
  seu$cell_type <- ifelse(seu$cell_type=="M-platelet", "MPA", seu$cell_type)
  seu$cell_type <- ifelse(seu$cell_type=="CD14 Mono", "cMono", seu$cell_type)
  seu$cell_type <- ifelse(seu$cell_type=="CD16 Mono", "nMono", seu$cell_type)
  Seurat::Idents(seu) <- seu$cell_type
  saveRDS(object = seu, file = seu_cache_path)
} else { # read cached object on re-run
  seu <- readRDS(file = seu_cache_path)
  # format names
  seu$cell_type <- ifelse(seu$cell_type=="M-platelet", "MPA", seu$cell_type)
  seu$cell_type <- ifelse(seu$cell_type=="CD14 Mono", "cMono", seu$cell_type)
  seu$cell_type <- ifelse(seu$cell_type=="CD16 Mono", "nMono", seu$cell_type)
}

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
                'doublet' = '#A9A9A9',
                'SimMPA' = '#d50000')

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

myeloid_obs_sim <- read.csv(file = here::here("intermediate/pbmc/anndata_elements/adata_pbmc_obs_mpa_sim.csv"), check.names = FALSE)
myeloid_umap_sim <- read.csv(file = here::here("intermediate/pbmc/anndata_elements/pbmc_mpa_sim_umap_coordinates.csv"),
                             check.names = FALSE, header = FALSE, row.names = NULL)
colnames(myeloid_umap_sim) <- c("UMAP1","UMAP2")

myeloid_obs_sim$sim_celltype <- myeloid_obs_sim$merged_type
myeloid_obs_sim$sim_celltype <- ifelse(myeloid_obs_sim$sim_celltype=="M-platelet", "MPA", myeloid_obs_sim$sim_celltype)
myeloid_obs_sim$sim_celltype <- ifelse(myeloid_obs_sim$sim_celltype=="cMo_Platelet_doublet", "SimMPA", myeloid_obs_sim$sim_celltype)
myeloid_obs_sim$sim_celltype <- ifelse(myeloid_obs_sim$sim_celltype=="CD14 Mono", "cMono", myeloid_obs_sim$sim_celltype)
myeloid_obs_sim$sim_celltype <- ifelse(myeloid_obs_sim$sim_celltype=="CD16 Mono", "nMono", myeloid_obs_sim$sim_celltype)

myeloid_umap_sim$cluster <- factor(myeloid_obs_sim$sim_celltype)

uclus <- unique(myeloid_umap_sim$cluster)
# ensure shuffled for random draw order
set.seed(123); myl_plot_df <- myeloid_umap_sim[sample(1:nrow(myeloid_umap_sim),nrow(myeloid_umap_sim),replace=F),]
clusx <- rep(NA, length=length(uclus)); names(clusx) <- uclus; clusy <- clusx
for(i in 1:length(uclus)) {
  clusx[i] <- median(myl_plot_df$UMAP1[myl_plot_df$cluster==names(clusx)[i]])
  clusy[i] <- median(myl_plot_df$UMAP2[myl_plot_df$cluster==names(clusy)[i]])
}
anno_df <- data.frame(xval = clusx, yval = clusy, lab = names(clusx)); anno_df$lab <- factor(anno_df$lab)

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
ggsave(filename = "pbmc_myeloid_sim_umap_16APR2025.pdf", plot = myl_umap, device = "pdf", path = here::here("output/figures"), 
       width = 8, height = 8, units = "in", dpi = 300, limitsize = F, bg = "white")

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
                                       condition_column = "condition_custom", 
                                       pid_column = "study_id",
                                       pseudobulk_test_mode = "cluster_by_category",
                                       filter_genes = "outer"), silent = TRUE)
  mast_res[[i]] <- seu_mast
}

mast_res_df <- do.call(rbind, lapply(X = mast_res, FUN = function(arg1) return(arg1[[1]][[1]][["raw_res"]]))) # if full object does not save out.. Error: C stack usage  7972852 is too close to the limit

saveRDS(object = mast_res_df, file = here::here("intermediate/mast/pbmc_mpa_mono_sd/pbmc_mpa_mono_mast_deg_sd1sd3_15JAN2025.rds"))

dir.create(here::here("intermediate/mast/platelet_genes"), showWarnings = FALSE, recursive = TRUE)

seu_test_plt <- subset(x = seu, subset = cell_type %in% c("cMono", "Platelet"))
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
                           test_categories = c("cMono","Platelet"),
                           test_condition = "all",
                           condition_column = "condition_custom", 
                           pid_column = "study_id",
                           pseudobulk_test_mode = "cluster_by_category",
                           filter_genes = "outer")
difftime(Sys.time(), start_mast, units = "mins")

write.csv(x = seu_mast_plt[["media"]][["cMo_Platelet"]][["filtered_res"]][order(seu_mast_plt[["media"]][["cMo_Platelet"]][["filtered_res"]]$logFC,decreasing=T),], 
          file = here::here("intermediate/mast/platelet_genes/mast_platelet_gene_result.csv"), row.names = FALSE)

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
                                     condition_column = "condition_custom", 
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