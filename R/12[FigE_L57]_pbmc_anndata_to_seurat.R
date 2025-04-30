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
write.csv(x = fmat2, file = "/media/MPEdge16/MM137/sc/py/py_out/pbmc/anndata_elements/mm_pbmc_named_cluster_cell_counts.csv")

uclus <- unique(seu$subset_cluster)
deseq_out <- vector("list", length = length(uclus)); names(deseq_out) <- uclus
for(i in 1:length(deseq_out)) {
  print(paste0("starting on [",uclus[i],"] which is [",i,"] of [",length(deseq_out),"] clusters at: ",Sys.time()))
  seu@meta.data$deseq_tmp <- ifelse(seu@meta.data$subset_cluster==names(deseq_out)[i], names(deseq_out)[i], "other")
  deseq_out[[i]] <- seutools::seurat_dge(seurat_object = seu, dge_method = "pseudobulk", mast_lane = NULL, 
                                         assay = "RNA", freq_expressed = 0.1, fc_threshold = log2(1.5), 
                                         test_clusters = "all", cluster_column = "deseq_tmp", 
                                         category_column = "study_day", test_categories = c("SD1","SD3"), 
                                         test_condition = "all", condition_column = "study_day", 
                                         pid_column = "study_id", pseudobulk_test_mode = "cluster_identity", 
                                         return_all_pseudobulk = TRUE)
}
saveRDS(object = deseq_out, file = "/media/MPEdge16/MM137/sc/py/py_out/pbmc/anndata_elements/pbmc_pseudobulk_res.rds")

seu_tnk <- subset(x = seu, cells = grep(pattern = "^TNK", x = seu$subset_cluster))
seu_myl <- subset(x = seu, cells = grep(pattern = "^M", x = seu$subset_cluster))
seu_b <- subset(x = seu, cells = grep(pattern = "^B", x = seu$subset_cluster))

# tnk
uclus <- unique(seu_tnk$subset_cluster)
deseq_out <- vector("list", length = length(uclus)); names(deseq_out) <- uclus
for(i in 1:length(deseq_out)) {
  print(paste0("starting on [",uclus[i],"] which is [",i,"] of [",length(deseq_out),"] clusters at: ",Sys.time()))
  seu_tnk@meta.data$deseq_tmp <- ifelse(seu_tnk@meta.data$subset_cluster==names(deseq_out)[i], names(deseq_out)[i], "other")
  deseq_out[[i]] <- seutools::seurat_dge(seurat_object = seu_tnk, dge_method = "pseudobulk", 
                                         assay = "RNA", freq_expressed = 0.1, fc_threshold = log2(1.5), 
                                         test_clusters = "all", cluster_column = "deseq_tmp", 
                                         category_column = "study_day", test_categories = c("SD1","SD3"), 
                                         test_condition = "all", condition_column = "study_day", 
                                         pid_column = "study_id", pseudobulk_test_mode = "cluster_identity", 
                                         return_all_pseudobulk = TRUE)
}
saveRDS(object = deseq_out, file = "/media/MPEdge16/MM137/sc/py/py_out/pbmc/anndata_elements/pbmc_tnk_pseudobulk_res.rds")

# myeloid
uclus <- unique(seu_myl$subset_cluster)
deseq_out <- vector("list", length = length(uclus)); names(deseq_out) <- uclus
for(i in 1:length(deseq_out)) {
  print(paste0("starting on [",uclus[i],"] which is [",i,"] of [",length(deseq_out),"] clusters at: ",Sys.time()))
  seu_myl@meta.data$deseq_tmp <- ifelse(seu_myl@meta.data$subset_cluster==names(deseq_out)[i], names(deseq_out)[i], "other")
  deseq_out[[i]] <- seutools::seurat_dge(seurat_object = seu_myl, dge_method = "pseudobulk", 
                                         assay = "RNA", freq_expressed = 0.1, fc_threshold = log2(1.5), 
                                         test_clusters = "all", cluster_column = "deseq_tmp", 
                                         category_column = "study_day", test_categories = c("SD1","SD3"), 
                                         test_condition = "all", condition_column = "study_day", 
                                         pid_column = "study_id", pseudobulk_test_mode = "cluster_identity", 
                                         return_all_pseudobulk = TRUE)
}
saveRDS(object = deseq_out, file = "/media/MPEdge16/MM137/sc/py/py_out/pbmc/anndata_elements/pbmc_myeloid_pseudobulk_res.rds")

# b
uclus <- unique(seu_b$subset_cluster)
deseq_out <- vector("list", length = length(uclus)); names(deseq_out) <- uclus
for(i in 1:length(deseq_out)) {
  print(paste0("starting on [",uclus[i],"] which is [",i,"] of [",length(deseq_out),"] clusters at: ",Sys.time()))
  seu_b@meta.data$deseq_tmp <- ifelse(seu_b@meta.data$subset_cluster==names(deseq_out)[i], names(deseq_out)[i], "other")
  deseq_out[[i]] <- seutools::seurat_dge(seurat_object = seu_b, dge_method = "pseudobulk", 
                                         assay = "RNA", freq_expressed = 0.1, fc_threshold = log2(1.5), 
                                         test_clusters = "all", cluster_column = "deseq_tmp", 
                                         category_column = "study_day", test_categories = c("SD1","SD3"), 
                                         test_condition = "all", condition_column = "study_day", 
                                         pid_column = "study_id", pseudobulk_test_mode = "cluster_identity", 
                                         return_all_pseudobulk = TRUE)
}
saveRDS(object = deseq_out, file = "/media/MPEdge16/MM137/sc/py/py_out/pbmc/anndata_elements/pbmc_bcell_pseudobulk_res.rds")

# volc
volcano_genes <- c("CD3E","CD8A","CD4","CD14","FCGR3A","NKG7","TRGV9","TRAV1-2","CD79A","JCHAIN","CXCR4","COCH","TCL1A","YBX3","BANK1","IGHA2","GZMH","GNLY","CCR7","LEF1","MKI67",
                   "TMSB10","TRAC","CCL5","RTKN2","FOXP3","IL2RA","TIGIT","CTLA4","LINC02446","CD3D","KLRD1","LILRA4","AXL","CLEC9A","IDO1","FCER1A","CLEC10A","CD1C","ITM2C","PLD4",
                   "SERPINF1","IRF4","S100A9","S100A8","LYZ","CDKN1C","FCER1G","TYROBP","PRF1","KLRF1","TOP2A","XCL1","XCL2","NCAM1","CD34","GATA2","KIT","SOX4","TNFRSF4","PPBP",
                   "PF4","NRFGN","GNG11","TUBB1","PTPN3","NUCB2","CAV1","MIR4422HG","TRDC","IL7R","CXCR6","CXCR4","GZMK","NCR3",)

tnk_dge <- readRDS(file = "/media/MPEdge16/MM137/sc/py/py_out/pbmc/anndata_elements/pbmc_tnk_pseudobulk_res.rds")
tnk_dge <- lapply(X = tnk_dge, FUN = function(arg1){
  df <- arg1[[names(arg1[grep(pattern = "(^TNK|^M|^B)", x = names(arg1))])]][["res"]]
  df$cluster <- names(arg1[grep(pattern = "(^TNK|^M|^B)", x = names(arg1))])
  df <- df[order(abs(df$log2FoldChange), decreasing = TRUE),]
  return(df)
})
tnk_dge <- tnk_dge[order(as.numeric(gsub(pattern = "(^TNK|^M|^B)_", replacement = "", x = names(tnk_dge))))]

kir_genes <- do.call(rbind, tnk_dge)$gene[grep(pattern = "^sKIR", x = do.call(rbind, tnk_dge)$gene)]; nk_genes <- c("NKG7","NCAM1","CD8A","XCL1",kir_genes,"CD3E","CD3D","CD3G")

tnk_volc <- lapply(X = tnk_dge, FUN = seutools::seu_plot_volcano, prio_top_genes = 5, de_method = "pseudobulk_py", gene_set = volcano_genes)
pdf(file = "/media/MPEdge16/MM137/sc/py/py_out/pbmc/anndata_elements/pbmc_tnk_pseudobulk_volcano.pdf", width = 16, height = 10.25)
lapply(X = tnk_volc, FUN = function(x) x[[1]])
dev.off()

tnk_kir_volc <- lapply(X = tnk_dge, FUN = seutools::seu_plot_volcano, prio_top_genes = 5, de_method = "pseudobulk_py", gene_set = nk_genes)
pdf(file = "/media/MPEdge16/MM137/sc/py/py_out/pbmc/anndata_elements/pbmc_tnk_KIR_pseudobulk_volcano.pdf", width = 16, height = 10.25)
lapply(X = tnk_kir_volc, FUN = function(x) x[[1]])
dev.off()

m_dge <- readRDS(file = "/media/MPEdge16/MM137/sc/py/py_out/pbmc/anndata_elements/pbmc_myeloid_pseudobulk_res.rds")
m_dge <- lapply(X = m_dge, FUN = function(arg1){
  df <- arg1[[names(arg1[grep(pattern = "(^TNK|^M|^B)", x = names(arg1))])]][["res"]]
  df$cluster <- names(arg1[grep(pattern = "(^TNK|^M|^B)", x = names(arg1))])
  df <- df[order(abs(df$log2FoldChange), decreasing = TRUE),]
  return(df)
})
m_dge <- m_dge[order(as.numeric(gsub(pattern = "(^TNK|^M|^B)_", replacement = "", x = names(m_dge))))]
m_volc <- lapply(X = m_dge, FUN = seutools::seu_plot_volcano, prio_top_genes = 5, de_method = "pseudobulk_py", gene_set = volcano_genes)
pdf(file = "/media/MPEdge16/MM137/sc/py/py_out/pbmc/anndata_elements/pbmc_myeloid_pseudobulk_volcano.pdf", width = 16, height = 10.25)
lapply(X = m_volc, FUN = function(x) x[[1]])
dev.off()

b_dge <- readRDS(file = "/media/MPEdge16/MM137/sc/py/py_out/pbmc/anndata_elements/pbmc_bcell_pseudobulk_res.rds")
b_dge <- lapply(X = b_dge, FUN = function(arg1){
  df <- arg1[[names(arg1[grep(pattern = "(^TNK|^M|^B)", x = names(arg1))])]][["res"]]
  df$cluster <- names(arg1[grep(pattern = "(^TNK|^M|^B)", x = names(arg1))])
  df <- df[order(abs(df$log2FoldChange), decreasing = TRUE),]
  return(df)
})
b_dge <- b_dge[order(as.numeric(gsub(pattern = "(^TNK|^M|^B)_", replacement = "", x = names(b_dge))))]
b_volc <- lapply(X = b_dge, FUN = seutools::seu_plot_volcano, prio_top_genes = 5, de_method = "pseudobulk_py", gene_set = volcano_genes)
pdf(file = "/media/MPEdge16/MM137/sc/py/py_out/pbmc/anndata_elements/pbmc_bcell_pseudobulk_volcano.pdf", width = 16, height = 10.25)
lapply(X = b_volc, FUN = function(x) x[[1]])
dev.off()

b_gene <- unique(c('MS4A1', 'TNFRSF13B', 'IGHM', 'IGHD', 'AIM2', 'CD79A', 'LINC01857', 'RALGPS2', 'BANK1', 'CD79B', # Bint
                   'MS4A1', 'COCH', 'AIM2', 'BANK1', 'SSPN', 'CD79A', 'TEX9', 'RALGPS2', 'TNFRSF13C', 'LINC01781', # Bmem
                   'IGHM', 'IGHD', 'CD79A', 'IL4R', 'MS4A1', 'CXCR4', 'BTG1', 'TCL1A', 'CD79B', 'YBX3', # Bnaive
                   'IGHA2', 'MZB1', 'TNFRSF17', 'DERL3', 'TXNDC5', 'TNFRSF13B', 'POU2AF1', 'CPNE5', 'HRASLS2', 'NT5DC2')) # plasmablast
b_b_volc <- lapply(X = b_dge, FUN = seutools::seu_plot_volcano, prio_top_genes = 0, de_method = "pseudobulk_py", gene_set = b_gene)
pdf(file = "/media/MPEdge16/MM137/sc/py/py_out/pbmc/anndata_elements/pbmc_bcell_b_pseudobulk_volcano.pdf", width = 16, height = 10.25)
lapply(X = b_b_volc, FUN = function(x) x[[1]])
dev.off()


if(F) { # run CIPR cluster ID by reference agreement of top markers
  library(SeuratData)
  library(CIPR)
  
  uclus <- unique(seu$subset_cluster)
  wilcox_dge <- vector("list", length = length(uclus)); names(wilcox_dge) <- uclus
  for(i in 1:length(wilcox_dge)) {
    print(paste0("starting on [",uclus[i],"] which is [",i,"] of [",length(wilcox_dge),"] clusters at: ",Sys.time()))
    seu@meta.data$tmp_dge <- ifelse(seu@meta.data$subset_cluster==uclus[i],uclus[i],"other")
    Idents(seu) <- seu@meta.data$tmp_dge
    wilcox_dge[[i]] <- SeuratWrappers::RunPresto(object = seu, assay = "RNA", ident.1 = uclus[i], ident.2 = "other",
                                                 test.use = "wilcox", only.pos = FALSE, min.pct = 0.05)
  }
  wilcox_dge <- lapply(X = wilcox_dge, FUN = function(x){
    x$gene <- row.names(x)
    return(x)
  })
  for(i in 1:length(wilcox_dge)) {
    wilcox_dge[[i]]$cluster <- names(wilcox_dge)[i]
  }
  allmarks <- do.call(rbind, wilcox_dge)
  
  seu_small <- subset(x = seu, cells = which(seu$subset_cluster %in% c("TNK_19","B_6")))
  Idents(seu_small) <- seu_small$subset_cluster
  allmarks <- Seurat::FindAllMarkers(object = seu_small, assay = "RNA")
  
  cipr_ref <- "hsrnaseq"
  cipr_method <- "logfc_spearman"
  
  cipr_out <- CIPR(input_dat = allmarks,
                   comp_method = cipr_method, 
                   reference = cipr_ref, 
                   select_ref_subsets = c("B cell","CD4+ T cell","CD8+ T cell","Dendritic cell","Monocyte","NK cell","Progenitor","MAIT-gdT"), 
                   top_num = 1, 
                   plot_ind = T,
                   plot_top = F, 
                   global_results_obj = T, 
                   global_plot_obj = T
                   # axis.text.x=element_text(color="red") # arguments to pass to ggplot2::theme() to change plotting parameters
  )
  
  new_cipr_plt <- lapply(X = ind_clu_plots, FUN = function(x) return(x + theme(panel.grid.major = element_line(colour = "grey", 
                                                                                                               linewidth = 0.25), 
                                                                               panel.grid.minor = element_line(colour = "grey", 
                                                                                                               linetype = "dashed", 
                                                                                                               linewidth = 0.25))))
  pdf(file = paste0("/media/MPEdge16/MM137/sc/py/py_out/pbmc/anndata_elements/CIPR_",cipr_ref,"_",cipr_method,"_plots.pdf"), 
      width = 25, height = 7)
  lapply(X = new_cipr_plt, FUN = function(x) x)
  dev.off()
}


# myeloid scvi dge
myl_dge <- read.csv(file = "/media/MPEdge16/MM137/sc/py/py_out/pbmc/pbmc_myeloid_scvi_dge_cluster.csv", check.names = FALSE)

# b scvi dge
bcell_dge <- read.csv(file = "/media/MPEdge16/MM137/sc/py/py_out/pbmc/pbmc_b_scvi_dge_cluster.csv", check.names = FALSE)
bcell_dge[bcell_dge$gene=="IGHM",]
bcell_dge[bcell_dge$gene=="IGHD",]
bcell_dge[bcell_dge$gene=="CD79A",]
bcell_dge[bcell_dge$gene=="IL4R",]
bcell_dge[bcell_dge$gene=="MS4A1",]
bcell_dge[bcell_dge$gene=="CXCR4",]
bcell_dge[bcell_dge$gene=="BTG1",]
bcell_dge[bcell_dge$gene=="TCL1A",]
bcell_dge[bcell_dge$gene=="CD79B",]
bcell_dge[bcell_dge$gene=="YBX3",]


