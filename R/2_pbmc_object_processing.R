rm(list = ls()); gc()

sink.reset <- function(){
  for(i in seq_len(sink.number())){
    sink(NULL)
  }
}

stop_print <- function() {
  on.exit(sink.reset())
  print("incorrect package versions for Seurat and/or SeuratObject")
}

sink(file = paste0("/media/MPEdge16/MM137/sc/seurat/logging/pbmc_integration_console_redirect_",strftime(Sys.time(),"%Y-%m-%d_%H%M%S"),".txt"), 
     type = "output")

library(Seurat) # https://satijalab.org/seurat/articles/announcements.html // https://satijalab.org/seurat/articles/essential_commands
library(SeuratWrappers)
library(Matrix)

setwd("/media/MPEdge16/MM137/sc/seurat")

print(sessionInfo())

seurat_version <- paste0(unlist(packageVersion("Seurat")), collapse = ".")
seuratobject_version <- paste0(unlist(packageVersion("SeuratObject")), collapse = ".")
if(seurat_version!="5.0.1" || seuratobject_version!="5.0.1") {
  stop_print()
  stop()
}

sink.reset()

# options(future.globals.maxSize=Inf)
# future::plan(multisession, workers = 20) # https://github.com/HenrikBengtsson/future/issues/420

cr_h5_fil <- list.files(path = "/media/MPEdge16/MM137/sc/10x_cloud_dl", pattern = "filtered_feature_bc_matrix.h5", 
                        all.files = TRUE, full.names = TRUE, recursive = TRUE, include.dirs = TRUE)
cr_h5_fil <- cr_h5_fil[!grepl(pattern = "_standard", x = cr_h5_fil)]

soc_fil <- list.files(path = "/media/MPEdge16/MM137/sc/souporcell_outs", pattern = "clusters\\.tsv$", 
                      all.files = TRUE, full.names = TRUE, recursive = TRUE, include.dirs = TRUE)
solo_fil <- list.files(path = "/media/MPEdge16/MM137/sc/py/py_out", pattern = "solo_scores", 
                       all.files = TRUE, full.names = TRUE, recursive = TRUE, include.dirs = TRUE)
snp_fil <- list.files(path = "/media/MPEdge16/MM137/sc/snp_correlation", pattern = "\\.txt$", 
                      all.files = TRUE, full.names = TRUE, recursive = TRUE, include.dirs = TRUE)
project_key <- read.csv(file = "/media/MPEdge16/MM137/sc/MM148_key.csv", check.names = FALSE)
project_key$Lane <- paste0("MM-",project_key$`Multiplex ID`)
project_key$`Sample ID` <- gsub(pattern = "SADIE", replacement = "Sadie-", x = project_key$`Sample ID`)

lanes <- paste0("MM-",1:8)

qc_metrics <- vector("list", length = length(lanes)); names(qc_metrics) <- lanes
seu_list <- qc_metrics

for(i in 1:length(cr_h5_fil)) {
  obj <- Seurat::Read10X_h5(filename = cr_h5_fil[grep(pattern = paste0(lanes[i],"/"), x = cr_h5_fil)])
  soc <- read.delim(file = soc_fil[grep(pattern = paste0(lanes[i],"_"), x = soc_fil)])
  snp <- read.delim(file = snp_fil[grep(pattern = paste0("\\/",gsub("-","",lanes[i]),"\\/"), x = snp_fil)])
  snp <- snp[snp$Cluster_ID!="unassigned",]
  snp$Genotype_ID <- stringr::str_extract(string = snp$Genotype_ID, pattern = "Sadie-[0-9]+")
  lane_key <- project_key[project_key$Lane==lanes[i],]
  solo <- read.csv(file = solo_fil[grep(pattern = paste0("solo_scores_",gsub("-","",lanes[i])), x = solo_fil)])
  colnames(solo) <- paste0("solo_",colnames(solo))
  num_solo_dbl <- sum(solo$solo_prediction=="doublet")
  seu <- CreateSeuratObject(counts = obj, assay = "RNA")
  seu <- PercentageFeatureSet(object = seu, pattern = "^MT\\-", col.name = "percent_mito", assay = "RNA")
  
  num_bc <- ncol(seu)
  num_soc_dbl <- sum(soc$status=="doublet"); rate_soc_dbl <- mean(soc$status[soc$status %in% c("singlet","doublet")]=="doublet") * 100
  num_solo_dbl <- sum(solo$solo_prediction=="doublet"); rate_solo_dbl <- (num_solo_dbl/sum(soc$status %in% c("singlet","doublet"))) * 100
  soc_assignments <- table(soc$assignment)
  
  seu <- AddMetaData(object = seu, metadata = colnames(seu), col.name = "barcode_1")
  seu <- AddMetaData(object = seu, metadata = rep(lanes[i], ncol(seu)), col.name = "lane")
  seu <- subset(x = seu, subset = barcode_1 %in% solo$solo_X)
  if(mean(solo$solo_X==colnames(seu))!=1) {
    stop("mismatch in barcode order; reorder solo prediction df..")
  } else {
    seu@meta.data <- cbind(seu@meta.data, solo)
  }
  pid <- rep(NA,ncol(seu)); names(pid) <- colnames(seu); day <- pid
  for(j in 1:nrow(snp)) {
    pid[which(names(pid) %in% soc$barcode[soc$assignment==snp$Cluster_ID[j]])] <- snp$Genotype_ID[j]
  }
  if(sum(is.na(pid))!=0) {
    stop("one or more barcodes not mapped to a study ID..")
  }
  for(j in 1:nrow(lane_key)) {
    day[which(pid==lane_key$`Sample ID`[j])] <- lane_key$`Study day`[j]
  }
  seu <- AddMetaData(object = seu, metadata = pid, col.name = "study_id")
  seu <- AddMetaData(object = seu, metadata = day, col.name = "study_day")
  
  seu <- subset(x = seu, subset = solo_prediction == "singlet")
  seu <- subset(x = seu, subset = percent_mito <= 10)
  pid_table <- table(seu$study_id)
  
  qc_metrics[[i]] <- list(cr_filtered_barcodes = num_bc, 
                          souporcell_doublets = num_soc_dbl, 
                          souporcell_doublet_rate = rate_soc_dbl, 
                          solo_doublets = num_solo_dbl, 
                          solo_doublet_rate = rate_solo_dbl, 
                          souporcell_assignment_table = soc_assignments, 
                          study_id_table = pid_table, 
                          thresholds = list(scanpy.pp.filter_cells = "min_genes = 200", 
                                            scanpy.pp.filter_genes = "min_cells = 10", 
                                            mito = "percent_mito <= 10"), 
                          other_methods = list(doublets = list(heterotypic = "souporcell", 
                                                               homotypic = "scvi-tools solo"), 
                                               pid_assignment = "souporcell snp_Illumina LCG Assay MEGA chip snp correlation"), 
                          date_time = Sys.time())
  
  colnames(seu) <- gsub(pattern = "\\-[0-9]+$", replacement = gsub(pattern = "MM", replacement = "", x = lanes[i]), x = colnames(seu))
  seu <- AddMetaData(object = seu, metadata = colnames(seu), col.name = "barcode_2")
  
  seu_list[[i]] <- seu
}
saveRDS(object = qc_metrics, file = "/media/MPEdge16/MM137/sc/seurat/pbmc/all_cells/mm_pbmc_lane_qc_metrics.rds")
saveRDS(object = seu_list, file = "/media/MPEdge16/MM137/sc/seurat/pbmc/all_cells/mm_pbmc_lane_list.rds")
if(F) {
  seu_list <- readRDS(file = "/media/MPEdge16/MM137/sc/seurat_analysis/pbmc/all_cells/mm_pbmc_lane_list.rds")
}

merged_obj <- merge(x = seu_list[[1]], y = seu_list[2:length(seu_list)])
merged_obj_join <- JoinLayers(object = merged_obj)
saveRDS(object = merged_obj_join, file = "/media/MPEdge16/MM137/sc/seurat/pbmc/all_cells/mm_pbmc_merged.rds")
# seurat to anndata: /media/MPEdge16/MM137/sc/seurat/pbmc/all_cells/mm_pbmc_lane_list.rds

counts_matrix <- GetAssayData(object = merged_obj_join, assay = "RNA", layer = "counts")
obj_meta <- merged_obj_join@meta.data
gene_names <- row.names(counts_matrix)
# barcode_names <- colnames(merged_obj_join)

writeMM(obj = counts_matrix, file = "/media/MPEdge16/MM137/sc/seurat/pbmc/RNA_assay_counts.mtx")
write.csv(x = obj_meta, file = "/media/MPEdge16/MM137/sc/seurat/pbmc/seurat_meta.csv", quote = FALSE, row.names = FALSE)
write.table(x = data.frame(gene = gene_names), file = "/media/MPEdge16/MM137/sc/seurat/pbmc/gene_names.csv", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

rm(list = ls()); gc()
q()
