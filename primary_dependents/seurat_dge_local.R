# Local copy of core DGE helper from jsim91/seutools (seurat_in.R).
# Sourced locally to avoid requiring seutools:: namespace for this manuscript workflow.
# Companion guide: see pipeline_docs/seurat_dge_companion.md for function behavior and repo-specific usage.

seurat_dge <- function(seurat_object,
                       dge_method = c("mast", "wilcox", "pseudobulk"),
                       assay = "RNA",
                       freq_expressed = 0.1,
                       filter_genes = "inner", # "inner" is more restrictive than "outer"
                       fc_threshold = log2(1.5),
                       test_clusters = "all", # one or more clusters to test, or "all"
                       mast_lane = NULL,
                       cluster_column = "cell_type",
                       category_column = "age_group",
                       test_categories = c("younger","older"), # order matters: translates into (test_categories[1]/test_categories[2]) for DESeq2 pseudobulk
                       test_condition = "all",
                       condition_column = "condition",
                       pid_column = "pid",
                       wilcox_only_positive = FALSE,
                       pseudobulk_test_mode = c("cluster_identity","cluster_by_category","cluster_by_condition"),
                       return_all_pseudobulk = TRUE,
                       gene_set_blacklist = NULL,
                       plot_n = 100) {
  suppressPackageStartupMessages({
    require(Seurat)
  })

  if(all(!is.null(test_categories), pseudobulk_test_mode!="cluster_by_category")) {
    stop("'test_categories' should be set to NULL when 'pseudobulk_test_mode' is not set to 'cluster_by_category'")
  }
  if(test_condition[1]!="all") {
    conditions <- test_condition[test_condition %in% seurat_object@meta.data[,condition_column]]
  } else {
    conditions <- unique(seurat_object@meta.data[,condition_column])
  }
  if(length(conditions)==0) {
    conditions <- NULL
    warning("None of the requested 'conditions' found in seurat_object@meta.data[,condition_column]. Continuing without subsetting on any condition(s).")
  } else {
    which_is_condition <- seurat_object@meta.data[,condition_column] %in% conditions
    seurat_object <- subset(x = seurat_object, cells = which(which_is_condition))
    print(paste0("seurat_object subset for: ",paste0(unique(seurat_object@meta.data[,condition_column]), collapse = ",")))
  }
  seu_conditions <- ifelse(test_condition[1]=="all","all",unique(seurat_object@meta.data[,condition_column]))
  if(length(seu_conditions)!=1) {
    seu_conditions <- paste0(seu_conditions, collapse = ".")
  }
  if(test_clusters[1]!="all") {
    test_clusters <- test_clusters[test_clusters %in% seurat_object@meta.data[,cluster_column]]
    if(length(test_clusters)==0) {
      annos <- NULL
      warning("None of 'test_clusters' found in seurat_object@meta.data[,cluster_column]. Continuing with test(s) on all cells.")
    } else {
      annos <- test_clusters
    }
  } else {
    annos <- unique(seurat_object@meta.data[,cluster_column])
  }
  if(is.null(test_categories)) {
    test_cats <- NULL
  } else if(test_categories[1]!="all") {
    test_categories <- test_categories[test_categories %in% seurat_object@meta.data[,category_column]]
    if(length(test_categories)==0) {
      test_cats <- NULL
      warning("no categories by 'test_categories' were found in seurat_object@meta.data[,category_column]")
    } else {
      test_cats <- test_categories
    }
  } else {
    warning("Setting 'test_categories' to 'all' or using 'all' as a tested category may produce errors or unexpected results. Categories should be explicitly defined.")
    test_categories <- unique(seurat_object@meta.data[,category_column])
    test_cats <- test_categories
  }
  if(all(is.null(annos), is.null(test_cats))) {
    stop("Nothing to compare")
  }

  if(tolower(dge_method)=="pseudobulk") {
    require(data.table)
    require(DESeq2)
    require(ashr)
    
    if(pseudobulk_test_mode!='cluster_identity') {
      seurat_object <- subset(seurat_object, cells = which(seurat_object@meta.data[,category_column] %in% test_cats))
    }
    capture_dir <- system.file(package = "seutools")
    Matrix::writeMM(obj = seurat_object@assays[[assay]]@layers$counts,
                    file = paste0(capture_dir,"/temp_files/__python_ct_matrix__.mtx"))
    write.csv(x = seurat_object@meta.data,
              file = paste0(capture_dir,"/temp_files/__python_obs_matrix__.csv"), row.names = FALSE)
    write.csv(x = data.frame(v1 = row.names(seurat_object)),
              file = paste0(capture_dir,"/temp_files/__python_feature_names__.csv"), row.names = FALSE)
    system(command = paste0("python3 ",
                            paste0(capture_dir,"/python/create_pseudobulk_profile.py")," ",
                            paste0(capture_dir,"/temp_files/__python_ct_matrix__.mtx")," ",
                            capture_dir,"/temp_files"," ",
                            paste0(capture_dir,"/temp_files/__python_obs_matrix__.csv")," ",
                            ifelse(pseudobulk_test_mode=="cluster_identity",cluster_column,category_column)," ",
                            paste0(capture_dir,"/temp_files/__python_feature_names__.csv")," ",
                            condition_column," ",
                            pid_column," ",
                            pseudobulk_test_mode))

    drop_novar2 <- function(arg1) {
      if(ncol(arg1)>nrow(arg1)) {
        check_var <- apply(X = arg1, MARGIN = 2, FUN = stats::var)
        which_novar <- which(check_var==0)
        if(length(which_novar)!=0) {
          arg1 <- arg1[,-which_novar]
        }
      } else {
        check_var <- apply(X = arg1, MARGIN = 1, FUN = stats::var)
        which_novar <- which(check_var==0)
        if(length(which_novar)!=0) {
          arg1 <- arg1[-which_novar,]
        }
      }
      if(ncol(arg1)>nrow(arg1)) {
        arg1 <- t(arg1)
      }
      return(arg1)
    }

    if(pseudobulk_test_mode=="cluster_identity") {
      ct_spl <- vector("list", length = length(unique(seurat_object@meta.data[,cluster_column])))
      names(ct_spl) <- gsub("( |\\-|\\/)","_",unique(seurat_object@meta.data[,cluster_column]))
      meta_list <- vector("list", length = length(ct_spl)); names(meta_list) <- names(ct_spl)
      for(i in 1:length(ct_spl)) {
        ct_spl[[i]] <- as.data.frame(data.table::fread(paste0(capture_dir,"/temp_files/__pseudobulk_sum_counts_", names(ct_spl)[i], "__.csv"), check.names = FALSE, header = FALSE))
        if(file.exists(paste0(capture_dir,"/temp_files/__pseudobulk_sum_counts_", names(ct_spl)[i], "__.csv"))) {
          file.remove(paste0(capture_dir,"/temp_files/__pseudobulk_sum_counts_", names(ct_spl)[i], "__.csv"))
        }
        meta_list[[i]] <- read.csv(paste0(capture_dir,"/temp_files/__pseudobulk_obs_", names(ct_spl)[i], "__.csv"), check.names = FALSE, row.names = 1)
        meta_list[[i]]$cell_group <- row.names(meta_list[[i]])
        if(file.exists(paste0(capture_dir,"/temp_files/__pseudobulk_obs_", names(ct_spl)[i], "__.csv"))) {
          file.remove(paste0(capture_dir,"/temp_files/__pseudobulk_obs_", names(ct_spl)[i], "__.csv"))
        }
        obj_var <- read.csv(paste0(capture_dir,"/temp_files/__pseudobulk_var_", names(ct_spl)[i], "__.csv"))[,1]
        if(file.exists(paste0(capture_dir,"/temp_files/__pseudobulk_var_", names(ct_spl)[i], "__.csv"))) {
          file.remove(paste0(capture_dir,"/temp_files/__pseudobulk_var_", names(ct_spl)[i], "__.csv"))
        }
        row.names(ct_spl[[i]]) <- meta_list[[i]]$cell_group
        colnames(ct_spl[[i]]) <- obj_var
        ct_spl[[i]] <- drop_novar2(ct_spl[[i]])
        ##
        meta_list[[i]]$cluster_id <- meta_list[[i]][,cluster_column]
        meta_list[[i]]$sample_id <- meta_list[[i]][,pid_column]
        meta_list[[i]]$group_id <- factor(meta_list[[i]][,cluster_column], levels = c("in_cluster","out"))
        row.names(meta_list[[i]]) <- meta_list[[i]]$cell_group
        drop_indices <- which(!meta_list[[i]][,pid_column] %in% names(table(meta_list[[i]][,pid_column])[table(meta_list[[i]][,pid_column])==2]))
        if(length(drop_indices)!=0) {
          meta_list[[i]] <- meta_list[[i]][which(!1:nrow(meta_list[[i]]) %in% drop_indices),]
        }
        ct_spl[[i]] <- ct_spl[[i]][,which(colnames(ct_spl[[i]]) %in% meta_list[[i]]$cell_group)]
      }
    } else {
      ct_data <- as.data.frame(data.table::fread(paste0(capture_dir,"/temp_files/__pseudobulk_sum_counts__.csv"), check.names = FALSE, header = FALSE))
      if(file.exists(paste0(capture_dir,"/temp_files/__pseudobulk_sum_counts__.csv"))) {
        file.remove(paste0(capture_dir,"/temp_files/__pseudobulk_sum_counts__.csv"))
      }
      obj_obs <- read.csv(paste0(capture_dir,"/temp_files/__pseudobulk_obs__.csv"), check.names = FALSE, row.names = 1)
      obj_obs$cell_group <- row.names(obj_obs)
      if(length(test_clusters)>1) {
        obj_obs[,cluster_column] <- paste0(test_clusters, collapse=",")
      } else if(length(test_clusters)==1) {
        obj_obs[,cluster_column] <- test_clusters
      } else {
        stop("'test_clusters' should be a single cluster or multiple clusters that will be joined when not testing for cluster identity")
      }
      if(file.exists(paste0(capture_dir,"/temp_files/__pseudobulk_obs__.csv"))) {
        file.remove(paste0(capture_dir,"/temp_files/__pseudobulk_obs__.csv"))
      }
      obj_var <- read.csv(paste0(capture_dir,"/temp_files/__pseudobulk_var__.csv"))[,1]
      if(file.exists(paste0(capture_dir,"/temp_files/__pseudobulk_var__.csv"))) {
        file.remove(paste0(capture_dir,"/temp_files/__pseudobulk_var__.csv"))
      }
      row.names(ct_data) <- obj_obs$cell_group
      colnames(ct_data) <- obj_var

      id_pattern <- paste0(gsub("-","\\\\-",paste0("(",paste0(unique(seurat_object@meta.data[,pid_column]), collapse = "|"),")")),"_")
      ct_spl <- split(x = ct_data, f = obj_obs[,cluster_column])
      ct_spl <- lapply(X = ct_spl, FUN = drop_novar2)
      obj_obs$cluster_id <- obj_obs[,cluster_column]
      obj_obs$sample_id <- obj_obs[,pid_column]
      obj_obs$group_id <- factor(obj_obs[,category_column], levels = test_categories)
      meta_list <- list(obj_obs)
    }

    deseq_input <- vector("list", length = length(ct_spl)); names(deseq_input) <- names(ct_spl)
    for(i in 1:length(deseq_input)) {
      deseq_input[[i]] <- list(ct_spl[[i]], meta_list[[i]])
    }

    use_adj_p <- TRUE # hard coding use adjusted p values; consider allowing use unadjusted for discovery
    padj_threshold <- 0.05 # hard coding 0.05; consider allowing to change
    do_deseq <- function(arg1, p_return_threshold = padj_threshold, stim = seu_conditions,
                         use_adj = use_adj_p) {

      count_data <- arg1[[1]]
      meta_data <- arg1[[2]]

      test_cats_order <- levels(meta_data$group_id)

      count_data <- round(count_data)

      if(nrow(meta_data)==0) {
        fn_flag <- 0
        res <- "Not enough cells from either of the tested groups. Min cells is 10."
      } else {
        group_id <- levels(meta_data$group_id)
        num_g1 <- sum(meta_data$group_id==group_id[1])
        num_g2 <- sum(meta_data$group_id==group_id[2])
        print(paste0("num_g1 is [",num_g1,"]"))
        print(paste0("group_id[1] is [",group_id[1],"]"))
        print(paste0("num_g2 is [",num_g2,"]"))
        print(paste0("group_id[2] is [",group_id[2],"]"))
        print(paste0("[",num_g1,"] data points available for [",group_id[1],"] and [",num_g2,"] data points available for [",group_id[2],"]"))
        if(any(num_g1<=1,num_g2<=1)) {
          fn_flag <- 0
          res <- paste0("[",num_g1,"] from group [",group_id[1],"] and [",num_g2,"] from group [",group_id[2],"] available for testing. Not enough to test. Min cells is 10.")
        } else {
          fn_flag <- 1
          dds <- DESeqDataSetFromMatrix(countData = count_data,
                                        colData = meta_data,
                                        design = ~ group_id)
          dds <- DESeq(dds)
          normalized_counts <- counts(dds, normalized = TRUE)
          contrast <- c("group_id", test_cats_order[1], test_cats_order[2])
          res <- results(dds,
                         contrast = contrast,
                         alpha = 0.05)
          res <- lfcShrink(dds,
                           contrast = contrast,
                           res=res,
                           type = "ashr")
          norm_ct <- DESeq2::counts(dds, normalized = TRUE)
          res <- data.frame(res); res <- res[!is.na(res$padj),]
          if(sum(is.na(res$padj))!=0) {
            res$padj <- p.adjust(p = res$pvalue, method = "fdr")
          }
          res$gene <- row.names(res)
          res$condition <- stim
          res$contrast_ident.1 <- test_cats_order[1]
          res$contrast_ident.2 <- test_cats_order[2]
          if(use_adj) {
            if(!is.null(p_return_threshold)) {
              res <- res[res$padj<=p_return_threshold,]
            }
          } else {
            if(!is.null(p_return_threshold)) {
              res <- res[res$pvalue<=p_return_threshold,]
            }
          }
        }
      }
      if(fn_flag==1) {
        return(list(res = res, norm_ct = norm_ct, meta = meta_data))#, hm = pheatm))
      } else {
        return(list(res = res, meta = meta_data))#, hm = pheatm))
      }
    }
    deseq_out <- lapply(X = deseq_input, FUN = do_deseq, p_return_threshold = padj_threshold,
                        stim = seu_conditions, use_adj = use_adj_p)

    if(return_all_pseudobulk) {
      return(deseq_out)
    } else {
      print("'return_all_pseudobulk' is set to FALSE; trimming non-significant results. Set 'return_all_pseudobulk' to TRUE (recommended) to return significant and non-significant results.")
    }

    drop_indices <- c()
    for(i in 1:length(deseq_out)) {
      if(nrow(deseq_out[[i]][[1]])!=0) {
        deseq_out[[i]][[1]]$cluster <- names(deseq_out)[i]
      } else {
        drop_indices <- append(drop_indices, i)
      }
    }
    if(length(drop_indices)==1) {
      deseq_out <- deseq_out[[-drop_indices]]
    } else if(length(drop_indices)>1) {
      deseq_out <- deseq_out[-drop_indices]
    }

    deseq_res <- lapply(X = deseq_out, FUN = function(arg1) return(arg1[[1]]))
    deseq_write <- do.call(rbind, deseq_res)
    return(deseq_write)
  }

  dge_outs <- vector("list", length = ifelse(is.null(conditions), 1, length(conditions)))
  if(is.null(conditions)) {
    names(dge_outs) <- "all"
  } else {
    names(dge_outs) <- conditions
  }
  dge_outs <- lapply(X = dge_outs, FUN = function(arg1, ann = annos){
    tmpv <- vector("list", length = ifelse(is.null(ann), 1, length(ann)))
    if(is.null(ann)) {
      names(tmpv) <- "all"
    } else {
      names(tmpv) <- ann
    }
    return(tmpv)
  })

  start_test <- Sys.time()
  for(i in 1:length(dge_outs)) {
    start_i <- Sys.time()
    if(names(dge_outs)[i]=="all") {
      subs1 <- seurat_object
    } else {
      seurat_object <- AddMetaData(object = seurat_object, metadata = seurat_object@meta.data[,condition_column], col.name = "c")
      subs1 <- subset(x = seurat_object, subset = c == names(dge_outs)[i])
    }
    for(j in 1:length(dge_outs[[i]])){
      start_j <- Sys.time()
      print(paste0("starting on [",names(dge_outs[[i]])[j],"] in [",names(dge_outs)[i],"] at ",Sys.time()))
      if(!is.null(annos)) {
        if(!is.null(test_cats)) {
          subs1 <- AddMetaData(object = subs1, metadata = subs1@meta.data[,cluster_column], col.name = "l")
          subs2 <- subset(x = subs1, subset = l == names(dge_outs[[i]])[j])
        } else {
          subs1 <- AddMetaData(object = subs1,
                               metadata = ifelse(subs1@meta.data[,cluster_column]==names(dge_outs[[i]])[j], names(dge_outs[[i]])[j], "other"),
                               col.name = "l")
          subs2 <- subs1
        }
      } else {
        subs2 <- subs1
      }
      if(!is.null(test_cats)) {
        ident1 <- test_categories[1]
        ident2 <- test_categories[2]
        Idents(subs2) <- subs2@meta.data[,category_column]
      } else {
        ident1 <- names(dge_outs[[i]])[j] # will be either "all" or a cluster
        if(is.null(ident1)) {
          stop("Nothing to test")
        }
        ident2 <- "other"
        Idents(subs2) <- ifelse(subs2@meta.data[,cluster_column]==ident1, ident1, "other")
      }
      print(paste0("for [",names(dge_outs[[i]])[j],"] in [",names(dge_outs)[i],"] ident.1 = ",ident1,", ident.2 = ",ident2))
      if(tolower(dge_method) == "mast") {
        suppressPackageStartupMessages({
          require(MAST)
          require(data.table)
          require(ggplot2)
        })
        test_idents <- c(ident1, ident2)
        print(paste0("dropping lowly expressed genes. Threshold >= ",freq_expressed*100,"% of cells."))
        keep_genes <- vector("list", length = length(test_idents)); names(keep_genes) <- test_idents; drop_blacklist_genes <- keep_genes
        num_unique_genes <- nrow(subs2)
        for(k in 1:2) {
          if(pseudobulk_test_mode=="cluster_by_category") {
            subs2_grp <- subset(x = subs2, cells = which(subs2@meta.data[,category_column]==test_idents[k]))
          } else if(pseudobulk_test_mode=="cluster_identity") {
            subs2_grp <- subset(x = subs2, cells = which(subs2@meta.data[,"l"]==test_idents[k]))
          }
          bem <- subs2_grp@assays[[assay]]@layers$counts > 0
          percent_expr <- rowSums(bem) / ncol(bem); names(percent_expr) <- row.names(subs2)
          keep_genes[[k]] <- row.names(subs2)[percent_expr >= freq_expressed]
          if(!is.null(gene_set_blacklist)) {
            drop_blacklist_genes[[k]] <- keep_genes[[k]][keep_genes[[k]] %in% gene_set_blacklist]
            keep_genes[[k]] <- keep_genes[[k]][!keep_genes[[k]] %in% gene_set_blacklist]
          } else {
            drop_blacklist_genes[[k]] <- c()
          }
        }
        grp_genes <- unlist(keep_genes); grp_genes <- grp_genes[!is.na(grp_genes)]
        num_blacklist_genes <- length(unique(unlist(drop_blacklist_genes)))
        if(filter_genes=="inner") {
          genes_to_keep <- names(table(grp_genes)[which(table(grp_genes)==2)])
        } else if(filter_genes=="outer") {
          genes_to_keep <- grp_genes
        }
        if(length(genes_to_keep)==0) {
          warning("no genes passed thresholding")
          next
        } else {
          print(paste0("number of dropped genes through blacklisting: ",num_blacklist_genes," genes."))
          print(paste0("total number of genes dropped: ",num_unique_genes - length(genes_to_keep)," genes. ",length(genes_to_keep)," genes left for MAST testing."))
        }
        subs2 <- subset(x = subs2, features = genes_to_keep)
        if(pseudobulk_test_mode=="cluster_by_category") {
          subs2$category <- subs2@meta.data[,category_column]
        } else if(pseudobulk_test_mode=="cluster_identity") {
          subs2$category <- subs2@meta.data[,"l"]
        }
        subs2$category <- ifelse(subs2$category==ident1,"Group1","Group2")
        Idents(subs2) <- subs2$category
        seu_as_sce <- as.SingleCellExperiment(subs2, assay = "RNA")
        #seu_as_sce <- scuttle::logNormCounts(seu_as_sce, log = TRUE)
        #logcounts(seu_as_sce) <- log2(seu_as_sce)
        no_dense <- counts(seu_as_sce)
        nzero_ct <- no_dense@x
        no_dense@x <- log2(nzero_ct+1)
        logcounts(seu_as_sce) <- no_dense # log2 counts without sparse->dense coercion
        cdr <- colSums(seu_as_sce@assays@data@listData[["logcounts"]]>0) # cellular detection rate; number of genes detected which is a well-known confounder (https://bioconductor.org/packages/release/bioc/vignettes/MAST/inst/doc/MAITAnalysis.html#22_Filtering)
        seu_as_sce$cngeneson <- scale(cdr)
        my_sca <- SceToSingleCellAssay(seu_as_sce, check_sanity = FALSE)
        if(!is.null(mast_lane)) {
          latv <- c("cngeneson",pid_column,mast_lane)
        } else {
          latv <- c("cngeneson",pid_column)
        }
        cond<-factor(colData(my_sca)$category)
        cond<-relevel(cond,"Group1")
        colData(my_sca)$category <- cond

        if (!is.null(pid_column)) {
          if (!pid_column %in% latv) {
            stop("Random effect variable (sample ID) should be included in latent variables! Specify sample ID using arg 'pid_column'")
          }
          latv <- latv[!latv %in% pid_column]
          fmla <- as.formula(object = paste0(" ~ category + ", paste(latv, collapse = "+"), glue::glue(" + (1|{pid_column})")))
          print(fmla)
          zlmCond <- MAST::zlm(formula = fmla,
                               sca = my_sca,
                               exprs_value = 'logcounts',
                               method="glmer",
                               ebayes=FALSE,
                               silent=T,
                               fitArgsD = list(nAGQ = 0),
                               strictConvergence = FALSE)

          # fitArgsD = list(nAGQ = 1)
        } else {
          stop("Trying to run without mixed effect model but this is not supported. Sample ID must be specified with arg 'pid_column'")
        }
        summaryCond <- MAST::summary(object = zlmCond, doLRT = 'categoryGroup2')
        summaryDt <- summaryCond$datatable

        fcHurdle <- merge(summaryDt[contrast=='categoryGroup2' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                          summaryDt[contrast=='categoryGroup2' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients

        fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]; fcHurdle_all <- fcHurdle
        fcHurdleSig <- merge(fcHurdle[fdr<.05 & abs(coef)>fc_threshold], data.table::as.data.table(mcols(my_sca)), by='primerid')
        setorder(fcHurdleSig, fdr); setorder(fcHurdle_all, fdr); fcHurdle_all <- fcHurdle_all[!is.na(fcHurdle_all$fdr),]
        mast_res <- as.data.frame(fcHurdleSig)
        mast_res2 <- as.data.frame(fcHurdle_all)

        lfcs <- MAST::logFC(zlmfit = zlmCond)

        lfc1 <- as.data.frame(lfcs$logFC[mast_res$primerid,])
        lfc1$primerid <- row.names(lfc1)

        lfc2 <- as.data.frame(lfcs$logFC[mast_res2$primerid,])
        lfc2$primerid <- row.names(lfc2)

        getlfcs <- as.data.frame(MAST::getLogFC(zlmfit = zlmCond))
        getlfcs <- getlfcs[which(getlfcs$contrast=="categoryGroup2"),]
        row.names(getlfcs) <- getlfcs$primerid

        getlfcs1 <- getlfcs[mast_res$primerid,]
        getlfcs2 <- getlfcs[mast_res2$primerid,]

        if(nrow(mast_res)!=0) {
          mast_res <- merge(x = mast_res, y = getlfcs1, by = "primerid", sort = FALSE, all.x = TRUE)
          mast_res$contrast <- paste0(category_column,".",ident2)
          mast_res$cluster <- names(dge_outs[[i]])[j]
        }

        mast_res2 <- merge(x = mast_res2, y = getlfcs2, by = "primerid", sort = FALSE, all.x = TRUE)
        mast_res2$contrast <- paste0(category_column,".",ident2)
        mast_res2$cluster <- names(dge_outs[[i]])[j]
        # plotting block: balanced, annotated per-gene plots
        if (nrow(fcHurdleSig) != 0) {
          # user requested number to plot (default 100)
          if (!exists("plot_n") || !is.numeric(plot_n) || length(plot_n) != 1) plot_n <- 100L

          # build stats table from mast_res2 (must contain primerid and logFC)
          stats_all <- mast_res2[, c("primerid", "logFC", "fdr")]
          stats_all$primerid <- as.character(stats_all$primerid)

          # significant primer ids
          sig_ids <- unique(as.character(fcHurdleSig$primerid))
          sig_stats <- stats_all[stats_all$primerid %in% sig_ids, , drop = FALSE]

          # split by sign
          pos_candidates <- sig_stats$primerid[!is.na(sig_stats$logFC) & sig_stats$logFC > 0]
          neg_candidates <- sig_stats$primerid[!is.na(sig_stats$logFC) & sig_stats$logFC < 0]

          # desired counts per side
          n_total <- as.integer(plot_n)
          n_pos_want <- ceiling(n_total / 2)
          n_neg_want <- floor(n_total / 2)

          # take available up to desired
          n_pos_take <- min(length(pos_candidates), n_pos_want)
          n_neg_take <- min(length(neg_candidates), n_neg_want)

          # fill shortage from remaining candidates by |logFC|
          if ((n_pos_take + n_neg_take) < n_total) {
            need_more <- n_total - (n_pos_take + n_neg_take)
            pos_remaining <- setdiff(pos_candidates, head(pos_candidates, n_pos_take))
            neg_remaining <- setdiff(neg_candidates, head(neg_candidates, n_neg_take))
            pool <- c(pos_remaining, neg_remaining)
            if (length(pool) > 0) {
              pool_stats <- sig_stats[match(pool, sig_stats$primerid), , drop = FALSE]
              pool_stats$absLFC <- abs(pool_stats$logFC)
              pool_stats <- pool_stats[order(-pool_stats$absLFC), , drop = FALSE]
              to_add <- head(pool_stats$primerid, need_more)
              pos_extra <- intersect(to_add, pos_remaining)
              neg_extra <- intersect(to_add, neg_remaining)
              n_pos_take <- n_pos_take + length(pos_extra)
              n_neg_take <- n_neg_take + length(neg_extra)
            }
          }

          # select top per side by abs(logFC)
          selected_pos <- character(0)
          selected_neg <- character(0)
          if (n_pos_take > 0 && length(pos_candidates) > 0) {
            pos_stats <- sig_stats[match(pos_candidates, sig_stats$primerid), , drop = FALSE]
            pos_stats$absLFC <- abs(pos_stats$logFC)
            pos_stats <- pos_stats[order(-pos_stats$absLFC), , drop = FALSE]
            selected_pos <- head(pos_stats$primerid, n_pos_take)
          }
          if (n_neg_take > 0 && length(neg_candidates) > 0) {
            neg_stats <- sig_stats[match(neg_candidates, sig_stats$primerid), , drop = FALSE]
            neg_stats$absLFC <- abs(neg_stats$logFC)
            neg_stats <- neg_stats[order(-neg_stats$absLFC), , drop = FALSE]
            selected_neg <- head(neg_stats$primerid, n_neg_take)
          }

          # union and order by abs(logFC)
          selected_set <- unique(c(selected_pos, selected_neg))
          if (length(selected_set) == 0) {
            # fallback: take top by abs(logFC) from all significant
            tmp <- sig_stats
            tmp$absLFC <- abs(tmp$logFC)
            tmp <- tmp[order(-tmp$absLFC), , drop = FALSE]
            selected_set <- head(tmp$primerid, n_total)
          }
          sel_stats <- sig_stats[match(selected_set, sig_stats$primerid), , drop = FALSE]
          # order by descending abs(logFC)
          ord <- order(-abs(sel_stats$logFC))
          selected_set <- selected_set[ord]

          # ensure present in my_sca rownames and cap at n_total
          present_ids <- intersect(selected_set, rownames(my_sca))
          if (length(present_ids) < n_total) {
            # pad with other significant genes by abs(logFC)
            remaining_pool <- setdiff(sig_stats$primerid, selected_set)
            if (length(remaining_pool) > 0) {
              rem_stats <- sig_stats[match(remaining_pool, sig_stats$primerid), , drop = FALSE]
              rem_stats$absLFC <- abs(rem_stats$logFC)
              rem_stats <- rem_stats[order(-rem_stats$absLFC), , drop = FALSE]
              add_ids <- intersect(rem_stats$primerid, rownames(my_sca))
              add_ids <- head(add_ids, n_total - length(present_ids))
              present_ids <- unique(c(present_ids, add_ids))
            }
          }
          present_ids <- unique(head(present_ids, n_total))

          if (length(present_ids) == 0) {
            warning("No significant genes (primerid) found in my_sca rownames after balanced selection; skipping plots.")
            gene_plots <- NA
          } else {
            # optional cap to avoid creating thousands of plots interactively
            max_plot_genes <- 200L
            if (length(present_ids) > max_plot_genes) {
              present_ids <- head(present_ids, max_plot_genes)
              message("Limiting plotting to first ", length(present_ids), " selected genes (set max_plot_genes to change).")
            }

            # build a long table of expression for present genes
            flat_dat <- as(my_sca[present_ids, ], "data.table")
            flat_dat <- data.table::as.data.table(flat_dat)
            if (!"primerid" %in% colnames(flat_dat)) flat_dat[, primerid := rownames(my_sca)[.I]]
            flat_dat[, category := ifelse(category == "Group1", ident1, ident2)]

            # stats lookup for annotations
            stats_dt <- mast_res2[match(present_ids, mast_res2$primerid), c("primerid", "logFC", "fdr")]

            # create list of ggplots with annotations
            gene_plots <- lapply(seq_along(present_ids), function(ii) {
              gid <- present_ids[ii]
              dat_g <- flat_dat[primerid == gid]
              if (nrow(dat_g) == 0) return(NULL)
              stat_row <- stats_dt[ii, , drop = FALSE]
              ann_logFC <- if (is.na(stat_row$logFC)) NA_real_ else stat_row$logFC
              ann_fdr   <- if (is.na(stat_row$fdr)) NA_real_ else stat_row$fdr

              y_max <- suppressWarnings(max(dat_g$logcounts, na.rm = TRUE))
              if (!is.finite(y_max)) y_max <- 0

              title_txt <- gid
              if (!is.null(rowData(seu_as_sce)$gene_symbol)) {
                sym <- rowData(seu_as_sce)$gene_symbol[rownames(seu_as_sce) == gid]
                if (!is.na(sym) && nzchar(sym)) title_txt <- paste0(sym, " (", gid, ")")
              }

              ggplot(dat_g, aes(x = category, y = logcounts, color = category)) +
                ggrastr::geom_jitter_rast(width = 0.2, alpha = 0.6, size = 1) +
                geom_violin(alpha = 0.3, trim = TRUE) +
                ggtitle(title_txt) +
                annotate("text",
                         x = 1.5,
                         y = y_max * 0.95,
                         label = paste0("logFC = ", ifelse(is.na(ann_logFC), "NA", round(ann_logFC, 2)),
                                        "\nFDR = ", ifelse(is.na(ann_fdr), "NA", signif(ann_fdr, 3))),
                         hjust = 0.5, vjust = 0, size = 6) +
                theme_bw(base_size = 14) +
                theme(legend.position = "none")
            })
            names(gene_plots) <- present_ids
            gene_plots <- gene_plots[!vapply(gene_plots, is.null, logical(1))]
          }

          dge_outs[[i]][[j]] <- list(
            filtered_res = mast_res,
            raw_res      = mast_res2,
            gene_plots   = gene_plots,
            zlmfit       = zlmCond,
            sca          = my_sca,
            formula      = fmla,
            fitArgsD     = list(nAGQ = 0),
            zlm_res      = zlmCond
          )

        } else {
          dge_outs[[i]][[j]] <- list(
            filtered_res = "no dge",
            raw_res      = mast_res2,
            gene_plots   = NA,
            zlmfit       = zlmCond,
            sca          = my_sca,
            formula      = fmla,
            fitArgsD     = list(nAGQ = 0),
            zlm_res      = zlmCond
          )
        }
      } else if(tolower(dge_method) == "wilcox") {
        suppressPackageStartupMessages({
          require(SeuratWrappers)
        })
        subs2 <- Seurat::NormalizeData(object = subs2, assay = assay, normalization.method = "LogNormalize")
        wilcox_res <- SeuratWrappers::RunPresto(object = subs2, assay = assay, ident.1 = ident1, ident.2 = ident2,
                                                test.use = "wilcox", only.pos = wilcox_only_positive, min.pct = 0.01)
        # colnames(wilcox_res)[which(colnames(wilcox_res)=="pct.1")] <- gsub(" ","",paste0("pct.",ident1))
        # colnames(wilcox_res)[which(colnames(wilcox_res)=="pct.2")] <- gsub(" ","",paste0("pct.",ident2))
        wilcox_res$gene <- row.names(wilcox_res)
        wilcox_res$cluster <- ident1
        wilcox_res <- wilcox_res[wilcox_res$p_val_adj<0.05,]
        wilcox_res <- wilcox_res[which(abs(wilcox_res$avg_log2FC)>fc_threshold),]

        dge_outs[[i]][[j]] <- wilcox_res
      } else {
        stop("no supported 'dge_method' requested: supported algorithms are: 'wilcox', 'mast', or 'pseudobulk'")
      }
      print(paste0("[",names(dge_outs[[i]])[j],"] in [",names(dge_outs)[i],"] testing took ",round(as.numeric(difftime(Sys.time(), start_j, units = "mins")),2)," mins"))
    }
    print(paste0("[",names(dge_outs)[i],"] testing took ",round(as.numeric(difftime(Sys.time(), start_i, units = "mins")),2)," mins"))
  }
  print(paste0("total test time: ",round(as.numeric(difftime(Sys.time(), start_test, units = "hours")),3)," hours"))
  return(dge_outs)
}
