# sadie fluor
rm(list = ls()); gc()

start_analysis <- Sys.time()

library(FCSimple)
library(flowCore)
library(ggplot2)
library(sccomp)

csv_inpath <- 'D:/fluorescence/sadie/downloads/BIS-76 SADIE 23&me/Flow Export'
result_outpath <- 'D:/fluorescence/sadie/result/limited_panel_multi_res_5NOV2025'
drop_markers <- c('CD56')
withhold_markers <- c('CD275', 'CD274(PD-L1)', 'CD162(PSGL-1)', 'HLA-DR', 'CD11c', 'CD86', 'CX3CR1', 'CD11b', 'CCR2', 'CD40', 'CD62P', 'CD63')
n_knn <- 20

if(!dir.exists(result_outpath)) {
  dir.create(result_outpath, recursive = TRUE)
}

csv_fil_pattern <- paste0('(',paste0(list.files(path = csv_inpath, pattern = '\\csv$'),collapse='|'),')')
csv_fil <- list.files(path = csv_inpath, pattern = '\\csv$', full.names = T);# csv_fil
csv_data <- lapply(X = csv_fil, FUN = data.table::fread, data.table = FALSE, nThread = 20)
ncell <- sapply(X = csv_data, FUN = nrow)

ds_size <- NA

if(!is.na(ds_size)) {
  data_list <- lapply(X = csv_data, FUN = function(arg1, downsample_to = ds_size) {
    if(nrow(arg1)>downsample_to) {
      set.seed(123); return(arg1[sample(1:nrow(arg1),downsample_to,replace=F),])
    } else {
      return(arg1)
    }
  })
} else {
  data_list <- csv_data
}
csv <- do.call(rbind, data_list)

ff <- read.FCS(filename = 'D:/fluorescence/sadie/downloads/BIS-76 SADIE 23&me/Unmixed/BIS-76/SADIE011_SD1.fcs', which.lines = 1:2)

fcs_obj <- list(data = as.matrix(csv), 
                source = rep(gsub(pattern = '_NK cells\\-\\.csv|export_', replacement = '', 
                                  x = stringr::str_extract(string = csv_fil, pattern = csv_fil_pattern)), 
                             sapply(data_list,nrow)), 
                pid = rep(stringr::str_extract(string = csv_fil, pattern = csv_fil_pattern), 
                          sapply(data_list,nrow)), 
                run_date = rep(ff@description$`$DATE`, nrow(csv)))
colnames(fcs_obj$data) <- gsub(pattern = '^.+\\:\\: ', replacement = '', x = colnames(fcs_obj$data))

fcs_obj <- FCSimple::fcs_update(fcs_join_obj = fcs_obj, instrument_type = 'flow')

if(length(drop_markers)!=0) {
  print('dropping:')
  print(drop_markers)
  fcs_obj$data <- fcs_obj$data[,!colnames(fcs_obj$data) %in% drop_markers]
}
head(fcs_obj$data, 3)

capture_full_exprs <- fcs_obj$data
if(length(withhold_markers)!=0) {
  if(sum(withhold_markers %in% colnames(fcs_obj$data))!=length(withhold_markers)) {
    stop('withhold markers incorrectly specified; are they all in the data, as written?')
  } else {
    print('withholding from umap and clustering:')
    print(withhold_markers)
    fcs_obj$data <- fcs_obj$data[,!colnames(fcs_obj$data) %in% withhold_markers]
    fcs_obj$misc <- list('markers_dropped' = drop_markers, 
                         'markers_withheld' = withhold_markers)
  }
}
head(fcs_obj$data, 3)
print(dim(fcs_obj$data))

fcs_obj <- FCSimple::fcs_reduce_dimensions(fcs_join_obj = fcs_obj, use_rep = 'data', algorithm = 'umap', 
                                           language = 'python', umap_nn = n_knn, umap_min_dist = 0.2, nthread = 26)
fcs_obj <- FCSimple::fcs_cluster(fcs_join_obj = fcs_obj, use_rep = 'data', language = 'python', algorithm = 'leiden', 
                                 leiden_louvain_resolution = 0.25, adjacency_knn = n_knn, 
                                 search_method = 'RANN', search_only = FALSE, num_cores = 28)
which_leiden <- which(names(fcs_obj)=='leiden')
names(fcs_obj)[which_leiden] <- paste0('leiden_res',gsub('\\.','p',as.character(fcs_obj[[which_leiden]]$settings$resolution_parameter)))
gc()

fcs_obj <- FCSimple::fcs_cluster(fcs_join_obj = fcs_obj, use_rep = 'data', language = 'python', algorithm = 'leiden', 
                                 leiden_louvain_resolution = 0.5, adjacency_knn = n_knn, 
                                 search_method = 'RANN', search_only = FALSE, num_cores = 28)
which_leiden <- which(names(fcs_obj)=='leiden')
names(fcs_obj)[which_leiden] <- paste0('leiden_res',gsub('\\.','p',as.character(fcs_obj[[which_leiden]]$settings$resolution_parameter)))
gc()

fcs_obj$data <- capture_full_exprs

# prepare testing
sd1 <- fcs_obj$metadata$patient_ID[intersect(grep(pattern = 'SADIE', x = fcs_obj$metadata$patient_ID), 
                                             grep(pattern = 'SD1', x = fcs_obj$metadata$patient_ID))]
sd3 <- fcs_obj$metadata$patient_ID[intersect(grep(pattern = 'SADIE', x = fcs_obj$metadata$patient_ID), 
                                             grep(pattern = 'SD3', x = fcs_obj$metadata$patient_ID))]
if(mean(gsub(pattern = '_(SD1|SD3)', replacement = '', x = sd1)==gsub(pattern = '_(SD1|SD3)', replacement = '', x = sd3))!=1) {
  stop('pid order mismatch')
}
test_groups <- list('V1_sadie' = sd1, 
                    'V3_sadie' = sd3)
test_colors <- list(c('V1_sadie' = 'turquoise3'), 
                    c('V3_sadie' = 'magenta3'))
test_comparisons <- list(c('V1_sadie', 'V3_sadie'))

# res 0.5
fcs_obj <- FCSimple::fcs_cluster_heatmap(fcs_join_obj = fcs_obj, algorithm = 'leiden_res0p5', override_correction = TRUE)
FCSimple::fcs_plot_heatmap(fcs_join_obj = fcs_obj, algorithm = 'leiden_res0p5', outdir = result_outpath, 
                           add_timestamp = FALSE, append_file_string = 'res0p5')
fcs_obj <- FCSimple::fcs_calculate_abundance(fcs_join_obj = fcs_obj, report_algorithm = 'leiden_res0p5', report_as = 'frequency')
FCSimple::fcs_report_abundance(fcs_join_obj = fcs_obj, report_algorithm = 'leiden_res0p5', outdir = result_outpath, 
                               add_timestamp = F, append_file_string = 'res0p5')
umap_p <- FCSimple::fcs_plot_reduction(fcs_join_obj = fcs_obj, algorithm = 'leiden_res0p5', point_size = 1.25, 
                                       reduction = 'umap', point_alpha = 0.5, return_plot = T, 
                                       annotate_text_size = 6) + theme_bw(base_size = 18) + theme(legend.position = 'none')
ggsave(filename = paste0('sadie_flow_umap_res_0p5.pdf'), plot = umap_p, device = 'pdf', path = result_outpath, 
       width = 8, height = 8, units = 'in', dpi = 300, limitsize = F)
FCSimple::fcs_plot_distribution(fcs_join_obj = fcs_obj, override_correction = TRUE, separate_by = 'cluster',
                                plot_algorithm = 'leiden_res0p5', outdir = result_outpath, trim_quantile = c(0.01,0.99),
                                add_timestamp = FALSE)
fcs_obj <- FCSimple::fcs_test_clusters(fcs_join_obj = fcs_obj, compare_list = test_groups, color_list = test_colors, 
                                       comparisons = test_comparisons, denominator_cell_type = 'Mono', 
                                       x_order = c('V1_sadie','V3_sadie'), algorithm = 'leiden_res0p5', dot_size = 2, 
                                       paired_test = T, paired_line_stroke = 0.2, heatmap_fontsize = 15, 
                                       hm_feature_text_size = 18, test_method = 'wilcox', boxplot_text_scaling = 1.5, 
                                       relative_heights = c(0.6,0.4), p_text_size = 8)
pdf(file = file.path(result_outpath, paste0('sadie_sd1_sd3_paired_test_res_0p5.pdf')), width = 10, height = 7)
lapply(X = fcs_obj$leiden_res0p5$cluster_test_results, FUN = function(x) x)
dev.off()

# res 0.25
fcs_obj <- FCSimple::fcs_cluster_heatmap(fcs_join_obj = fcs_obj, algorithm = 'leiden_res0p25', override_correction = TRUE)
FCSimple::fcs_plot_heatmap(fcs_join_obj = fcs_obj, algorithm = 'leiden_res0p25', outdir = result_outpath, 
                           add_timestamp = FALSE, append_file_string = 'res0p25')
fcs_obj <- FCSimple::fcs_calculate_abundance(fcs_join_obj = fcs_obj, report_algorithm = 'leiden_res0p25', report_as = 'frequency')
FCSimple::fcs_report_abundance(fcs_join_obj = fcs_obj, report_algorithm = 'leiden_res0p25', outdir = result_outpath, 
                               add_timestamp = F, append_file_string = 'res0p25')
umap_p <- FCSimple::fcs_plot_reduction(fcs_join_obj = fcs_obj, algorithm = 'leiden_res0p25', point_size = 1.25, 
                                       reduction = 'umap', point_alpha = 0.5, return_plot = T, 
                                       annotate_text_size = 6) + theme_bw(base_size = 18) + theme(legend.position = 'none')
ggsave(filename = paste0('sadie_flow_umap_res_0p25.pdf'), plot = umap_p, device = 'pdf', path = result_outpath, 
       width = 8, height = 8, units = 'in', dpi = 300, limitsize = F)
FCSimple::fcs_plot_distribution(fcs_join_obj = fcs_obj, override_correction = TRUE, separate_by = 'cluster',
                                plot_algorithm = 'leiden_res0p25', outdir = result_outpath, trim_quantile = c(0.01,0.99),
                                add_timestamp = FALSE)
fcs_obj <- FCSimple::fcs_test_clusters(fcs_join_obj = fcs_obj, compare_list = test_groups, color_list = test_colors, 
                                       comparisons = test_comparisons, denominator_cell_type = 'Mono', 
                                       x_order = c('V1_sadie','V3_sadie'), algorithm = 'leiden_res0p25', dot_size = 2, 
                                       paired_test = T, paired_line_stroke = 0.2, heatmap_fontsize = 15, 
                                       hm_feature_text_size = 18, test_method = 'wilcox', boxplot_text_scaling = 1.5, 
                                       relative_heights = c(0.6,0.4), p_text_size = 8)
pdf(file = file.path(result_outpath, paste0('sadie_sd1_sd3_paired_test_res_0p25.pdf')), width = 10, height = 7)
lapply(X = fcs_obj$leiden_res0p25$cluster_test_results, FUN = function(x) x)
dev.off()

FCSimple::fcs_project_parameters(fcs_join_obj = fcs_obj, override_correction = TRUE, reduction = 'umap', sample_size = 150000, 
                                 outdir = result_outpath, point_size = 0.75, point_alpha = 0.3, trim_quantile = 0.03)
FCSimple::fcs_plot_distribution(fcs_join_obj = fcs_obj, override_correction = TRUE, separate_by = 'none', 
                                outdir = result_outpath, trim_quantile = c(0.01,0.99), add_timestamp = FALSE)

# save(fcs_obj, file = file.path(result_outpath, paste0('sadie_fcs_obj_multi_res.RData')))

fcs_outpath <- file.path(result_outpath, 'fcs')
if(!dir.exists(fcs_outpath)) {
  dir.create(path = fcs_outpath, showWarnings = FALSE, recursive = TRUE)
}
unique_src <- unique(fcs_obj$metadata$patient_ID)
for(i in 1:length(unique_src)) {
  print(paste0('writing: ',i,' of ',length(unique_src)))
  row_index <- fcs_obj$source==unique_src[i]
  tmp_exprs <- fcs_obj$data[row_index,]
  tmp_coord <- fcs_obj$umap$coordinates[row_index,]
  tmp_leid_0p25 <- fcs_obj$leiden_res0p25$clusters[row_index]
  tmp_leid_0p5 <- fcs_obj$leiden_res0p5$clusters[row_index]
  
  ff_data <- do.call(cbind, list(as.matrix(tmp_exprs), as.matrix(tmp_coord), as.matrix(data.frame(leiden_0p25 = tmp_leid_0p25, 
                                                                                                  leiden_0p5 = tmp_leid_0p5))))
  new_ff <- new('flowFrame', exprs = ff_data)
  new_ff@description$DATE <- fcs_obj$run_date[row_index][1]
  new_ff@description$SAMPLE <- unique_src[i]
  new_ff@description$WRITE_TIME <- Sys.time()
  
  write.FCS(x = new_ff, filename = file.path(fcs_outpath, paste0(unique_src[i],'_R_multi_res.fcs')))
}

# sccomp
load(file.path(result_outpath, paste0('sadie_fcs_obj_multi_res.RData')))

logit_inv <- function(z) 1/(1+exp(-z))
simplify_sccomp_result <- function(sccomp_res, fdr_cutoff = 0.05) {
  sccomp_res_df <- as.data.frame(sccomp_res)
  n <- nrow(sccomp_res_df)
  
  out <- data.frame(
    celltype = character(n),
    oddsRatio_credibleInterval = character(n),
    bayesian_FDR = numeric(n),
    Interpretation_short = character(n),
    Interpretation_long = character(n),
    stringsAsFactors = FALSE
  )
  compare_str <- sccomp_res_df$factor[1]
  enriched_str <- gsub(pattern = compare_str, replacement = '', x = sccomp_res_df$parameter[1])
  
  for (i in seq_len(n)) {
    row <- sccomp_res_df[i, ]
    ct <- row$celltype
    delta <- row$c_effect
    lo <- row$c_lower
    hi <- row$c_upper
    fdr <- row$c_FDR
    
    # odds ratio and interval
    OR <- exp(delta)
    OR_lo <- exp(lo)
    OR_hi <- exp(hi)
    
    # format OR string
    OR_CI_str <- sprintf("%.2f (%.2f-%.2f)", OR, OR_lo, OR_hi)
    
    # short interpretation
    if (!is.na(fdr) && fdr < fdr_cutoff) {
      if (OR > 1) {
        short <- paste0("Enriched in ",enriched_str)
      } else {
        short <- paste0("Depleted in ",enriched_str)
      }
      long <- sprintf(
        paste0("For cell type %s, the odds of abundance were significantly %s in ",enriched_str," compared to the reference (OR = %.2f, 95%% credible interval: %.2f-%.2f, Bayesian FDR = %.3f)."),
        ct, ifelse(OR > 1, "higher", "lower"), OR, OR_lo, OR_hi, fdr
      )
    } else {
      short <- "No evidence of difference"
      long <- sprintf(
        paste0("For cell type %s, there was no strong evidence of a difference between ",enriched_str," and the reference (OR = %.2f, 95%% credible interval: %.2f-%.2f, Bayesian FDR = %.3f)."),
        ct, OR, OR_lo, OR_hi, fdr
      )
    }
    
    out[i, ] <- list(ct, OR_CI_str, fdr, short, long)
  }
  
  return(out)
}
interval_plots <- function(x, contrasts_string = 'comparison') {
  x$parameter <- gsub(pattern = contrasts_string, replacement = '', x = x$parameter)
  
  attr_x <- attributes(x)
  res_x <- as.data.frame(x)
  
  model_str <- paste0(as.character(attr_x$formula_composition),collapse=' ')
  inference_method <- attr_x$inference_method
  noise_model <- attr_x$noise_model
  
  param_split <- split(x = res_x, f = res_x$parameter)
  
  iterate_parameters <- function(x2) {
    factor_order <- x2$celltype[order(x2$c_effect, decreasing = FALSE)]
    x2$celltype <- factor(x2$celltype, levels = factor_order)
    
    param <- x2$parameter[1]
    param_str <- gsub(pattern = '\\) \\- \\(', replacement = ') -\n(', x = param)
    param_str <- gsub(pattern = 'cmv', replacement = 'pp65', x = param_str)
    param_str <- gsub(pattern = 'HIVp', replacement = 'HIV+', x = param_str)
    param_str <- gsub(pattern = 'HIVn', replacement = 'HIV-', x = param_str)
    
    int_pl <- ggplot(data = x2, mapping = aes(x = c_effect, y = celltype, color = c_FDR < 0.05)) + 
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") + 
      geom_errorbar(aes(xmin = c_lower, xmax = c_upper, color = c_FDR < 0.05), linewidth = 0.8, width = 0.5) + 
      geom_point(mapping = aes(fill = c_FDR < 0.05), size = 4, pch = 21, stroke = 0.4, color = 'black') + 
      scale_color_manual(values = c("grey30", "red")) + 
      scale_fill_manual(values = c("grey30", "red")) + 
      labs(x = "Posterior mean and 95% credible interval (log-odds scale)", 
           y = 'Cell Group',
           title = param_str) + 
      theme_bw(base_size = 14) + 
      theme(axis.text.y = element_text(size = 16), 
            axis.text.x = element_text(size = 15), 
            axis.title.x = element_text(size = 16), 
            axis.title.y = element_blank(), 
            plot.title = element_text(size = 18, hjust = 0.5), 
            legend.position = 'bottom', 
            legend.text = element_text(size = 14), 
            legend.title = element_text(size = 15))
    return(int_pl)
  }
  iterate_pl <- lapply(X = param_split, FUN = iterate_parameters)
  return(iterate_pl)
}
tile_plots <- function(plotlist, n_row = 2, n_col = 2, rm_legend = FALSE, common_leg = FALSE) {
  require(ggplot2)
  require(ggpubr)
  iseq <- seq(from = 1, to = length(plotlist), by = n_row * n_col)
  start_seq <- iseq
  end_seq <- iseq[2:length(iseq)]-1
  if(max(end_seq)<length(plotlist)) {
    end_seq <- append(end_seq, length(plotlist))
  }
  outplots <- vector("list", length = length(start_seq))
  for(i in 1:length(start_seq)) {
    outplots[[i]] <- ggpubr::ggarrange(plotlist = plotlist[start_seq[i]:end_seq[i]], nrow = n_row,
                                       ncol = n_col, legend = ifelse(rm_legend, "none","bottom"), common.legend = common_leg)
    print(start_seq[i]:end_seq[i])
  }
  return(outplots)
}
add_caption <- function(plot_obj) {
  combined_plot <- plot_obj + 
    patchwork::plot_annotation(caption = paste0("Bayesian FDR: Stephens' method (doi: 10.1093/biostatistics/kxw041)\n", 
                                                "FDR-significant populations may cross fold change thresholds because Bayesian FDR considers posterior probabilities rather than p-values.\n", 
                                                "The method sorts null hypothesis probabilities in ascending order and calculates cumulative averages for robust false discovery control."))
  return(combined_plot)
}
get_iqr <- function(df) {
  v3_iqr <- round(quantile(x = df$frequency[df$visit=='V3'], probs = c(0.25,0.5,0.75)),2)
  v3_mean <- round(mean(df$frequency[df$visit=='V3']),2)
  v1_iqr <- round(quantile(x = df$frequency[df$visit=='V1'], probs = c(0.25,0.5,0.75)),2)
  v1_mean <- round(mean(df$frequency[df$visit=='V1']),2)
  report_df <- data.frame('V3' = paste0('median=',v3_iqr[2], ' (IQR=',v3_iqr[1],'-',v3_iqr[3],'; mean=',v3_mean,')'), 
                          'V1' = paste0('median=',v1_iqr[2], ' (IQR=',v1_iqr[1],'-',v1_iqr[3],'; mean=',v1_mean,')'))
  return(list(res = report_df, 
              should_drop = c('V3' = v3_mean, 
                              'V1' = v1_mean)))
}

# res 0p5
# counts <- FCSimple::fcs_calculate_abundance(fcs_join_obj = fcs_obj, report_algorithm = 'leiden_res0p5', 
#                                             report_as = 'count', return_abundance = T)
# write.csv(x = counts, file = file.path(result_outpath, paste0('counts_table_res_0p5.csv')), row.names = TRUE)
counts <- read.csv(file = file.path(result_outpath, paste0('counts_table_res_0p5.csv')), row.names = 1)

leiden_counts <- as.data.frame(counts)
leiden_counts$source <- row.names(leiden_counts)
leiden_counts$pid <- stringr::str_extract(string = row.names(leiden_counts), pattern = 'SADIE[0-9]+|VR[0-9]+')
leiden_counts$visit <- gsub('SD','V',stringr::str_extract(string = row.names(leiden_counts), pattern = 'V[0-9]|SD[0-9]'))
ct_table <- leiden_counts
ct_table$sample_id <- paste0(ct_table$pid,'_',ct_table$visit)
sadie_table <- ct_table[grep(pattern = 'SADIE', x = ct_table$source),]

numcol <- c()
for(i in 1:ncol(sadie_table)) {
  getclass <- class(sadie_table[,i])
  if(getclass %in% c('integer','numeric')) {
    numcol <- c(numcol,i)
    sadie_table[,i] <- as.integer(sadie_table[,i])
  }
}
sadie_table_long <- reshape2::melt(sadie_table[,c(numcol, which(colnames(sadie_table)=='sample_id'))])

colnames(sadie_table_long) <- c('sccomp_sample','celltype','count')
sadie_table_long$celltype <- as.factor(sadie_table_long$celltype)
sadie_table_long$visit <- gsub('_','',stringr::str_extract(string = sadie_table_long$sccomp_sample, pattern = '_.+$'))
sadie_table_long$sccomp_sample <- gsub(pattern = '_.+$', replacement = '', x = sadie_table_long$sccomp_sample)

sccomp_table <- tibble::as_tibble(sadie_table_long)
sccomp_table$celltype <- factor(sccomp_table$celltype)
sccomp_table$pid <- factor(sccomp_table$sccomp_sample)
sccomp_table$sccomp_sample <- factor(paste0(sccomp_table$sccomp_sample,'_',sccomp_table$visit))
sccomp_table$visit <- factor(sccomp_table$visit, levels = c('V1','V3'))

test_table <- as.data.frame(sccomp_table)
test_table_spl <- split(x = test_table, f = test_table$sccomp_sample)
test_table_spl <- lapply(X = test_table_spl, FUN = function(x){
  x$frequency <- (x$count/sum(x$count)) * 100
  return(x)
})
freq_long <- do.call(rbind, test_table_spl)
freq_spl <- split(x = freq_long, f = freq_long$celltype)
freq_spl <- freq_spl[order(as.numeric(names(freq_spl)), decreasing = F)]

get_iqr <- function(df) {
  v3_iqr <- round(quantile(x = df$frequency[df$visit=='V3'], probs = c(0.25,0.5,0.75)),2)
  v3_mean <- round(mean(df$frequency[df$visit=='V3']),2)
  v1_iqr <- round(quantile(x = df$frequency[df$visit=='V1'], probs = c(0.25,0.5,0.75)),2)
  v1_mean <- round(mean(df$frequency[df$visit=='V1']),2)
  report_df <- data.frame('V3' = paste0('median=',v3_iqr[2], ' (IQR=',v3_iqr[1],'-',v3_iqr[3],'; mean=',v3_mean,')'), 
                          'V1' = paste0('median=',v1_iqr[2], ' (IQR=',v1_iqr[1],'-',v1_iqr[3],'; mean=',v1_mean,')'))
  return(list(res = report_df, 
              should_drop = c('V3' = v3_mean, 
                              'V1' = v1_mean)))
}
summary_stats <- lapply(X = freq_spl, FUN = get_iqr)
summary_stat <- lapply(X = summary_stats, FUN = function(x) return(x[[1]]))
summary_stat <- do.call(rbind, summary_stat)
summary_stat <- cbind(data.frame(celltype = names(freq_spl)), 
                      summary_stat)
write.csv(x = summary_stat, file = file.path(result_outpath,
                                             paste0('sadie_visit_sccomp_summary_stats_res_0p5.csv')), row.names = FALSE)
fmla_visit <- as.formula(object = " ~ visit + (1 | pid)")
sccomp::clear_stan_model_cache()
sccomp_est_visit <- sccomp_estimate(.data = sccomp_table, 
                                    formula_composition = fmla_visit, 
                                    sample = 'sccomp_sample', 
                                    cell_group = 'celltype', 
                                    abundance = 'count', 
                                    inference_method = 'hmc',
                                    bimodal_mean_variability_association = F, # recommended to be TRUE for RNAseq
                                    mcmc_seed = 123,
                                    chains = 10, # can increase to 6–8 for more cross‐chain checks
                                    adapt_delta = 0.99, # tighter step‐size target
                                    max_treedepth = 20, # default is 10
                                    cores = 30,
                                    verbose = F)
clean_visit <- sccomp_remove_outliers(.estimate = sccomp_est_visit, cores = 28)
sccomp_res_visit <- sccomp_test(clean_visit)
visit_contrast <- subset(sccomp_res_visit, parameter == "visitV3")

report_res_visit <- simplify_sccomp_result(visit_contrast)
report_res_visit <- report_res_visit[order(as.numeric(report_res_visit$celltype),decreasing=FALSE),]
intervals_1d_visit <- interval_plots(x = visit_contrast, contrasts_string = 'visit')[[1]]
itervals_w_caption_visit <- add_caption(intervals_1d_visit)
write.csv(x = report_res_visit, file = file.path(result_outpath, paste0('sadie_visit_sccomp_res_0p5.csv')), row.names = FALSE)
ggsave(filename = paste0('sadie_visit_sccomp_intervals_res_0p5.pdf'), plot = itervals_w_caption_visit,
       device = 'pdf', path = result_outpath, width = 8, height = 8, units = 'in', dpi = 300, limitsize = F)

# res 0p5, AHA
counts <- read.csv(file = file.path(result_outpath, paste0('counts_table_res_0p5.csv')), row.names = 1, check.names = F)

leiden_counts <- as.data.frame(counts)
leiden_counts$source <- row.names(leiden_counts)
leiden_counts$pid <- stringr::str_extract(string = row.names(leiden_counts), pattern = 'SADIE[0-9]+|VR[0-9]+')
leiden_counts$visit <- gsub('SD','V',stringr::str_extract(string = row.names(leiden_counts), pattern = 'V[0-9]|SD[0-9]'))
ct_table <- leiden_counts
ct_table$sample_id <- paste0(ct_table$pid,'_',ct_table$visit)
aha_table <- ct_table[grep(pattern = '^VR', x = ct_table$source),]

numcol <- c()
for(i in 1:ncol(aha_table)) {
  getclass <- class(aha_table[,i])
  if(getclass %in% c('integer','numeric')) {
    numcol <- c(numcol,i)
    aha_table[,i] <- as.integer(aha_table[,i])
  }
}
aha_table_long <- reshape2::melt(aha_table[,c(numcol, which(colnames(aha_table)=='sample_id'))])

colnames(aha_table_long) <- c('sccomp_sample','celltype','count')
aha_table_long$celltype <- as.factor(aha_table_long$celltype)
aha_table_long$visit <- gsub('_','',stringr::str_extract(string = aha_table_long$sccomp_sample, pattern = '_.+$'))
aha_table_long$sccomp_sample <- gsub(pattern = '_.+$', replacement = '', x = aha_table_long$sccomp_sample)

sccomp_table <- tibble::as_tibble(aha_table_long)
sccomp_table$celltype <- factor(sccomp_table$celltype)
sccomp_table$pid <- factor(sccomp_table$sccomp_sample)
sccomp_table$sccomp_sample <- factor(paste0(sccomp_table$sccomp_sample,'_',sccomp_table$visit))
sccomp_table$visit <- factor(sccomp_table$visit, levels = c('V1','V2'))

test_table <- as.data.frame(sccomp_table)
test_table_spl <- split(x = test_table, f = test_table$sccomp_sample)
test_table_spl <- lapply(X = test_table_spl, FUN = function(x){
  x$frequency <- (x$count/sum(x$count)) * 100
  return(x)
})
freq_long <- do.call(rbind, test_table_spl)
freq_spl <- split(x = freq_long, f = freq_long$celltype)
freq_spl <- freq_spl[order(as.numeric(names(freq_spl)), decreasing = F)]

get_iqr <- function(df, reference_group = 'V1', test_group = 'V2') {
  test_iqr <- round(quantile(x = df$frequency[df$visit==test_group], probs = c(0.25,0.5,0.75)),2)
  test_mean <- round(mean(df$frequency[df$visit==test_group]),2)
  reference_iqr <- round(quantile(x = df$frequency[df$visit==reference_group], probs = c(0.25,0.5,0.75)),2)
  reference_mean <- round(mean(df$frequency[df$visit==reference_group]),2)
  report_df <- data.frame('x1' = paste0('median=',test_iqr[2], ' (IQR=',test_iqr[1],'-',test_iqr[3],'; mean=',test_mean,')'), 
                          'x2' = paste0('median=',reference_iqr[2], ' (IQR=',reference_iqr[1],'-',reference_iqr[3],'; mean=',reference_mean,')'))
  colnames(report_df) <- c(test_group, reference_group)
  return(list(res = report_df, 
              should_drop = c('test' = test_mean, 
                              'reference' = reference_mean)))
}
summary_stats <- lapply(X = freq_spl, FUN = get_iqr)
summary_stat <- lapply(X = summary_stats, FUN = function(x) return(x[[1]]))
summary_stat <- do.call(rbind, summary_stat)
summary_stat <- cbind(data.frame(celltype = names(freq_spl)), 
                      summary_stat)
write.csv(x = summary_stat, file = file.path(result_outpath, 
                                             paste0('aha_visit_sccomp_summary_stats_res_0p5.csv')), row.names = FALSE)
fmla_visit <- as.formula(object = " ~ visit + (1 | pid)")
sccomp::clear_stan_model_cache()
sccomp_est_visit <- sccomp_estimate(.data = sccomp_table, 
                                    formula_composition = fmla_visit, 
                                    sample = 'sccomp_sample', 
                                    cell_group = 'celltype', 
                                    abundance = 'count', 
                                    inference_method = 'hmc', 
                                    bimodal_mean_variability_association = FALSE, # recommended to be TRUE for RNAseq
                                    mcmc_seed = 123, 
                                    chains = 10, # can increase to 6–8 for more cross‐chain checks
                                    adapt_delta = 0.99, # tighter step‐size target  
                                    max_treedepth = 20, # default is 10
                                    cores = 30,
                                    verbose = F)
clean_visit <- sccomp_remove_outliers(.estimate = sccomp_est_visit, cores = 28)
sccomp_res_visit <- sccomp_test(clean_visit)
visit_contrast <- subset(sccomp_res_visit, parameter == "visitV2")

report_res_visit <- simplify_sccomp_result(visit_contrast)
report_res_visit <- report_res_visit[order(as.numeric(report_res_visit$celltype),decreasing=FALSE),]
intervals_1d_visit <- interval_plots(x = visit_contrast, contrasts_string = 'visit')[[1]]
itervals_w_caption_visit <- add_caption(intervals_1d_visit)
write.csv(x = report_res_visit, file = file.path(result_outpath, paste0('aha_visit_sccomp_res_0p5.csv')), row.names = FALSE)
ggsave(filename = paste0('aha_visit_sccomp_intervals_res_0p5.pdf'), plot = itervals_w_caption_visit, 
       device = 'pdf', path = result_outpath, width = 8, height = 8, units = 'in', dpi = 300, limitsize = F)

# res 0p25
counts <- FCSimple::fcs_calculate_abundance(fcs_join_obj = fcs_obj, report_algorithm = 'leiden_res0p25', 
                                            report_as = 'count', return_abundance = T)
write.csv(x = counts, file = file.path(result_outpath, paste0('counts_table_res_0p25.csv')), row.names = TRUE)

leiden_counts <- as.data.frame(counts)
leiden_counts$source <- row.names(leiden_counts)
leiden_counts$pid <- stringr::str_extract(string = row.names(leiden_counts), pattern = 'SADIE[0-9]+|VR[0-9]+')
leiden_counts$visit <- gsub('SD','V',stringr::str_extract(string = row.names(leiden_counts), pattern = 'V[0-9]|SD[0-9]'))
ct_table <- leiden_counts
ct_table$sample_id <- paste0(ct_table$pid,'_',ct_table$visit)
sadie_table <- ct_table[grep(pattern = 'SADIE', x = ct_table$source),]

numcol <- c()
for(i in 1:ncol(sadie_table)) {
  getclass <- class(sadie_table[,i])
  if(getclass %in% c('integer','numeric')) {
    numcol <- c(numcol,i)
    sadie_table[,i] <- as.integer(sadie_table[,i])
  }
}
sadie_table_long <- reshape2::melt(sadie_table[,c(numcol, which(colnames(sadie_table)=='sample_id'))])

colnames(sadie_table_long) <- c('sccomp_sample','celltype','count')
sadie_table_long$celltype <- as.factor(sadie_table_long$celltype)
sadie_table_long$visit <- gsub('_','',stringr::str_extract(string = sadie_table_long$sccomp_sample, pattern = '_.+$'))
sadie_table_long$sccomp_sample <- gsub(pattern = '_.+$', replacement = '', x = sadie_table_long$sccomp_sample)

sccomp_table <- tibble::as_tibble(sadie_table_long)
sccomp_table$celltype <- factor(sccomp_table$celltype)
sccomp_table$pid <- factor(sccomp_table$sccomp_sample)
sccomp_table$sccomp_sample <- factor(paste0(sccomp_table$sccomp_sample,'_',sccomp_table$visit))
sccomp_table$visit <- factor(sccomp_table$visit, levels = c('V1','V3'))

test_table <- as.data.frame(sccomp_table)
test_table_spl <- split(x = test_table, f = test_table$sccomp_sample)
test_table_spl <- lapply(X = test_table_spl, FUN = function(x){
  x$frequency <- (x$count/sum(x$count)) * 100
  return(x)
})
freq_long <- do.call(rbind, test_table_spl)
freq_spl <- split(x = freq_long, f = freq_long$celltype)
freq_spl <- freq_spl[order(as.numeric(names(freq_spl)), decreasing = F)]

summary_stats <- lapply(X = freq_spl, FUN = get_iqr)
summary_stat <- lapply(X = summary_stats, FUN = function(x) return(x[[1]]))
summary_stat <- do.call(rbind, summary_stat)
summary_stat <- cbind(data.frame(celltype = names(freq_spl)), 
                      summary_stat)
write.csv(x = summary_stat, file = file.path(result_outpath, 
                                             paste0('sadie_visit_sccomp_summary_stats_res_0p25.csv')), row.names = FALSE)
fmla_visit <- as.formula(object = " ~ visit + (1 | pid)")
sccomp::clear_stan_model_cache()
sccomp_est_visit <- sccomp_estimate(.data = sccomp_table, 
                                    formula_composition = fmla_visit, 
                                    sample = 'sccomp_sample', 
                                    cell_group = 'celltype', 
                                    abundance = 'count', 
                                    inference_method = 'hmc', 
                                    bimodal_mean_variability_association = FALSE, # recommended to be TRUE for RNAseq
                                    mcmc_seed = 123, 
                                    chains = 10, # can increase to 6–8 for more cross‐chain checks
                                    adapt_delta = 0.99, # tighter step‐size target  
                                    max_treedepth = 20, # default is 10
                                    cores = 30,
                                    verbose = F)
clean_visit <- sccomp_remove_outliers(.estimate = sccomp_est_visit, cores = 28)
sccomp_res_visit <- sccomp_test(clean_visit)
visit_contrast <- subset(sccomp_res_visit, parameter == "visitV3")

report_res_visit <- simplify_sccomp_result(visit_contrast)
report_res_visit <- report_res_visit[order(as.numeric(report_res_visit$celltype),decreasing=FALSE),]
intervals_1d_visit <- interval_plots(x = visit_contrast, contrasts_string = 'visit')[[1]]
itervals_w_caption_visit <- add_caption(intervals_1d_visit)
write.csv(x = report_res_visit, file = file.path(result_outpath, paste0('sadie_visit_sccomp_res_0p25.csv')), row.names = FALSE)
ggsave(filename = paste0('sadie_visit_sccomp_intervals_res_0p25.pdf'), plot = itervals_w_caption_visit, 
       device = 'pdf', path = result_outpath, width = 8, height = 8, units = 'in', dpi = 300, limitsize = F)

print(paste0('total runtime: ',round(as.numeric(difftime(Sys.time(), start_analysis, units = 'mins')),2), ' mins'))

