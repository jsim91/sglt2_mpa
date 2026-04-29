# corr platelet mpa
rm(list = ls()); gc()

# install.packages('pspearman')
library(pspearman)

inpath <- 'D:/10x/sadie/MPA_pilot/reports'

correlate_celltype <- 'Platelet' # CD14 Mono, Platelet

counts_df <- read.csv(file = file.path(inpath, 'mm_pbmc_named_cluster_cell_counts.csv'), check.names = FALSE)
pid_coln <- if('pt_day_id' %in% colnames(counts_df)) 'pt_day_id' else colnames(counts_df)[1]
row.names(counts_df) <- counts_df[,pid_coln]
counts_df <- counts_df[,setdiff(colnames(counts_df), pid_coln), drop = FALSE]

if(!('M-platelet' %in% colnames(counts_df) && correlate_celltype %in% colnames(counts_df))) {
  stop("Expected 'M-platelet' and target correlate_celltype columns in mm_pbmc_named_cluster_cell_counts.csv")
}

freq_df <- apply(X = counts_df, MARGIN = 1, FUN = function(arg){
  return((arg/sum(arg))*100)
})
freq_df <- t(freq_df)
freq_df <- as.data.frame(freq_df)

sd_coln <- which(colnames(freq_df) %in% c('M-platelet',correlate_celltype))
sd1 <- freq_df[grep(pattern = 'SD1$', x = row.names(freq_df)),sd_coln]
sd2 <- freq_df[grep(pattern = 'SD2$', x = row.names(freq_df)),sd_coln]
sd3 <- freq_df[grep(pattern = 'SD3$', x = row.names(freq_df)),sd_coln]

sd1_corr <- pspearman::spearman.test(x = sd1$`M-platelet`, y = sd1[,correlate_celltype])
sd2_corr <- pspearman::spearman.test(x = sd2$`M-platelet`, y = sd2[,correlate_celltype])
sd3_corr <- pspearman::spearman.test(x = sd3$`M-platelet`, y = sd3[,correlate_celltype])
sdn_corr <- pspearman::spearman.test(x = freq_df$`M-platelet`, y = freq_df[,correlate_celltype])

print(paste0('SD1: pvalue = ',sd1_corr$p.value,', rho = ',sd1_corr$estimate))
print(paste0('SD2: pvalue = ',sd2_corr$p.value,', rho = ',sd2_corr$estimate))
print(paste0('SD3: pvalue = ',sd3_corr$p.value,', rho = ',sd3_corr$estimate))
print(paste0('all SD: pvalue = ',sdn_corr$p.value,', rho = ',sdn_corr$estimate))

# Platelet result:
# "SD1: pvalue = 0.664583333333333, rho = 0.19047619047619"
# "SD2: pvalue = 0.299206349206349, rho = 0.428571428571429"
# "SD3: pvalue = 0.582142857142857, rho = 0.238095238095238"
# "all SD: pvalue = 0.119103019512185, rho = 0.32695652173913"