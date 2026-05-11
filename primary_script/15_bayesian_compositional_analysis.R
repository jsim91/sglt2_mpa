rm(list=ls())
library(sccomp)

mydf <- readRDS(here::here('intermediate/sccomp/pbmc2data.rds'))
mydf$sample <- paste(mydf$subject_id, mydf$study_day, sep='_')
mydf$time <- factor(mydf$study_day, levels=1:3, labels=c('Baseline', '2 Weeks', '12 Weeks'))
mydf$cell_type <- factor(mydf$cell_type)
mydf$subject_id <- factor(mydf$subject_id)
mydf$drug <- ifelse(mydf$subject_id %in% c("Sadie-016", "Sadie-018", "Sadie-023",  "Sadie-029"), 1, 0)

cleanup_dirs <- c(
  here::here('intermediate/sccomp/pbmc2mod_output'),
  here::here('intermediate/sccomp/pbmc2modc_output')
)

on.exit({
  for (d in cleanup_dirs) {
    if (dir.exists(d)) {
      unlink(d, recursive = TRUE, force = TRUE)
    }
  }
}, add = TRUE)


pbmc2mod <- sccomp_estimate(mydf, formula_composition = ~ time + (1|subject_id), .sample = sample, 
.cell_group = cell_type, .abundance = count, bimodal_mean_variability_association = TRUE, cores = 1, 
verbose = FALSE, mcmc_seed=7, output_directory=here::here('intermediate/sccomp/pbmc2mod_output')) 
saveRDS(pbmc2mod, file=here::here('intermediate/sccomp/pbmc2mod.rds'))

pbmc2modc <- sccomp_estimate(mydf, formula_composition = ~ time + drug + (1|subject_id), .sample = sample, 
.cell_group = cell_type, .abundance = count, bimodal_mean_variability_association = TRUE, cores = 1, 
verbose = FALSE, mcmc_seed=7, output_directory=here::here('intermediate/sccomp/pbmc2modc_output')) 
saveRDS(pbmc2modc, file=here::here('intermediate/sccomp/pbmc2modc.rds'))


### not including drug as a covariate ###

fit1 <- readRDS(here::here('intermediate/sccomp/pbmc2mod.rds'))
rst1 <- sccomp_test(fit1, percent_false_positive=5, test_composition_above_logit_fold_change=0.2)

## Logit fold change
subset(as.data.frame(rst1), parameter %in% c('time2 Weeks'))
subset(as.data.frame(rst1), parameter %in% c('time12 Weeks'))

## Change in cell proportions
pdfun <- function(mod=fit1) {
  newdf <- data.frame(sample=c('a', 'b', 'c'), time=c('Baseline', '2 Weeks', '12 Weeks'))
  newdf$time <- factor(newdf$time, levels=c('Baseline', '2 Weeks', '12 Weeks'))
  
  tmp1 <- sccomp_predict(mod, formula_composition = ~ time, new_data=newdf, summary_instead_of_draws=FALSE, number_of_draws=1000, mcmc_seed=7)
  
  dt1 <- subset(tmp1, time=='Baseline')
  dt2 <- subset(tmp1, time=='2 Weeks')
  dt3 <- subset(tmp1, time=='12 Weeks')
  dt2$proportion <- dt2$proportion - dt1$proportion
  dt3$proportion <- dt3$proportion - dt1$proportion
  
  dt <- rbind(dt2[,c('time', 'cell_type', 'proportion')], dt3[,c('time', 'cell_type', 'proportion')])
  
  tmp2 <- aggregate(dt$proportion, by=list(dt$time, dt$cell_type), FUN=mean)
  tmp3 <- aggregate(dt$proportion, by=list(dt$time, dt$cell_type), FUN=function(x) quantile(x, 0.025))
  tmp4 <- aggregate(dt$proportion, by=list(dt$time, dt$cell_type), FUN=function(x) quantile(x, 0.975))
  
  rstdf <- Reduce(f=function(x,y) merge(x,y, by=c('Group.1', 'Group.2')), x=list(tmp2, tmp3, tmp4))
  colnames(rstdf) <- c('time', 'cell_group', 'm', 'lower', 'upper')
  return(rstdf)
}

pdfun(mod=fit1)


### including drug as a covariate ###

fit2 <- readRDS(here::here('intermediate/sccomp/pbmc2modc.rds'))
rst2 <- sccomp_test(fit2, percent_false_positive=5, test_composition_above_logit_fold_change=0.2)

## Logit fold change
subset(as.data.frame(rst2), parameter %in% c('time2 Weeks'))
subset(as.data.frame(rst2), parameter %in% c('time12 Weeks'))

## Change in cell proportions

pdfun(mod=fit2)
