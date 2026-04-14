##############
#INTRODUCTION#
##############

#This code performs MR analyses between VATadjBMI and PCOS - replicating the results from Liu et al.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

###################
#Loading functions#
###################

compute_neff_gsem <- function(eaf, se, Ncases, Ncontrols) {
  maf <- pmin(eaf, 1 - eaf)  # Compute MAF from EAF
  Neff <- 4 / (2 * maf * (1 - maf) * se^2)  # Compute Neff using GSEM formula
  
  # Calculate total effective N
  v <- Ncases / (Ncases + Ncontrols)
  TotalNeff <- 4 * v * (1 - v) * (Ncases + Ncontrols)
  
  # Apply upper and lower limits to Neff
  Neff <- ifelse(Neff > 1.1 * TotalNeff, 1.1 * TotalNeff, Neff)
  Neff <- ifelse(Neff < 0.5 * TotalNeff, 0.5 * TotalNeff, Neff)
  
  return(Neff)
}


mr_plots <- function(dat)
{
  require(TwoSampleMR)
  require(ggplot2)
  require(ggrepel)
  require(dplyr)
  
  temp <- subset(dat, outcome == outcome[1] & exposure == exposure[1])
  exposure_name <- temp$exposure[1]
  outcome_name <- temp$outcome[1]
  
  if(! "labels" %in% names(dat)) dat$labels <- NA
  
  exposure_units <- temp$units.exposure[1]
  outcome_units <- temp$units.outcome[1]
  
  mrs <- mr_singlesnp(temp, all_method=c("mr_egger_regression", "mr_ivw", "mr_weighted_median", "mr_weighted_mode"))
  mrl <- mr_leaveoneout(temp)
  
  mrs$index <- 1:nrow(mrs)
  mrl$index <- 1:nrow(mrl)
  
  mrs <- dplyr::arrange(merge(mrs, select(temp, SNP, labels), all.x=TRUE), index)
  mrl <- dplyr::arrange(merge(mrl, select(temp, SNP, labels), all.x=TRUE), index)
  
  mrres <- mr(temp, method_list=c("mr_egger_regression", "mr_ivw", "mr_weighted_median", "mr_weighted_mode"))
  
  p <- gridExtra::grid.arrange(
    mr_forest_plot(mrs)[[1]] +
      ggplot2::labs(
        title="a)",
        x=paste0("MR estimate (", outcome_units, " per ", exposure_units, ")")
      ),
    mr_scatter_plot(mrres, temp)[[1]] +
      ggplot2::labs(
        title="b)",
        x=paste0("SNP effect on ", exposure_name),
        y=paste0("SNP effect on ", outcome_name)
      ) +
      geom_label_repel(ggplot2::aes(label=labels), point.padding = unit(0.5, "lines")),
    mr_leaveoneout_plot(mrl)[[1]] +
      ggplot2::labs(
        title="c)",
        x=paste0("MR estimate (", outcome_units, " per ", exposure_units, ")"),
        y="Excluded variant"
      ),
    mr_funnel_plot(mrs)[[1]] +
      ggplot2::labs(title="d)") +
      ggplot2::theme(legend.position="none") +
      ggrepel::geom_label_repel(ggplot2::aes(label=labels), point.padding = unit(0.5, "lines")),
    ncol=2
  )
  
  return(p)
}


remove_outlier <- function(mr_tmp_df, rucker){
  #OUT.0: let's make vectors for the models:
  ivw_vect <- c("A", "B")
  egger_vect <- c("C", "D")
  #OUT.1: Let's check the model that fits the data better: IVW (A, B) or Egger (C, D).
  rucker_model <- rucker[[1]]$res
  #Let's format the data:
  df_4_radial <- RadialMR::format_radial(BXG = mr_tmp_df$beta.exposure, BYG = mr_tmp_df$beta.outcome,
                                         seBXG = mr_tmp_df$se.exposure, seBYG = mr_tmp_df$se.outcome,
                                         RSID = mr_tmp_df$SNP)
  if(rucker_model%in%ivw_vect){
    radial_output <- tryCatch(RadialMR::ivw_radial(r_input = df_4_radial, alpha = 0.05, weights = 3),  error = function(e){
      return(NA)
    })
  } else {
    radial_output <- tryCatch(RadialMR::egger_radial(r_input = df_4_radial, alpha = 0.05, weights = 3), error = function(e){
      return(NA)
    })
  }
  if(length(radial_output) == 1){
    return(NA)
  }
  #IF WE DON'T HAVE OUTLIERS RETURN ORIGINAL
  if(length(radial_output$outlier) ==  1){
    mr_wo_outliers <- mr_tmp_df
  } else {
    outliers <- radial_output$outliers$SNP
    mr_wo_outliers <- mr_tmp_df[which(!(mr_tmp_df$SNP%in%outliers)),]
  }
  return(mr_wo_outliers)
}


compute_strict_het <- function(dat_1_filt_1){
  
  dat_1_filt_1$mr <- dat_1_filt_1$beta.outcome/dat_1_filt_1$beta.exposure
  dat_1_filt_1$mr_se <- ((dat_1_filt_1$mr*((dat_1_filt_1$se.exposure/dat_1_filt_1$beta.exposure)^2+(dat_1_filt_1$se.outcome/dat_1_filt_1$beta.outcome)^2)^0.5)^2)^0.5
  het_Q_Isq <- meta::metagen(dat_1_filt_1$mr, dat_1_filt_1$mr_se)
  
  return(het_Q_Isq)
  
}

#####################################################
#STEP 1: let's obtain the data for each associations#
#####################################################

#Set pathway to find data:

path_2_input <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/PCOS_2026"

setwd(path_2_input)

#Let's load the data for all WHRadjBMI mathches

liu_pcos_doctor <- fread("output/2_replicating_liu_et_al/4_mr_fat_depots/1_expo_outcome_df/vatadj_PCOS.txt")
liu_pcos_broad <- fread("output/2_replicating_liu_et_al/4_mr_fat_depots/1_expo_outcome_df/vatadj_PCOS (Broad).txt")
liu_pcos_consortium <- fread("output/2_replicating_liu_et_al/4_mr_fat_depots/1_expo_outcome_df/vatadj_PCOS (Consortium).txt")
liu_pcos_meta_analysis <- fread("output/2_replicating_liu_et_al/4_mr_fat_depots/1_expo_outcome_df/vatadj_PCOS (Day et al).txt")
liu_pcos_venkatesh <- fread("output/2_replicating_liu_et_al/4_mr_fat_depots/1_expo_outcome_df/vatadj_PCOS (Venkatesh et al).txt")
liu_pcos_adj_age <- fread("output/2_replicating_liu_et_al/4_mr_fat_depots/1_expo_outcome_df/vatadj_PCOS (adj age).txt")
liu_pcos_adj_age_bmi <- fread("output/2_replicating_liu_et_al/4_mr_fat_depots/1_expo_outcome_df/vatadj_PCOS (adj age+BMI).txt")

##############################################################
#STEP 2: let's first rerun the analyses for the meta-analysis#
##############################################################

mr_df <- liu_pcos_meta_analysis[which(((as.numeric(liu_pcos_meta_analysis$beta.exposure)^2)/(as.numeric(liu_pcos_meta_analysis$se.exposure)^2)) > 10),] #8/8

#Let's make all numeric or it will go crazy

mr_df$beta.exposure <- as.numeric(mr_df$beta.exposure)
mr_df$se.exposure <- as.numeric(mr_df$se.exposure)
mr_df$eaf.exposure <- as.numeric(mr_df$eaf.exposure)

mr_df$beta.outcome <- as.numeric(mr_df$beta.outcome)
mr_df$se.outcome <- as.numeric(mr_df$se.outcome)
mr_df$eaf.outcome <- as.numeric(mr_df$eaf.outcome)

#STEP 2: let's compute MR steiger:

mr_df$r.exposure <- TwoSampleMR::get_r_from_pn(
  p = as.numeric(mr_df$pval.exposure),n = as.numeric(mr_df$samplesize.exposure))

mr_df$r.outcome <- TwoSampleMR::get_r_from_lor(
  lor = as.numeric(mr_df$beta.outcome),
  af = as.numeric(mr_df$eaf.outcome),
  ncase = as.numeric(mr_df$ncase.outcome),
  ncontrol = as.numeric(mr_df$ncontrol.outcome),
  prevalence = mr_df$prevalence.outcome[1] 
)

#Adding additional columns to make steiger run...

mr_df$id.exposure <- "vatadj"
mr_df$id.outcome <- "PCOS (meta-analysis)"
mr_df$exposure <- "vatadj"
mr_df$outcome <- "PCOS (meta-analysis)"
mr_df$units.exposure <- "SD"
mr_df$units.outcome <- "LOR"

mr_steiger_df <- TwoSampleMR::steiger_filtering(mr_df)
mr_steiger_df <- mr_steiger_df[which(mr_steiger_df$steiger_dir == TRUE),] #1
mr_steiger_df$mr_keep <- TRUE

#Let's save the data:

dir.create("output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/")
dir.create("output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/vatadj/")
dir.create("output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/vatadj/meta_analysis/")
saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/vatadj/meta_analysis/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure           id.outcome              outcome exposure                    method nsnp           b        se      pval
# 1      vatadj PCOS (meta-analysis) PCOS (meta-analysis)   vatadj                  MR Egger    8 -1.18657297 1.4780834 0.4527207
# 2      vatadj PCOS (meta-analysis) PCOS (meta-analysis)   vatadj           Weighted median    8 -0.02185381 0.2305092 0.9244684
# 3      vatadj PCOS (meta-analysis) PCOS (meta-analysis)   vatadj Inverse variance weighted    8  0.02115538 0.1790141 0.9059273
# 4      vatadj PCOS (meta-analysis) PCOS (meta-analysis)   vatadj               Simple mode    8 -0.16507726 0.3563556 0.6572521
# 5      vatadj PCOS (meta-analysis) PCOS (meta-analysis)   vatadj             Weighted mode    8 -0.09718880 0.3177546 0.7686061

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method   Estimate        SE      CI_low   CI_upp         P
# 1  Egger fixed effects 0.09396366 0.0864737 -0.07552168 0.263449 0.2772068
# 2 Egger random effects 0.09396366 0.0864737 -0.07552168 0.263449 0.1386034
# 
# [[1]]$Q
# Method         Q df         P
# 1   Q_ivw 4.1207396  7 0.7657703
# 2 Q_egger 3.4431635  6 0.7515156
# 3  Q_diff 0.6775762  1 0.4104226
# 
# [[1]]$res
# [1] "A"

strict_het <- compute_strict_het(mr_steiger_df)

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0 [0.0000; 0.3372]; tau = 0 [0.0000; 0.5807]
# I^2 = 0.0% [0.0%; 67.6%]; H = 1.00 [1.00; 1.76]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 3.99    7  0.7815

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/vatadj/meta_analysis/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/vatadj/meta_analysis/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/vatadj/meta_analysis//sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/vatadj/meta_analysis//strict_het_test_before_outlier_extraction.RDS")

#Egger is quite biased, these two variants have no effect...
#Let's check if RadialMR removes this...

#mr_post <- remove_outlier(mr_steiger_df, rucker) #no significant outliers!

#saveRDS(mr_post, "output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/vatadj/meta_analysis/instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

#mr_res <- TwoSampleMR::mr(mr_post)

#Let's check the plot:

#rucker <- TwoSampleMR::mr_rucker(mr_post)

#Here there is heterogeneity and pleitropy

#strict_het <- compute_strict_het(mr_post) #they disagree, but still

#That is quite good:

# mr_post$labels <- NA
# mr_post$units.exposure <- "SD"
# mr_post$units.outcome <- "LOR"
# 
# plotio=mr_plots(mr_post)
# ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/meta_analysis/plots_after_outlier_extraction.svg", height=16, width=16, units="in")
# 
# saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/meta_analysis/mr_res_after_outlier_extraction.RDS")
# saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/meta_analysis/sensitivity_test_after_outlier_extraction.RDS")
# saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/meta_analysis/strict_het_test_after_outlier_extraction.RDS")

###############################################################
#Let's first analyse the data for the whole set for vatadj - PCOS#
###############################################################

#STEP 1: let's use only those that have F>10

mr_df <- liu_pcos_doctor[which(((liu_pcos_doctor$beta.exposure^2)/(liu_pcos_doctor$se.exposure^2)) > 10),] #9/9

#STEP 2: let's compute MR steiger:

mr_df$r.exposure <- TwoSampleMR::get_r_from_pn(
  p = mr_df$pval.exposure,n = mr_df$samplesize.exposure)

mr_df$r.outcome <- TwoSampleMR::get_r_from_lor(
  lor = mr_df$beta.outcome,
  af = mr_df$eaf.outcome,
  ncase = mr_df$ncase.outcome,
  ncontrol = mr_df$ncontrol.outcome,
  prevalence = mr_df$prevalence.outcome[1] #whole index population unadjusted period prevalence in % from https://risteys.finngen.fi/endpoints/E4_PCOS
)

mr_steiger_df <- TwoSampleMR::steiger_filtering(mr_df)
mr_steiger_df <- mr_steiger_df[which(mr_steiger_df$steiger_dir == TRUE),]
mr_steiger_df$mr_keep <- TRUE

#Let's save the data:

dir.create("output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/")
dir.create("output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/doctor")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/doctor/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp          b        se      pval
# 1    exposure    outcome outcome exposure                  MR Egger    9 1.33764855 0.8685785 0.1674458
# 2    exposure    outcome outcome exposure           Weighted median    9 0.10378676 0.2225947 0.6410303
# 3    exposure    outcome outcome exposure Inverse variance weighted    9 0.03643231 0.2167546 0.8665196
# 4    exposure    outcome outcome exposure               Simple mode    9 0.12981062 0.3888485 0.7470876
# 5    exposure    outcome outcome exposure             Weighted mode    9 0.13803264 0.3815663 0.7269055

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method   Estimate         SE     CI_low       CI_upp          P
# 1  Egger fixed effects -0.1116509 0.05367328 -0.2168485 -0.006453171 0.03750769
# 2 Egger random effects -0.1116509 0.07252004 -0.2537875  0.030485811 0.93816938
# 
# [[1]]$Q
# Method        Q df          P
# 1   Q_ivw 17.10624  8 0.02902162
# 2 Q_egger 12.77903  7 0.07767881
# 3  Q_diff  4.32721  1 0.03750769
# 
# [[1]]$res
# [1] "C"

strict_het <- compute_strict_het(mr_steiger_df)

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0.1930 [0.0000; 1.2216]; tau = 0.4393 [0.0000; 1.1052]
# I^2 = 47.8% [0.0%; 75.7%]; H = 1.38 [1.00; 2.03]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 15.33    8  0.0530

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/doctor/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/doctor/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/doctor/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/doctor/strict_het_test_before_outlier_extraction.RDS")

mr_post <- remove_outlier(mr_steiger_df, rucker) #no significant outliers!

saveRDS(mr_post, "output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/vatadj/doctor/instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure                    method nsnp           b        se      pval
# 1    exposure    outcome outcome exposure                  MR Egger    8 -0.84036809 1.5019992 0.5960726
# 2    exposure    outcome outcome exposure           Weighted median    8  0.02312718 0.2200889 0.9163114
# 3    exposure    outcome outcome exposure Inverse variance weighted    8 -0.14591351 0.1827307 0.4245706
# 4    exposure    outcome outcome exposure               Simple mode    8  0.14031391 0.3739609 0.7186147
# 5    exposure    outcome outcome exposure             Weighted mode    8  0.14708202 0.3736991 0.7056013

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method   Estimate         SE     CI_low    CI_upp         P
# 1  Egger fixed effects 0.05489361 0.09802065 -0.1372233 0.2470105 0.5754652
# 2 Egger random effects 0.05489361 0.11773314 -0.1758591 0.2856463 0.3205166
# 
# [[1]]$Q
# Method         Q df         P
# 1   Q_ivw 8.9695497  7 0.2548497
# 2 Q_egger 8.6559262  6 0.1938731
# 3  Q_diff 0.3136234  1 0.5754652
# 
# [[1]]$res
# [1] "A"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_post) #they disagree, but still

#Quantifying heterogeneity (with 95%-CIs):
#  tau^2 = 0.0382 [0.0000; 0.9085]; tau = 0.1954 [0.0000; 0.9532]
#I^2 = 16.3% [0.0%; 59.3%]; H = 1.09 [1.00; 1.57]

#Test of heterogeneity:
#  Q d.f. p-value
#8.36    7  0.3016

#That is quite good:

mr_post$labels <- NA
mr_post$units.exposure <- "SD"
mr_post$units.outcome <- "LOR"
 
plotio=mr_plots(mr_post)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/doctor/plots_after_outlier_extraction.svg", height=16, width=16, units="in")
 
saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/doctor/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/doctor/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/doctor/strict_het_test_after_outlier_extraction.RDS")


###############################################################
#Let's first analyse the data for the whole set for bmi - PCOS# broad definition - anovulation
###############################################################

#STEP 1: let's use only those that have F>10

mr_df <- liu_pcos_broad[which(((liu_pcos_broad$beta.exposure^2)/(liu_pcos_broad$se.exposure^2)) > 10),] #2

#STEP 2: let's compute MR steiger:

mr_df$r.exposure <- TwoSampleMR::get_r_from_pn(
  p = mr_df$pval.exposure,n = mr_df$samplesize.exposure)

mr_df$r.outcome <- TwoSampleMR::get_r_from_lor(
  lor = mr_df$beta.outcome,
  af = mr_df$eaf.outcome,
  ncase = mr_df$ncase.outcome,
  ncontrol = mr_df$ncontrol.outcome,
  prevalence = mr_df$prevalence.outcome[1] #whole index population unadjusted period prevalence in % from https://risteys.finngen.fi/endpoints/E4_PCOS_BROAD
)

mr_steiger_df <- TwoSampleMR::steiger_filtering(mr_df)
mr_steiger_df <- mr_steiger_df[which(mr_steiger_df$steiger_dir == TRUE),]
mr_steiger_df$mr_keep <- TRUE

#Let's save the data:

dir.create("output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/broad")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/broad/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp          b        se      pval
# 1    exposure    outcome outcome exposure                  MR Egger    9  0.9993988 0.6071359 0.1437408
# 2    exposure    outcome outcome exposure           Weighted median    9 -0.1143044 0.1531785 0.4555364
# 3    exposure    outcome outcome exposure Inverse variance weighted    9 -0.1355168 0.1604056 0.3982011
# 4    exposure    outcome outcome exposure               Simple mode    9 -0.1823151 0.2592564 0.5018658
# 5    exposure    outcome outcome exposure             Weighted mode    9 -0.1410282 0.2454956 0.5814436

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method    Estimate         SE     CI_low       CI_upp           P
# 1  Egger fixed effects -0.09716229 0.03702180 -0.1697237 -0.024600895 0.008678609
# 2 Egger random effects -0.09716229 0.05060186 -0.1963401  0.002015531 0.972579438
# 
# [[1]]$Q
# Method         Q df           P
# 1   Q_ivw 19.965033  8 0.010469153
# 2 Q_egger 13.077235  7 0.070248786
# 3  Q_diff  6.887798  1 0.008678609
# 
# [[1]]$res
# [1] "C"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_steiger_df) #they disagree, but still

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0.1188 [0.0010; 0.6550]; tau = 0.3446 [0.0324; 0.8093]
# I^2 = 54.9% [4.6%; 78.7%]; H = 1.49 [1.02; 2.17]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 17.74    8  0.0232

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/broad/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/broad/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/broad/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/broad/strict_het_test_before_outlier_extraction.RDS")

#A single variant seems to be introducing quite a bit of heterogeneity.
#I do not think this might affect results, but let's see.

mr_post <- remove_outlier(mr_steiger_df, rucker)
 
saveRDS(mr_post, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/broad/instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure                    method nsnp           b        se       pval
# 1    exposure    outcome outcome exposure                  MR Egger    7 -0.07160627 1.1080893 0.95097951
# 2    exposure    outcome outcome exposure           Weighted median    7 -0.10796796 0.1588881 0.49680712
# 3    exposure    outcome outcome exposure Inverse variance weighted    7 -0.22029098 0.1219284 0.07080522
# 4    exposure    outcome outcome exposure               Simple mode    7 -0.03096208 0.2445295 0.90337856
# 5    exposure    outcome outcome exposure             Weighted mode    7 -0.04262127 0.2107952 0.84644751

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method    Estimate         SE     CI_low    CI_upp         P
# 1  Egger fixed effects -0.01135331 0.08178786 -0.1716546 0.1489479 0.8895970
# 2 Egger random effects -0.01135331 0.08406889 -0.1761253 0.1534187 0.5537129
# 
# [[1]]$Q
# Method          Q df         P
# 1   Q_ivw 5.30205500  6 0.5056961
# 2 Q_egger 5.28278563  5 0.3823521
# 3  Q_diff 0.01926937  1 0.8895970
# 
# [[1]]$res
# [1] "A"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_post) #they disagree, but still

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0 [0.0000; 0.3688]; tau = 0 [0.0000; 0.6073]
# I^2 = 0.0% [0.0%; 70.8%]; H = 1.00 [1.00; 1.85]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 5.02    6  0.5417

#That is quite good:

mr_post$labels <- NA
mr_post$units.exposure <- "SD"
mr_post$units.outcome <- "LOR"
 
plotio=mr_plots(mr_post)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/broad/plots_after_outlier_extraction.svg", height=16, width=16, units="in")
 
saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/broad/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/broad/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/broad/strict_het_test_after_outlier_extraction.RDS")

################################################################
#Let's first analyse the data for the whole set for vatadj - PCOS# consortium definition
################################################################

#STEP 1: let's use only those that have F>10

mr_df <- liu_pcos_consortium[which(((liu_pcos_consortium$beta.exposure^2)/(liu_pcos_consortium$se.exposure^2)) > 10),] #3

#STEP 2: let's compute MR steiger:

mr_df$r.exposure <- TwoSampleMR::get_r_from_pn(
  p = mr_df$pval.exposure,n = mr_df$samplesize.exposure)

mr_df$r.outcome <- TwoSampleMR::get_r_from_lor(
  lor = mr_df$beta.outcome,
  af = mr_df$eaf.outcome,
  ncase = mr_df$ncase.outcome,
  ncontrol = mr_df$ncontrol.outcome,
  prevalence = mr_df$prevalence.outcome[1] #whole index population unadjusted period prevalence in % from https://risteys.finngen.fi/endpoints/E4_PCOS_CONCORTIUM
)

mr_steiger_df <- TwoSampleMR::steiger_filtering(mr_df)
mr_steiger_df <- mr_steiger_df[which(mr_steiger_df$steiger_dir == TRUE),]
mr_steiger_df$mr_keep <- TRUE

#Let's save the data:

dir.create("output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/consortium")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/consortium/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp           b         se      pval
# 1    exposure    outcome outcome exposure                  MR Egger    9  0.09892835 0.23926164 0.6916328
# 2    exposure    outcome outcome exposure           Weighted median    9 -0.04910894 0.05572106 0.3781364
# 3    exposure    outcome outcome exposure Inverse variance weighted    9 -0.04483902 0.05142286 0.3832266
# 4    exposure    outcome outcome exposure               Simple mode    9 -0.04695949 0.10055106 0.6529403
# 5    exposure    outcome outcome exposure             Weighted mode    9 -0.04154984 0.10146364 0.6929145

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method    Estimate         SE      CI_low     CI_upp         P
# 1  Egger fixed effects -0.01224799 0.01370053 -0.03910054 0.01460456 0.3713332
# 2 Egger random effects -0.01224799 0.01986657 -0.05118574 0.02668977 0.7312218
# 
# [[1]]$Q
# Method          Q df          P
# 1   Q_ivw 15.5178731  8 0.04982424
# 2 Q_egger 14.7186749  7 0.03977950
# 3  Q_diff  0.7991981  1 0.37133322
# 
# [[1]]$res
# [1] "B"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_steiger_df) #they disagree, but still

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0.0104 [0.0000; 0.0623]; tau = 0.1022 [0.0000; 0.2496]
# I^2 = 43.8% [0.0%; 74.1%]; H = 1.33 [1.00; 1.96]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 14.25    8  0.0755

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/consortium/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/consortium/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/consortium/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/consortium/strict_het_test_before_outlier_extraction.RDS")

#There is pleiotropy.

#STEP 5: removing outliers!

mr_post <- remove_outlier(mr_steiger_df, rucker) #but not significant outliers!!
  
saveRDS(mr_post, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/consortium/instruments_after_outlier_removal_df.RDS")
  
#STEP 6: let's run the results for real now:
  
mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure                    method nsnp           b         se      pval
# 1    exposure    outcome outcome exposure                  MR Egger    7 -0.50943115 0.38945807 0.2477607
# 2    exposure    outcome outcome exposure           Weighted median    7 -0.04699529 0.05539038 0.3961943
# 3    exposure    outcome outcome exposure Inverse variance weighted    7 -0.04261239 0.04412442 0.3341781
# 4    exposure    outcome outcome exposure               Simple mode    7 -0.05463733 0.08722031 0.5541031
# 5    exposure    outcome outcome exposure             Weighted mode    7 -0.04611156 0.08421237 0.6037482

# #Let's check the plot:
  
rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method   Estimate         SE       CI_low     CI_upp          P
# 1  Egger fixed effects 0.03565323 0.01919406 -0.001966427 0.07327289 0.06323803
# 2 Egger random effects 0.03565323 0.01919406 -0.001966427 0.07327289 0.03161902
# 
# [[1]]$Q
# Method        Q df         P
# 1   Q_ivw 3.564485  6 0.7353714
# 2 Q_egger 2.109073  5 0.8338561
# 3  Q_diff 1.455412  1 0.2276615
# 
# [[1]]$res
# [1] "A"

# #Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_post) #they disagree, but still
  
#That is quite good:
  
mr_post$labels <- NA
mr_post$units.exposure <- "SD"
mr_post$units.outcome <- "LOR"
  
plotio=mr_plots(mr_post)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/consortium/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/consortium/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/consortium/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/consortium/strict_het_test_after_outlier_extraction.RDS")

###############################################################
#Let's first analyse the data for the whole set for vatadj - PCOS# age adjusted
###############################################################

#STEP 1: let's use only those that have F>10

mr_df <- liu_pcos_adj_age[which(((liu_pcos_adj_age$beta.exposure^2)/(liu_pcos_adj_age$se.exposure^2)) > 10),] #9
liu_pcos_doctor_match=liu_pcos_doctor[which(liu_pcos_doctor$SNP%in%mr_df$SNP),]
mr_df=mr_df[which(mr_df$SNP%in%liu_pcos_doctor_match$SNP),]
mr_df=mr_df[order(match(mr_df$SNP, liu_pcos_doctor_match$SNP)),]
print(length(mr_df$SNP==liu_pcos_doctor_match$SNP))

mr_df$eaf.outcome=liu_pcos_doctor_match$eaf.outcome

#STEP 2: let's compute MR steiger:

mr_df$r.exposure <- TwoSampleMR::get_r_from_pn(
  p = mr_df$pval.exposure,n = mr_df$samplesize.exposure)

mr_df$r.outcome <- TwoSampleMR::get_r_from_lor(
  lor = mr_df$beta.outcome,
  af = mr_df$eaf.outcome,
  ncase = mr_df$ncase.outcome,
  ncontrol = mr_df$ncontrol.outcome,
  prevalence = mr_df$prevalence.outcome[1] 
)

mr_steiger_df <- TwoSampleMR::steiger_filtering(mr_df)
mr_steiger_df <- mr_steiger_df[which(mr_steiger_df$steiger_dir == TRUE),]
mr_steiger_df$mr_keep <- TRUE

#Let's save the data:

dir.create("output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/age_adjusted")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/age_adjusted/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp         b        se      pval
# 1    exposure    outcome outcome exposure                  MR Egger    9 0.4099381 0.6907147 0.5715184
# 2    exposure    outcome outcome exposure           Weighted median    9 0.1951359 0.1734474 0.2605707
# 3    exposure    outcome outcome exposure Inverse variance weighted    9 0.1645174 0.1405333 0.2417335
# 4    exposure    outcome outcome exposure               Simple mode    9 0.2616675 0.2946362 0.4003979
# 5    exposure    outcome outcome exposure             Weighted mode    9 0.2358937 0.3069777 0.4643038

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method    Estimate         SE     CI_low     CI_upp         P
# 1  Egger fixed effects -0.02064505 0.04697519 -0.1127147 0.07142464 0.6603078
# 2 Egger random effects -0.02064505 0.05673870 -0.1318509 0.09056076 0.6420194
# 
# [[1]]$Q
# Method          Q df         P
# 1   Q_ivw 10.4053587  8 0.2377193
# 2 Q_egger 10.2122088  7 0.1768634
# 3  Q_diff  0.1931499  1 0.6603078
# 
# [[1]]$res
# [1] "A"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_steiger_df) #they disagree, but still

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0.0468 [0.0000; 0.4504]; tau = 0.2164 [0.0000; 0.6711]
# I^2 = 19.3% [0.0%; 60.9%]; H = 1.11 [1.00; 1.60]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 9.92    8  0.2710

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/age_adjusted//plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/age_adjusted//mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/age_adjusted//sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/age_adjusted//strict_het_test_before_outlier_extraction.RDS")

#We have a weird outlier! It does not affect things much, but let's see!

mr_post <- remove_outlier(mr_steiger_df, rucker) #no significant outliers. 

saveRDS(mr_post, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/age_adjusted//instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure                    method nsnp         b        se       pval
# 1    exposure    outcome outcome exposure                  MR Egger    8 0.8367417 0.5963944 0.21018895
# 2    exposure    outcome outcome exposure           Weighted median    8 0.2134985 0.1746961 0.22166463
# 3    exposure    outcome outcome exposure Inverse variance weighted    8 0.2896501 0.1339769 0.03062276
# 4    exposure    outcome outcome exposure               Simple mode    8 0.2133970 0.2866024 0.48077265
# 5    exposure    outcome outcome exposure             Weighted mode    8 0.2053189 0.2742651 0.47848333

# #Let's check the plot:
  
rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method    Estimate         SE     CI_low     CI_upp         P
# 1  Egger fixed effects -0.04515979 0.03846389 -0.1205476 0.03022805 0.2403619
# 2 Egger random effects -0.04515979 0.03846389 -0.1205476 0.03022805 0.8798191
# 
# [[1]]$Q
# Method         Q df         P
# 1   Q_ivw 4.7436336  7 0.6912154
# 2 Q_egger 3.8574122  6 0.6959653
# 3  Q_diff 0.8862214  1 0.3465033
# 
# [[1]]$res
# [1] "A"

# #Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_post) #they disagree, but still

#That is quite good:

mr_post$labels <- NA
mr_post$units.exposure <- "SD"
mr_post$units.outcome <- "LOR"

plotio=mr_plots(mr_post)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/age_adjusted/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/age_adjusted/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/age_adjusted/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/age_adjusted/strict_het_test_after_outlier_extraction.RDS")

###############################################################
#Let's first analyse the data for the whole set for bmi - PCOS# age and BMI adjusted
###############################################################

#STEP 1: let's use only those that have F>10

mr_df <- liu_pcos_adj_age_bmi[which(((liu_pcos_adj_age_bmi$beta.exposure^2)/(liu_pcos_adj_age_bmi$se.exposure^2)) > 10),] #3
liu_pcos_doctor_match=liu_pcos_doctor[which(liu_pcos_doctor$SNP%in%mr_df$SNP),]
mr_df=mr_df[which(mr_df$SNP%in%liu_pcos_doctor_match$SNP),]
mr_df=mr_df[order(match(mr_df$SNP, liu_pcos_doctor_match$SNP)),]
print(length(mr_df$SNP==liu_pcos_doctor_match$SNP))

mr_df$eaf.outcome=liu_pcos_doctor_match$eaf.outcome

#STEP 2: let's compute MR steiger:

mr_df$r.exposure <- TwoSampleMR::get_r_from_pn(
  p = mr_df$pval.exposure,n = mr_df$samplesize.exposure)

mr_df$r.outcome <- TwoSampleMR::get_r_from_lor(
  lor = mr_df$beta.outcome,
  af = mr_df$eaf.outcome,
  ncase = mr_df$ncase.outcome,
  ncontrol = mr_df$ncontrol.outcome,
  prevalence = mr_df$prevalence.outcome[1] 
)

mr_steiger_df <- TwoSampleMR::steiger_filtering(mr_df)
mr_steiger_df <- mr_steiger_df[which(mr_steiger_df$steiger_dir == TRUE),]
mr_steiger_df$mr_keep <- TRUE

#Let's save the data:

dir.create("output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/age_and_bmi_adjusted")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/age_and_bmi_adjusted/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp            b        se      pval
# 1    exposure    outcome outcome exposure                  MR Egger    9  0.448863570 0.8623435 0.6187573
# 2    exposure    outcome outcome exposure           Weighted median    9 -0.116474888 0.1987026 0.5577567
# 3    exposure    outcome outcome exposure Inverse variance weighted    9 -0.003371093 0.1772629 0.9848272
# 4    exposure    outcome outcome exposure               Simple mode    9 -0.292785291 0.3095985 0.3719882
# 5    exposure    outcome outcome exposure             Weighted mode    9 -0.265820364 0.3066402 0.4112528

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method    Estimate         SE     CI_low     CI_upp         P
# 1  Egger fixed effects -0.03802222 0.05521793 -0.1462474 0.07020294 0.4910847
# 2 Egger random effects -0.03802222 0.07080135 -0.1767903 0.10074588 0.7043754
# 
# [[1]]$Q
# Method          Q df         P
# 1   Q_ivw 11.9827051  8 0.1519772
# 2 Q_egger 11.5085562  7 0.1179231
# 3  Q_diff  0.4741489  1 0.4910847
# 
# [[1]]$res
# [1] "A"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_steiger_df) #they disagree, but still

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0.0762 [0.0000; 0.8326]; tau = 0.2760 [0.0000; 0.9125]
# I^2 = 27.9% [0.0%; 66.5%]; H = 1.18 [1.00; 1.73]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 11.10    8  0.1960

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/age_and_bmi_adjusted/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/age_and_bmi_adjusted/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/age_and_bmi_adjusted/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/age_and_bmi_adjusted/strict_het_test_before_outlier_extraction.RDS")

#We have pleiotropy, let's check this out:

mr_post <- remove_outlier(mr_steiger_df, rucker) #no significant outliers. 

saveRDS(mr_post, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/age_and_bmi_adjusted//instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure                    method nsnp          b        se      pval
# 1    exposure    outcome outcome exposure                  MR Egger    8  0.6924648 0.6777545 0.3463363
# 2    exposure    outcome outcome exposure           Weighted median    8 -0.1915558 0.1932666 0.3216134
# 3    exposure    outcome outcome exposure Inverse variance weighted    8 -0.1450069 0.1540412 0.3465250
# 4    exposure    outcome outcome exposure               Simple mode    8 -0.3732993 0.3202855 0.2819820
# 5    exposure    outcome outcome exposure             Weighted mode    8 -0.3732993 0.3100874 0.2677716

# #Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method    Estimate         SE     CI_low      CI_upp          P
# 1  Egger fixed effects -0.07158272 0.04041226 -0.1507893 0.007623863 0.07650887
# 2 Egger random effects -0.07158272 0.04041226 -0.1507893 0.007623863 0.96174556
# 
# [[1]]$Q
# Method        Q df         P
# 1   Q_ivw 4.688883  7 0.6978693
# 2 Q_egger 3.078868  6 0.7988849
# 3  Q_diff 1.610015  1 0.2044897
# 
# [[1]]$res
# [1] "A"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_post) #they disagree, but still

#That is quite good:

mr_post$labels <- NA
mr_post$units.exposure <- "SD"
mr_post$units.outcome <- "LOR"

plotio=mr_plots(mr_post)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/age_and_bmi_adjusted/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/age_and_bmi_adjusted/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/age_and_bmi_adjusted/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/age_and_bmi_adjusted/strict_het_test_after_outlier_extraction.RDS")

###############################################################
#Let's first analyse the data for the whole set for bmi - PCOS# Venkatesh
###############################################################

#STEP 1: let's use only those that have F>10

mr_df <- liu_pcos_venkatesh[which(((liu_pcos_venkatesh$beta.exposure^2)/(liu_pcos_venkatesh$se.exposure^2)) > 10),] #5

#STEP 2: let's compute MR steiger:

mr_df$r.exposure <- TwoSampleMR::get_r_from_pn(
  p = mr_df$pval.exposure,n = mr_df$samplesize.exposure)

mr_df$r.outcome <- TwoSampleMR::get_r_from_lor(
  lor = mr_df$beta.outcome,
  af = mr_df$eaf.outcome,
  ncase = mr_df$ncase.outcome,
  ncontrol = mr_df$ncontrol.outcome,
  prevalence = mr_df$prevalence.outcome[1] 
)

mr_steiger_df <- TwoSampleMR::steiger_filtering(mr_df)
mr_steiger_df <- mr_steiger_df[which(mr_steiger_df$steiger_dir == TRUE),]
mr_steiger_df$mr_keep <- TRUE

#Let's save the data:

dir.create("output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/venkatesh")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/venkatesh/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp         b        se      pval
# 1    exposure    outcome outcome exposure                  MR Egger   10 0.5101747 0.5756814 0.4013580
# 2    exposure    outcome outcome exposure           Weighted median   10 0.2126890 0.1689565 0.2080887
# 3    exposure    outcome outcome exposure Inverse variance weighted   10 0.1564273 0.1345632 0.2450398
# 4    exposure    outcome outcome exposure               Simple mode   10 0.3713493 0.2985195 0.2449336
# 5    exposure    outcome outcome exposure             Weighted mode   10 0.3892955 0.3042031 0.2326454

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method    Estimate         SE     CI_low     CI_upp         P
# 1  Egger fixed effects -0.02963324 0.04012898 -0.1082846 0.04901811 0.4602411
# 2 Egger random effects -0.02963324 0.04679191 -0.1213437 0.06207721 0.7367306
# 
# [[1]]$Q
# Method          Q df         P
# 1   Q_ivw 11.4224610  9 0.2478586
# 2 Q_egger 10.8771527  8 0.2087575
# 3  Q_diff  0.5453083  1 0.4602411
# 
# [[1]]$res
# [1] "A"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_steiger_df) #they disagree, but still

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0.0450 [0.0000; 0.4896]; tau = 0.2120 [0.0000; 0.6997]
# I^2 = 17.8% [0.0%; 58.8%]; H = 1.10 [1.00; 1.56]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 10.94    9  0.2796

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/venkatesh/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/venkatesh/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/venkatesh/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/venkatesh/strict_het_test_before_outlier_extraction.RDS")

#Asymmetry! Do we have outliers?

#mr_post <- remove_outlier(mr_steiger_df, rucker) #no outliers

#saveRDS(mr_post, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/venkatesh/instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

#mr_res <- TwoSampleMR::mr(mr_post)

# #Let's check the plot:
  
#rucker <- TwoSampleMR::mr_rucker(mr_post)

# #Here there is heterogeneity and pleitropy

#strict_het <- compute_strict_het(mr_post) #they disagree, but still

#That is quite good:

# mr_post$labels <- NA
# mr_post$units.exposure <- "SD"
# mr_post$units.outcome <- "LOR"
# 
# plotio=mr_plots(mr_post)
# ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/venkatesh/plots_after_outlier_extraction.svg", height=16, width=16, units="in")
# 
# saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/venkatesh/mr_res_after_outlier_extraction.RDS")
# saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/venkatesh/sensitivity_test_after_outlier_extraction.RDS")
# saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vatadj/venkatesh/strict_het_test_after_outlier_extraction.RDS")