##############
#INTRODUCTION#
##############

#This code performs MR analyses between asat_gsat_ratioBMI and PCOS - replicating the results from Liu et al.

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

liu_pcos_doctor <- fread("output/2_replicating_liu_et_al/4_mr_fat_depots/1_expo_outcome_df/asat_gsat_ratio_PCOS.txt")
liu_pcos_broad <- fread("output/2_replicating_liu_et_al/4_mr_fat_depots/1_expo_outcome_df/asat_gsat_ratio_PCOS (Broad).txt")
liu_pcos_consortium <- fread("output/2_replicating_liu_et_al/4_mr_fat_depots/1_expo_outcome_df/asat_gsat_ratio_PCOS (Consortium).txt")
liu_pcos_meta_analysis <- fread("output/2_replicating_liu_et_al/4_mr_fat_depots/1_expo_outcome_df/asat_gsat_ratio_PCOS (Day et al).txt")
liu_pcos_venkatesh <- fread("output/2_replicating_liu_et_al/4_mr_fat_depots/1_expo_outcome_df/asat_gsat_ratio_PCOS (Venkatesh et al).txt")
liu_pcos_adj_age <- fread("output/2_replicating_liu_et_al/4_mr_fat_depots/1_expo_outcome_df/asat_gsat_ratio_PCOS (adj age).txt")
liu_pcos_adj_age_bmi <- fread("output/2_replicating_liu_et_al/4_mr_fat_depots/1_expo_outcome_df/asat_gsat_ratio_PCOS (adj age+BMI).txt")

##############################################################
#STEP 2: let's first rerun the analyses for the meta-analysis#
##############################################################

mr_df <- liu_pcos_meta_analysis[which(((as.numeric(liu_pcos_meta_analysis$beta.exposure)^2)/(as.numeric(liu_pcos_meta_analysis$se.exposure)^2)) > 10),] #7/7

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

mr_df$id.exposure <- "asat_gsat_ratio"
mr_df$id.outcome <- "PCOS (meta-analysis)"
mr_df$exposure <- "asat_gsat_ratio"
mr_df$outcome <- "PCOS (meta-analysis)"
mr_df$units.exposure <- "SD"
mr_df$units.outcome <- "LOR"

mr_steiger_df <- TwoSampleMR::steiger_filtering(mr_df)
mr_steiger_df <- mr_steiger_df[which(mr_steiger_df$steiger_dir == TRUE),] #19
mr_steiger_df$mr_keep <- TRUE

#Let's save the data:

dir.create("output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/")
dir.create("output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/asat_gsat_ratio/")
dir.create("output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/asat_gsat_ratio/meta_analysis/")
saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/asat_gsat_ratio/meta_analysis/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure           id.outcome              outcome        exposure                    method nsnp          b        se      pval
# 1 asat_gsat_ratio PCOS (meta-analysis) PCOS (meta-analysis) asat_gsat_ratio                  MR Egger    7 0.18054298 0.6730031 0.7992115
# 2 asat_gsat_ratio PCOS (meta-analysis) PCOS (meta-analysis) asat_gsat_ratio           Weighted median    7 0.04212673 0.2263782 0.8523742
# 3 asat_gsat_ratio PCOS (meta-analysis) PCOS (meta-analysis) asat_gsat_ratio Inverse variance weighted    7 0.05576550 0.1775482 0.7534553
# 4 asat_gsat_ratio PCOS (meta-analysis) PCOS (meta-analysis) asat_gsat_ratio               Simple mode    7 0.10440111 0.3148758 0.7514943
# 5 asat_gsat_ratio PCOS (meta-analysis) PCOS (meta-analysis) asat_gsat_ratio             Weighted mode    7 0.09325777 0.2840505 0.7538294

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method    Estimate         SE     CI_low     CI_upp         P
# 1  Egger fixed effects -0.01007739 0.04666161 -0.1015325 0.08137768 0.8290131
# 2 Egger random effects -0.01007739 0.05205209 -0.1120976 0.09194284 0.5767563
# 
# [[1]]$Q
# Method          Q df         P
# 1   Q_ivw 6.26859874  6 0.3937840
# 2 Q_egger 6.22195678  5 0.2852168
# 3  Q_diff 0.04664197  1 0.8290131
# 
# [[1]]$res
# [1] "A"

strict_het <- compute_strict_het(mr_steiger_df)

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 < 0.0001 [0.0000; 1.1138]; tau = 0.0002 [0.0000; 1.0554]
# I^2 = 0.0% [0.0%; 70.8%]; H = 1.00 [1.00; 1.85]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 5.83    6  0.4422

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/asat_gsat_ratio/meta_analysis/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/asat_gsat_ratio/meta_analysis/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/asat_gsat_ratio/meta_analysis//sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/asat_gsat_ratio/meta_analysis//strict_het_test_before_outlier_extraction.RDS")

#Slight asymmetry, but very tame - results are what they are!

#mr_post <- remove_outlier(mr_steiger_df, rucker) #no significant outliers!

#saveRDS(mr_post, "output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/asat_gsat_ratio/meta_analysis/instruments_after_outlier_removal_df.RDS")

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
# ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/meta_analysis/plots_after_outlier_extraction.svg", height=16, width=16, units="in")
# 
# saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/meta_analysis/mr_res_after_outlier_extraction.RDS")
# saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/meta_analysis/sensitivity_test_after_outlier_extraction.RDS")
# saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/meta_analysis/strict_het_test_after_outlier_extraction.RDS")

###############################################################
#Let's first analyse the data for the whole set for asat_gsat_ratio - PCOS#
###############################################################

#STEP 1: let's use only those that have F>10

mr_df <- liu_pcos_doctor[which(((liu_pcos_doctor$beta.exposure^2)/(liu_pcos_doctor$se.exposure^2)) > 10),] #19/19

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

dir.create("output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/")
dir.create("output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/doctor")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/doctor/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp           b        se      pval
# 1    exposure    outcome outcome exposure                  MR Egger    7  1.44230371 0.8920443 0.1668332
# 2    exposure    outcome outcome exposure           Weighted median    7  0.28097408 0.2449441 0.2513425
# 3    exposure    outcome outcome exposure Inverse variance weighted    7  0.35338134 0.2630925 0.1792128
# 4    exposure    outcome outcome exposure               Simple mode    7 -0.09040908 0.4765145 0.8557766
# 5    exposure    outcome outcome exposure             Weighted mode    7  0.35059104 0.3506885 0.3560417

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method    Estimate         SE     CI_low      CI_upp          P
# 1  Egger fixed effects -0.08644971 0.04745907 -0.1794678 0.006568359 0.06852125
# 2 Egger random effects -0.08644971 0.06796957 -0.2196676 0.046768200 0.89829362
# 
# [[1]]$Q
# Method         Q df          P
# 1   Q_ivw 13.573682  6 0.03477836
# 2 Q_egger 10.255588  5 0.06830883
# 3  Q_diff  3.318094  1 0.06852125
# 
# [[1]]$res
# [1] "B"

strict_het <- compute_strict_het(mr_steiger_df)

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0.2565 [0.0000; 2.5423]; tau = 0.5064 [0.0000; 1.5945]
# I^2 = 51.4% [0.0%; 79.4%]; H = 1.43 [1.00; 2.20]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 12.35    6  0.0546

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/doctor/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/doctor/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/doctor/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/doctor/strict_het_test_before_outlier_extraction.RDS")

#Weird funnel plot, but overall great. Increased heterogeneity due to small SNPs

mr_post <- remove_outlier(mr_steiger_df, rucker) #outliers detected!

saveRDS(mr_post, "output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/asat_gsat_ratio/doctor/instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure                    method nsnp          b        se      pval
# 1    exposure    outcome outcome exposure                  MR Egger    5 1.46469761 0.6524453 0.1104550
# 2    exposure    outcome outcome exposure           Weighted median    5 0.27453851 0.2616673 0.2940911
# 3    exposure    outcome outcome exposure Inverse variance weighted    5 0.31470281 0.2095973 0.1332356
# 4    exposure    outcome outcome exposure               Simple mode    5 0.09224507 0.4289616 0.8402531
# 5    exposure    outcome outcome exposure             Weighted mode    5 0.29336775 0.3437756 0.4415433

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method    Estimate         SE     CI_low      CI_upp           P
# 1  Egger fixed effects -0.09806606 0.03038373 -0.1576171 -0.03851504 0.001248401
# 2 Egger random effects -0.09806606 0.03038373 -0.1576171 -0.03851504 0.999375799
# 
# [[1]]$Q
# Method         Q df          P
# 1   Q_ivw 4.4142178  4 0.35284055
# 2 Q_egger 0.9869833  3 0.80440157
# 3  Q_diff 3.4272344  1 0.06412941
# 
# [[1]]$res
# [1] "A"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_post) #they disagree, but still

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0.0425 [0.0000; 1.8061]; tau = 0.2061 [0.0000; 1.3439]
# I^2 = 7.0% [0.0%; 80.7%]; H = 1.04 [1.00; 2.27]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 4.30    4  0.3670

#That is quite good:

mr_post$labels <- NA
mr_post$units.exposure <- "SD"
mr_post$units.outcome <- "LOR"
 
plotio=mr_plots(mr_post)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/doctor/plots_after_outlier_extraction.svg", height=16, width=16, units="in")
 
saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/doctor/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/doctor/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/doctor/strict_het_test_after_outlier_extraction.RDS")

#Worsened due to less SNPs - we have to be careful!

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

dir.create("output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/broad")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/broad/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp          b        se       pval
# 1    exposure    outcome outcome exposure                  MR Egger    7  1.0810051 0.7350957 0.20136955
# 2    exposure    outcome outcome exposure           Weighted median    7  0.3750075 0.1856779 0.04341795
# 3    exposure    outcome outcome exposure Inverse variance weighted    7  0.3691428 0.2065662 0.07393044
# 4    exposure    outcome outcome exposure               Simple mode    7 -0.1497423 0.3841786 0.71016448
# 5    exposure    outcome outcome exposure             Weighted mode    7  0.4938661 0.2469281 0.09242123

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method    Estimate         SE     CI_low      CI_upp          P
# 1  Egger fixed effects -0.05651996 0.03241071 -0.1200438 0.007003876 0.08118238
# 2 Egger random effects -0.05651996 0.05602004 -0.1663172 0.053277294 0.84349445
# 
# [[1]]$Q
# Method        Q df           P
# 1   Q_ivw 17.97863  6 0.006285833
# 2 Q_egger 14.93756  5 0.010632511
# 3  Q_diff  3.04107  1 0.081182385
# 
# [[1]]$res
# [1] "B"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_steiger_df) #they disagree, but still

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0.2138 [0.0157; 1.6189]; tau = 0.4624 [0.1253; 1.2724]
# I^2 = 63.1% [16.2%; 83.7%]; H = 1.65 [1.09; 2.48]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 16.25    6  0.0125

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/broad/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/broad/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/broad/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/broad/strict_het_test_before_outlier_extraction.RDS")

#Asymmetry and heterogeneity introduced!

mr_post <- remove_outlier(mr_steiger_df, rucker)
 
saveRDS(mr_post, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/broad/instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure                    method nsnp           b        se       pval
# 1    exposure    outcome outcome exposure                  MR Egger    5  0.69459851 0.7074014 0.39855379
# 2    exposure    outcome outcome exposure           Weighted median    5  0.30650897 0.1840518 0.09584483
# 3    exposure    outcome outcome exposure Inverse variance weighted    5  0.33481493 0.1981251 0.09104388
# 4    exposure    outcome outcome exposure               Simple mode    5 -0.04424897 0.3436428 0.90375882
# 5    exposure    outcome outcome exposure             Weighted mode    5  0.32938999 0.2280081 0.22206520

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method    Estimate         SE     CI_low     CI_upp         P
# 1  Egger fixed effects -0.03039064 0.03703729 -0.1029824 0.04220110 0.4119073
# 2 Egger random effects -0.03039064 0.05682925 -0.1417739 0.08099265 0.7035959
# 
# [[1]]$Q
# Method         Q df          P
# 1   Q_ivw 7.7362503  4 0.10173170
# 2 Q_egger 7.0629615  3 0.06991782
# 3  Q_diff 0.6732889  1 0.41190731
# 
# [[1]]$res
# [1] "A"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_post) #they disagree, but still

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0.0858 [0.0000; 2.0459]; tau = 0.2930 [0.0000; 1.4303]
# I^2 = 43.4% [0.0%; 79.2%]; H = 1.33 [1.00; 2.19]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 7.07    4  0.1325

#That is quite good:

mr_post$labels <- NA
mr_post$units.exposure <- "SD"
mr_post$units.outcome <- "LOR"
 
plotio=mr_plots(mr_post)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/broad/plots_after_outlier_extraction.svg", height=16, width=16, units="in")
 
saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/broad/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/broad/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/broad/strict_het_test_after_outlier_extraction.RDS")

################################################################
#Let's first analyse the data for the whole set for asat_gsat_ratio - PCOS# consortium definition
################################################################

#STEP 1: let's use only those that have F>10

mr_df <- liu_pcos_consortium[which(((liu_pcos_consortium$beta.exposure^2)/(liu_pcos_consortium$se.exposure^2)) > 10),] #7

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

dir.create("output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/consortium")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/consortium/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp         b         se        pval
# 1    exposure    outcome outcome exposure                  MR Egger    7 0.2123583 0.16497133 0.254373880
# 2    exposure    outcome outcome exposure           Weighted median    7 0.1642153 0.05638041 0.003583984
# 3    exposure    outcome outcome exposure Inverse variance weighted    7 0.1354352 0.04324741 0.001738343
# 4    exposure    outcome outcome exposure               Simple mode    7 0.1404029 0.08431864 0.146937093
# 5    exposure    outcome outcome exposure             Weighted mode    7 0.1567015 0.06294584 0.047196380

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method    Estimate         SE      CI_low     CI_upp         P
# 1  Egger fixed effects -0.00610708 0.01171014 -0.02905854 0.01684438 0.6020042
# 2 Egger random effects -0.00610708 0.01257111 -0.03074601 0.01853185 0.6864465
# 
# [[1]]$Q
# Method         Q df         P
# 1   Q_ivw 6.0342466  6 0.4193647
# 2 Q_egger 5.7622629  5 0.3300434
# 3  Q_diff 0.2719837  1 0.6020042
# 
# [[1]]$res
# [1] "A"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_steiger_df) #they disagree, but still

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 < 0.0001 [0.0000; 0.0843]; tau = 0.0014 [0.0000; 0.2904]
# I^2 = 0.0% [0.0%; 70.8%]; H = 1.00 [1.00; 1.85]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 5.71    6  0.4560

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/consortium/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/consortium/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/consortium/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/consortium/strict_het_test_before_outlier_extraction.RDS")

#There is a weird outlier...

#STEP 5: removing outliers!

#mr_post <- remove_outlier(mr_steiger_df, rucker) #but not significant outliers!!
  
#saveRDS(mr_post, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/consortium/instruments_after_outlier_removal_df.RDS")
  
#STEP 6: let's run the results for real now:
  
#mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure                    method nsnp           b         se         pval
# 1    exposure    outcome outcome exposure                  MR Egger   17 -0.14328950 0.08244983 0.1027083675
# 2    exposure    outcome outcome exposure           Weighted median   17 -0.13485038 0.03563801 0.0001543949
# 3    exposure    outcome outcome exposure Inverse variance weighted   17 -0.08847733 0.02558887 0.0005449109
# 4    exposure    outcome outcome exposure               Simple mode   17 -0.14458503 0.06459219 0.0397612828
# 5    exposure    outcome outcome exposure             Weighted mode   17 -0.14160290 0.05607439 0.0224993239

# #Let's check the plot:
  
#rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method    Estimate          SE       CI_low     CI_upp         P
# 1  Egger fixed effects 0.005067152 0.007188351 -0.009021757 0.01915606 0.4808653
# 2 Egger random effects 0.005067152 0.007188351 -0.009021757 0.01915606 0.2404326
# 
# [[1]]$Q
# Method          Q df         P
# 1   Q_ivw 15.2523250 16 0.5062413
# 2 Q_egger 14.7632668 15 0.4686002
# 3  Q_diff  0.4890582  1 0.4843477
# 
# [[1]]$res
# [1] "A"

# #Here there is heterogeneity and pleitropy

#strict_het <- compute_strict_het(mr_post) #they disagree, but still

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 < 0.0001 [0.0000; 0.0167]; tau = 0.0017 [0.0000; 0.1293]
# I^2 = 0.0% [0.0%; 51.1%]; H = 1.00 [1.00; 1.43]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 14.60   16  0.5542
#   
#That is quite good:
  
#mr_post$labels <- NA
#mr_post$units.exposure <- "SD"
#mr_post$units.outcome <- "LOR"
  
#plotio=mr_plots(mr_post)
#ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/consortium/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

#saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/consortium/mr_res_after_outlier_extraction.RDS")
#saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/consortium/sensitivity_test_after_outlier_extraction.RDS")
#saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/consortium/strict_het_test_after_outlier_extraction.RDS")

#I think we might have now more asymmetry, despite heterogeneity not popping up...

#mr_post_2 <- remove_outlier(mr_post, rucker) #but not significant outliers!! - we are fine

###############################################################
#Let's first analyse the data for the whole set for asat_gsat_ratio - PCOS# age adjusted
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

dir.create("output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/age_adjusted")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/age_adjusted/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp         b        se       pval
# 1    exposure    outcome outcome exposure                  MR Egger    7 0.7816093 0.6055845 0.25327420
# 2    exposure    outcome outcome exposure           Weighted median    7 0.2448390 0.2071693 0.23727282
# 3    exposure    outcome outcome exposure Inverse variance weighted    7 0.3785354 0.1587913 0.01713238
# 4    exposure    outcome outcome exposure               Simple mode    7 0.2442107 0.3265429 0.48281195
# 5    exposure    outcome outcome exposure             Weighted mode    7 0.4352001 0.2918372 0.18649202

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method    Estimate         SE     CI_low     CI_upp         P
# 1  Egger fixed effects -0.03303176 0.04203600 -0.1154208 0.04935729 0.4319865
# 2 Egger random effects -0.03303176 0.04772261 -0.1265664 0.06050284 0.7555821
# 
# [[1]]$Q
# Method         Q df         P
# 1   Q_ivw 7.0617748  6 0.3151721
# 2 Q_egger 6.4442981  5 0.2653539
# 3  Q_diff 0.6174767  1 0.4319865
# 
# [[1]]$res
# [1] "A"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_steiger_df) #they disagree, but still

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0.0431 [0.0000; 0.9083]; tau = 0.2077 [0.0000; 0.9530]
# I^2 = 11.7% [0.0%; 74.2%]; H = 1.06 [1.00; 1.97]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 6.79    6  0.3404

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/age_adjusted//plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/age_adjusted//mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/age_adjusted//sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/age_adjusted//strict_het_test_before_outlier_extraction.RDS")

#We have a weird outlier! It does not affect things much, but let's see...

mr_post <- remove_outlier(mr_steiger_df, rucker) #no significant outliers. 

#saveRDS(mr_post, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/age_adjusted//instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

#mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure                    method nsnp             b         se      pval
# 1    exposure    outcome outcome exposure                  MR Egger   17 -0.2653283742 0.28216070 0.3619384
# 2    exposure    outcome outcome exposure           Weighted median   17 -0.0008582422 0.11885157 0.9942384
# 3    exposure    outcome outcome exposure Inverse variance weighted   17 -0.1159567239 0.08407172 0.1678147
# 4    exposure    outcome outcome exposure               Simple mode   17  0.0640096156 0.22443312 0.7791454
# 5    exposure    outcome outcome exposure             Weighted mode   17  0.0745336637 0.21030649 0.7276631

# #Let's check the plot:
  
#rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method   Estimate         SE      CI_low     CI_upp         P
# 1  Egger fixed effects 0.01399889 0.02432745 -0.03368205 0.06167982 0.5649966
# 2 Egger random effects 0.01399889 0.02432745 -0.03368205 0.06167982 0.2824983
# 
# [[1]]$Q
# Method          Q df         P
# 1   Q_ivw 14.2396565 16 0.5808651
# 2 Q_egger 13.9321040 15 0.5306847
# 3  Q_diff  0.3075525  1 0.5791860
# 
# [[1]]$res
# [1] "A"

# #Here there is heterogeneity and pleitropy

#strict_het <- compute_strict_het(mr_post) 

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0 [0.0000; 0.1319]; tau = 0 [0.0000; 0.3631]
# I^2 = 0.0% [0.0%; 51.1%]; H = 1.00 [1.00; 1.43]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 13.69   16  0.6221

#That is quite good:

#mr_post$labels <- NA
#mr_post$units.exposure <- "SD"
#mr_post$units.outcome <- "LOR"

#plotio=mr_plots(mr_post)
#ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/age_adjusted/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

#saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/age_adjusted/mr_res_after_outlier_extraction.RDS")
#saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/age_adjusted/sensitivity_test_after_outlier_extraction.RDS")
#saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/age_adjusted/strict_het_test_after_outlier_extraction.RDS")

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

dir.create("output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/age_and_bmi_adjusted")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/age_and_bmi_adjusted/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)
 
# id.exposure id.outcome outcome exposure                    method nsnp         b        se        pval
# 1    exposure    outcome outcome exposure                  MR Egger    7 0.7191835 1.1316279 0.553027349
# 2    exposure    outcome outcome exposure           Weighted median    7 0.7451571 0.2503578 0.002916842
# 3    exposure    outcome outcome exposure Inverse variance weighted    7 0.4381075 0.2854198 0.124794071
# 4    exposure    outcome outcome exposure               Simple mode    7 0.7376079 0.3875349 0.105680390
# 5    exposure    outcome outcome exposure             Weighted mode    7 0.7484847 0.3145003 0.054774291

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method    Estimate         SE     CI_low     CI_upp         P
# 1  Egger fixed effects -0.02309069 0.04949426 -0.1200977 0.07391627 0.6408342
# 2 Egger random effects -0.02309069 0.08939423 -0.1983002 0.15211878 0.6019130
# 
# [[1]]$Q
# Method          Q df           P
# 1   Q_ivw 16.5286075  6 0.011181109
# 2 Q_egger 16.3109547  5 0.006010146
# 3  Q_diff  0.2176528  1 0.640834222
# 
# [[1]]$res
# [1] "B"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_steiger_df) #they disagree, but still

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0.3199 [0.0108; 4.4344]; tau = 0.5656 [0.1041; 2.1058]
# I^2 = 60.7% [9.9%; 82.8%]; H = 1.59 [1.05; 2.41]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 15.26    6  0.0183

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/age_and_bmi_adjusted/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/age_and_bmi_adjusted/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/age_and_bmi_adjusted/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/age_and_bmi_adjusted/strict_het_test_before_outlier_extraction.RDS")

#We have a weird outlier

mr_post <- remove_outlier(mr_steiger_df, rucker) #no significant outliers. 

saveRDS(mr_post, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/age_and_bmi_adjusted//instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure                    method nsnp         b        se         pval
# 1    exposure    outcome outcome exposure                  MR Egger    6 0.8578994 0.7756736 0.3307646640
# 2    exposure    outcome outcome exposure           Weighted median    6 0.7614862 0.2545060 0.0027714144
# 3    exposure    outcome outcome exposure Inverse variance weighted    6 0.7284371 0.2152925 0.0007157443
# 4    exposure    outcome outcome exposure               Simple mode    6 0.7445048 0.3363783 0.0777812737
# 5    exposure    outcome outcome exposure             Weighted mode    6 0.7535259 0.2780267 0.0422628342

# #Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method    Estimate         SE     CI_low     CI_upp         P
# 1  Egger fixed effects -0.01076057 0.04964446 -0.1080619 0.08654079 0.8284011
# 2 Egger random effects -0.01076057 0.06131416 -0.1309341 0.10941298 0.5696562
# 
# [[1]]$Q
# Method          Q df         P
# 1   Q_ivw 6.14852824  5 0.2920336
# 2 Q_egger 6.10154651  4 0.1916920
# 3  Q_diff 0.04698173  1 0.8284011
# 
# [[1]]$res
# [1] "A"

strict_het <- compute_strict_het(mr_post) 

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 < 0.0001 [0.0000; 4.2532]; tau = 0.0015 [0.0000; 2.0623]
# I^2 = 10.7% [0.0%; 77.3%]; H = 1.06 [1.00; 2.10]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 5.60    5  0.3470

#That is quite good:

mr_post$labels <- NA
mr_post$units.exposure <- "SD"
mr_post$units.outcome <- "LOR"

plotio=mr_plots(mr_post)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/age_and_bmi_adjusted/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/age_and_bmi_adjusted/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/age_and_bmi_adjusted/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/age_and_bmi_adjusted/strict_het_test_after_outlier_extraction.RDS")

###############################################################
#Let's first analyse the data for the whole set for bmi - PCOS# Venkatesh
###############################################################

#STEP 1: let's use only those that have F>10

mr_df <- liu_pcos_venkatesh[which(((liu_pcos_venkatesh$beta.exposure^2)/(liu_pcos_venkatesh$se.exposure^2)) > 10),] #7

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

dir.create("output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/venkatesh")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/venkatesh/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp          b        se      pval
# 1    exposure    outcome outcome exposure                  MR Egger    7 -0.2459888 0.9287558 0.8016904
# 2    exposure    outcome outcome exposure           Weighted median    7  0.3054048 0.2882472 0.2893613
# 3    exposure    outcome outcome exposure Inverse variance weighted    7  0.2062848 0.2231624 0.3552934
# 4    exposure    outcome outcome exposure               Simple mode    7 -0.1403455 0.5019740 0.7891823
# 5    exposure    outcome outcome exposure             Weighted mode    7  0.5417476 0.3703298 0.1938252

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method   Estimate         SE     CI_low    CI_upp         P
# 1  Egger fixed effects 0.03496663 0.06889114 -0.1000575 0.1699908 0.6117595
# 2 Egger random effects 0.03496663 0.06965371 -0.1015521 0.1714854 0.3078314
# 
# [[1]]$Q
# Method         Q df         P
# 1   Q_ivw 5.3689246  6 0.4974383
# 2 Q_egger 5.1113039  5 0.4024482
# 3  Q_diff 0.2576207  1 0.6117595
# 
# [[1]]$res
# [1] "A"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_steiger_df) #they disagree, but still

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0.0653 [0.0000; 1.0423]; tau = 0.2555 [0.0000; 1.0209]
# I^2 = 0.0% [0.0%; 70.8%]; H = 1.00 [1.00; 1.85]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 5.12    6  0.5288

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/venkatesh/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/venkatesh/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/venkatesh/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/venkatesh/strict_het_test_before_outlier_extraction.RDS")

#Asymmetry! And some outliers!

mr_post <- remove_outlier(mr_steiger_df, rucker) #no significant outlier

#saveRDS(mr_post, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/venkatesh/instruments_after_outlier_removal_df.RDS")

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
# ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/venkatesh/plots_after_outlier_extraction.svg", height=16, width=16, units="in")
# 
# saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/venkatesh/mr_res_after_outlier_extraction.RDS")
# saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/venkatesh/sensitivity_test_after_outlier_extraction.RDS")
# saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asat_gsat_ratio/venkatesh/strict_het_test_after_outlier_extraction.RDS")