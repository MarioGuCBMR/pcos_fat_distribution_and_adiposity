##############
#INTRODUCTION#
##############

#This code performs MR analyses between WHRadjBMI and PCOS - replicating the results from Liu et al.

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

liu_pcos_doctor <- fread("output/2_replicating_liu_et_al/4_mr_fat_depots/1_expo_outcome_df/asatadj_PCOS.txt")
liu_pcos_broad <- fread("output/2_replicating_liu_et_al/4_mr_fat_depots/1_expo_outcome_df/asatadj_PCOS (Broad).txt")
liu_pcos_consortium <- fread("output/2_replicating_liu_et_al/4_mr_fat_depots/1_expo_outcome_df/asatadj_PCOS (Consortium).txt")
liu_pcos_meta_analysis <- fread("output/2_replicating_liu_et_al/4_mr_fat_depots/1_expo_outcome_df/asatadj_PCOS (Day et al).txt")
liu_pcos_venkatesh <- fread("output/2_replicating_liu_et_al/4_mr_fat_depots/1_expo_outcome_df/asatadj_PCOS (Venkatesh et al).txt")
liu_pcos_adj_age <- fread("output/2_replicating_liu_et_al/4_mr_fat_depots/1_expo_outcome_df/asatadj_PCOS (adj age).txt")
liu_pcos_adj_age_bmi <- fread("output/2_replicating_liu_et_al/4_mr_fat_depots/1_expo_outcome_df/asatadj_PCOS (adj age+BMI).txt")

##############################################################
#STEP 2: let's first rerun the analyses for the meta-analysis#
##############################################################

mr_df <- liu_pcos_meta_analysis[which(((as.numeric(liu_pcos_meta_analysis$beta.exposure)^2)/(as.numeric(liu_pcos_meta_analysis$se.exposure)^2)) > 10),] #3/3

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

mr_df$id.exposure <- "asatadj"
mr_df$id.outcome <- "PCOS (meta-analysis)"
mr_df$exposure <- "asatadj"
mr_df$outcome <- "PCOS (meta-analysis)"
mr_df$units.exposure <- "SD"
mr_df$units.outcome <- "LOR"

mr_steiger_df <- TwoSampleMR::steiger_filtering(mr_df)
mr_steiger_df <- mr_steiger_df[which(mr_steiger_df$steiger_dir == TRUE),] #1
mr_steiger_df$mr_keep <- TRUE

#Let's save the data:

dir.create("output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/")
dir.create("output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/asatadj/")
dir.create("output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/asatadj/meta_analysis/")
saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/asatadj/meta_analysis/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure           id.outcome              outcome exposure                    method nsnp        b        se      pval
# 1     asatadj PCOS (meta-analysis) PCOS (meta-analysis)  asatadj                  MR Egger    3 9.032750 5.9993123 0.3732339
# 2     asatadj PCOS (meta-analysis) PCOS (meta-analysis)  asatadj           Weighted median    3 0.313918 0.3628970 0.3870204
# 3     asatadj PCOS (meta-analysis) PCOS (meta-analysis)  asatadj Inverse variance weighted    3 0.192141 0.3311177 0.5617257
# 4     asatadj PCOS (meta-analysis) PCOS (meta-analysis)  asatadj               Simple mode    3 0.423612 0.4855518 0.4749649
# 5     asatadj PCOS (meta-analysis) PCOS (meta-analysis)  asatadj             Weighted mode    3 0.423612 0.4927838 0.4805801

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method   Estimate       SE    CI_low     CI_upp           P
# 1  Egger fixed effects -0.5801258 0.221635 -1.014522 -0.1457292 0.008858089
# 2 Egger random effects -0.5801258 0.221635 -1.014522 -0.1457292 0.995570955
# 
# [[1]]$Q
# Method         Q df         P
# 1   Q_ivw 2.4945526  2 0.2872862
#2 Q_egger 0.3177281  1 0.5729762
#3  Q_diff 2.1768245  1 0.1401033

#[[1]]$res
#[1] "A"

strict_het <- compute_strict_het(mr_steiger_df)

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0.0368 [0.0000; 13.8718]; tau = 0.1919 [0.0000; 3.7245]
# I^2 = 17.5% [0.0%; 91.4%]; H = 1.10 [1.00; 3.41]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 2.42    2  0.2978

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/asatadj/meta_analysis/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/asatadj/meta_analysis/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/asatadj/meta_analysis//sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/asatadj/meta_analysis//strict_het_test_before_outlier_extraction.RDS")

#Egger is quite biased, these two variants have no effect...
#Let's check if RadialMR removes this...

#mr_post <- remove_outlier(mr_steiger_df, rucker) #it is not working

#saveRDS(mr_post, "output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/asatadj/meta_analysis/instruments_after_outlier_removal_df.RDS")

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
# ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/meta_analysis/plots_after_outlier_extraction.svg", height=16, width=16, units="in")
# 
# saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/meta_analysis/mr_res_after_outlier_extraction.RDS")
# saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/meta_analysis/sensitivity_test_after_outlier_extraction.RDS")
# saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/meta_analysis/strict_het_test_after_outlier_extraction.RDS")

###############################################################
#Let's first analyse the data for the whole set for asatadj - PCOS#
###############################################################

#STEP 1: let's use only those that have F>10

mr_df <- liu_pcos_doctor[which(((liu_pcos_doctor$beta.exposure^2)/(liu_pcos_doctor$se.exposure^2)) > 10),] #2/2

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

dir.create("output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/")
dir.create("output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/doctor")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/doctor/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp           b        se      pval
# 1    exposure    outcome outcome exposure                  MR Egger    3 -2.02369827 6.3965287 0.8049329
# 2    exposure    outcome outcome exposure           Weighted median    3 -0.01478666 0.3239004 0.9635877
# 3    exposure    outcome outcome exposure Inverse variance weighted    3 -0.16963892 0.2716792 0.5323591
# 4    exposure    outcome outcome exposure               Simple mode    3  0.06136722 0.4186801 0.8969094
# 5    exposure    outcome outcome exposure             Weighted mode    3  0.07521477 0.4424394 0.8806510

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method  Estimate        SE     CI_low    CI_upp         P
# 1  Egger fixed effects 0.1216959 0.3613350 -0.5865077 0.8298995 0.7362712
# 2 Egger random effects 0.1216959 0.4193421 -0.7001996 0.9435914 0.3858290
# 
# [[1]]$Q
# Method        Q df         P
# 1   Q_ivw 1.460274  2 0.4818430
# 2 Q_egger 1.346843  1 0.2458308
# 3  Q_diff 0.113431  1 0.7362712
# 
# [[1]]$res
# [1] "A"

strict_het <- compute_strict_het(mr_steiger_df)

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0 [0.0000; 5.7807]; tau = 0 [0.0000; 2.4043]
# I^2 = 0.0% [0.0%; 89.6%]; H = 1.00 [1.00; 3.10]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 1.43    2  0.4902

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/doctor/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/doctor/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/doctor/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/doctor/strict_het_test_before_outlier_extraction.RDS")

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

dir.create("output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/broad")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/broad/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp          b         se      pval
# 1    exposure    outcome outcome exposure                  MR Egger    3  7.8846085 10.9080738 0.6015525
# 2    exposure    outcome outcome exposure           Weighted median    3 -0.1661824  0.2723815 0.5417897
# 3    exposure    outcome outcome exposure Inverse variance weighted    3  0.2067391  0.4648997 0.6565394
# 4    exposure    outcome outcome exposure               Simple mode    3 -0.2994467  0.2844662 0.4029074
# 5    exposure    outcome outcome exposure             Weighted mode    3 -0.2782800  0.3087414 0.4625379

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method  Estimate       SE     CI_low      CI_upp          P
# 1  Egger fixed effects -0.503987 0.246510 -0.9871377 -0.02083631 0.04090524
# 2 Egger random effects -0.503987 0.715153 -1.9056611  0.89768709 0.75950966
# 
# [[1]]$Q
# Method         Q df           P
# 1   Q_ivw 12.596384  2 0.001839628
# 2 Q_egger  8.416448  1 0.003718415
# 3  Q_diff  4.179936  1 0.040905245
# 
# [[1]]$res
# [1] "D"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_steiger_df) #they disagree, but still

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0.5109 [0.0530; 24.9577]; tau = 0.7148 [0.2302; 4.9958]
# I^2 = 81.3% [41.7%; 94.0%]; H = 2.31 [1.31; 4.08]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 10.68    2  0.0048

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/broad/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/broad/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/broad/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/broad/strict_het_test_before_outlier_extraction.RDS")

#A single variant seems to be introducing quite a bit of heterogeneity.
#I do not think this might affect results, but let's see.

# mr_post <- remove_outlier(mr_steiger_df, rucker)
# 
# saveRDS(mr_post, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/broad/instruments_after_outlier_removal_df.RDS")

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
# ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/broad/plots_after_outlier_extraction.svg", height=16, width=16, units="in")
# 
# saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/broad/mr_res_after_outlier_extraction.RDS")
# saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/broad/sensitivity_test_after_outlier_extraction.RDS")
# saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/broad/strict_het_test_after_outlier_extraction.RDS")

################################################################
#Let's first analyse the data for the whole set for asatadj - PCOS# consortium definition
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

dir.create("output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/consortium")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/consortium/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp           b         se      pval
# 1    exposure    outcome outcome exposure                  MR Egger    3 -0.47595366 2.18058536 0.8631914
# 2    exposure    outcome outcome exposure           Weighted median    3  0.12968077 0.08823796 0.1416509
# 3    exposure    outcome outcome exposure Inverse variance weighted    3  0.08938105 0.07854976 0.2551661
# 4    exposure    outcome outcome exposure               Simple mode    3  0.16697656 0.11851047 0.2942102
# 5    exposure    outcome outcome exposure             Weighted mode    3  0.17078327 0.12322998 0.3000814

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method   Estimate         SE     CI_low    CI_upp         P
# 1  Egger fixed effects 0.03710582 0.08914174 -0.1376088 0.2118204 0.6772224
# 2 Egger random effects 0.03710582 0.14294892 -0.2430689 0.3172806 0.3975962
# 
# [[1]]$Q
# Method         Q df         P
# 1   Q_ivw 2.7448463  2 0.2534920
# 2 Q_egger 2.5715769  1 0.1087992
# 3  Q_diff 0.1732694  1 0.6772224
# 
# [[1]]$res
# [1] "A"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_steiger_df) #they disagree, but still

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0.0054 [0.0000; 0.6723]; tau = 0.0734 [0.0000; 0.8200]
# I^2 = 25.1% [0.0%; 92.2%]; H = 1.16 [1.00; 3.58]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 2.67    2  0.2633

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/consortium/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/consortium/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/consortium/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/consortium/strict_het_test_before_outlier_extraction.RDS")

#There is pleiotropy.

#STEP 5: removing outliers!

#mr_post <- remove_outlier(mr_steiger_df, rucker) #but not significant outliers!!
  
#saveRDS(mr_post, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/consortium/instruments_after_outlier_removal_df.RDS")
  
#STEP 6: let's run the results for real now:
  
#mr_res <- TwoSampleMR::mr(mr_post)

# #Let's check the plot:
  
#rucker <- TwoSampleMR::mr_rucker(mr_post)

# #Here there is heterogeneity and pleitropy

#strict_het <- compute_strict_het(mr_post) #they disagree, but still
  
#That is quite good:
  
#mr_post$labels <- NA
#mr_post$units.exposure <- "SD"
#mr_post$units.outcome <- "LOR"
  
#plotio=mr_plots(mr_post)
#ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/consortium/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

#saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/consortium/mr_res_after_outlier_extraction.RDS")
#saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/consortium/sensitivity_test_after_outlier_extraction.RDS")
#saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/consortium/strict_het_test_after_outlier_extraction.RDS")

###############################################################
#Let's first analyse the data for the whole set for asatadj - PCOS# age adjusted
###############################################################

#STEP 1: let's use only those that have F>10

mr_df <- liu_pcos_adj_age[which(((liu_pcos_adj_age$beta.exposure^2)/(liu_pcos_adj_age$se.exposure^2)) > 10),] #3
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

dir.create("output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/age_adjusted")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/age_adjusted/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp           b        se      pval
# 1    exposure    outcome outcome exposure                  MR Egger    3 -4.84159764 4.5471258 0.4800393
# 2    exposure    outcome outcome exposure           Weighted median    3 -0.19547901 0.2706089 0.4700688
# 3    exposure    outcome outcome exposure Inverse variance weighted    3 -0.05920374 0.2221557 0.7898567
# 4    exposure    outcome outcome exposure               Simple mode    3 -0.20127525 0.3234403 0.5972392
# 5    exposure    outcome outcome exposure             Weighted mode    3 -0.20127525 0.3004542 0.5719070

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method  Estimate         SE    CI_low    CI_upp             P
# 1  Egger fixed effects 0.3142084 0.01252995 0.2896502 0.3387667 8.955628e-139
# 2 Egger random effects 0.3142084 0.01252995 0.2896502 0.3387667 4.477814e-139
# 
# [[1]]$Q
# Method           Q df         P
# 1   Q_ivw 1.110566857  2 0.5739096
# 2 Q_egger 0.001763265  1 0.9665057
# 3  Q_diff 1.108803592  1 0.2923422
# 
# [[1]]$res
# [1] "A"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_steiger_df) #they disagree, but still

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0 [0.0000; 3.4925]; tau = 0 [0.0000; 1.8688]
# I^2 = 0.0% [0.0%; 89.6%]; H = 1.00 [1.00; 3.10]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 1.09    2  0.5785

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/age_adjusted//plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/age_adjusted//mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/age_adjusted//sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/age_adjusted//strict_het_test_before_outlier_extraction.RDS")

#There is heterogeneity and a bit of asymmetry. 
#Not enough, but let's check RadialMR.

#mr_post <- remove_outlier(mr_steiger_df, rucker) #no significant outliers. 

#saveRDS(mr_post, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/age_adjusted//instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

#mr_res <- TwoSampleMR::mr(mr_post)

# #Let's check the plot:
  
#rucker <- TwoSampleMR::mr_rucker(mr_post)

# #Here there is heterogeneity and pleitropy

#strict_het <- compute_strict_het(mr_post) #they disagree, but still

#That is quite good:

#mr_post$labels <- NA
#mr_post$units.exposure <- "SD"
#mr_post$units.outcome <- "LOR"

#plotio=mr_plots(mr_post)
#ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/age_adjusted/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

#saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/age_adjusted/mr_res_after_outlier_extraction.RDS")
#saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/age_adjusted/sensitivity_test_after_outlier_extraction.RDS")
#saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/age_adjusted/strict_het_test_after_outlier_extraction.RDS")

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

dir.create("output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/age_and_bmi_adjusted")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/age_and_bmi_adjusted/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp          b        se      pval
# 1    exposure    outcome outcome exposure                  MR Egger    3 -7.1347310 5.3486117 0.4095264
# 2    exposure    outcome outcome exposure           Weighted median    3 -0.3897719 0.3262077 0.2321425
# 3    exposure    outcome outcome exposure Inverse variance weighted    3 -0.3345332 0.2611034 0.2001133
# 4    exposure    outcome outcome exposure               Simple mode    3 -0.4490406 0.4209693 0.3978263
# 5    exposure    outcome outcome exposure             Weighted mode    3 -0.4722239 0.4174847 0.3753905

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method  Estimate        SE     CI_low    CI_upp          P
# 1  Egger fixed effects 0.4468157 0.2091162 0.03695553 0.8566759 0.03262350
# 2 Egger random effects 0.4468157 0.2091162 0.03695553 0.8566759 0.01631175
# 
# [[1]]$Q
# Method         Q df         P
# 1   Q_ivw 1.9752142  2 0.3724669
# 2 Q_egger 0.3549078  1 0.5513480
# 3  Q_diff 1.6203064  1 0.2030491
# 
# [[1]]$res
# [1] "A"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_steiger_df) #they disagree, but still

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 < 0.0001 [0.0000; 8.4407]; tau = 0.0008 [0.0000; 2.9053]
# I^2 = 0.0% [0.0%; 89.6%]; H = 1.00 [1.00; 3.10]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 1.92    2  0.3826

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/age_and_bmi_adjusted/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/age_and_bmi_adjusted/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/age_and_bmi_adjusted/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/age_and_bmi_adjusted/strict_het_test_before_outlier_extraction.RDS")

#We have pleiotropy, let's check this out:

#mr_post <- remove_outlier(mr_steiger_df, rucker) #no significant outliers. 

#saveRDS(mr_post, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/age_and_bmi_adjusted//instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

#mr_res <- TwoSampleMR::mr(mr_post)

# #Let's check the plot:

#rucker <- TwoSampleMR::mr_rucker(mr_post)

# #Here there is heterogeneity and pleitropy

#strict_het <- compute_strict_het(mr_post) #they disagree, but still

#That is quite good:

#mr_post$labels <- NA
#mr_post$units.exposure <- "SD"
#mr_post$units.outcome <- "LOR"

#plotio=mr_plots(mr_post)
#ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/age_and_bmi_adjusted/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

#saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/age_and_bmi_adjusted/mr_res_after_outlier_extraction.RDS")
#saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/age_and_bmi_adjusted/sensitivity_test_after_outlier_extraction.RDS")
#saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/age_and_bmi_adjusted/strict_het_test_after_outlier_extraction.RDS")

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

dir.create("output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/venkatesh")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/venkatesh/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp           b        se      pval
# 1    exposure    outcome outcome exposure                  MR Egger    3 -4.65264095 7.6037937 0.6504244
# 2    exposure    outcome outcome exposure           Weighted median    3  0.18294445 0.3633596 0.6146265
# 3    exposure    outcome outcome exposure Inverse variance weighted    3  0.01328892 0.3153496 0.9663869
# 4    exposure    outcome outcome exposure               Simple mode    3  0.26134304 0.4758622 0.6379963
# 5    exposure    outcome outcome exposure             Weighted mode    3  0.20716724 0.4095435 0.6632073

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method  Estimate        SE     CI_low   CI_upp         P
# 1  Egger fixed effects 0.2939233 0.3935964 -0.4775115 1.065358 0.4552065
# 2 Egger random effects 0.2939233 0.4783912 -0.6437062 1.231553 0.2694757
# 
# [[1]]$Q
# Method         Q df         P
# 1   Q_ivw 2.0349397  2 0.3615085
# 2 Q_egger 1.4772844  1 0.2241999
# 3  Q_diff 0.5576553  1 0.4552065
# 
# [[1]]$res
# [1] "A"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_steiger_df) #they disagree, but still

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 < 0.0001 [0.0000; 19.2984]; tau = 0.0005 [0.0000; 4.3930]
# I^2 = 0.0% [0.0%; 89.6%]; H = 1.00 [1.00; 3.10]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 1.98    2  0.3722

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/venkatesh/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/venkatesh/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/venkatesh/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/venkatesh/strict_het_test_before_outlier_extraction.RDS")

#Pleiotropy! Let's remove the outliers

#mr_post <- remove_outlier(mr_steiger_df, rucker) #does not work!

#saveRDS(mr_post, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/venkatesh/instruments_after_outlier_removal_df.RDS")

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
# ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/venkatesh/plots_after_outlier_extraction.svg", height=16, width=16, units="in")
# 
# saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/venkatesh/mr_res_after_outlier_extraction.RDS")
# saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/venkatesh/sensitivity_test_after_outlier_extraction.RDS")
# saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/asatadj/venkatesh/strict_het_test_after_outlier_extraction.RDS")