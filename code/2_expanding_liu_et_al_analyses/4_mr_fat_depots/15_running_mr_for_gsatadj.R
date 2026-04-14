##############
#INTRODUCTION#
##############

#This code performs MR analyses between gsatadjBMI and PCOS - replicating the results from Liu et al.

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

liu_pcos_doctor <- fread("output/2_replicating_liu_et_al/4_mr_fat_depots/1_expo_outcome_df/gsatadj_PCOS.txt")
liu_pcos_broad <- fread("output/2_replicating_liu_et_al/4_mr_fat_depots/1_expo_outcome_df/gsatadj_PCOS (Broad).txt")
liu_pcos_consortium <- fread("output/2_replicating_liu_et_al/4_mr_fat_depots/1_expo_outcome_df/gsatadj_PCOS (Consortium).txt")
liu_pcos_meta_analysis <- fread("output/2_replicating_liu_et_al/4_mr_fat_depots/1_expo_outcome_df/gsatadj_PCOS (Day et al).txt")
liu_pcos_venkatesh <- fread("output/2_replicating_liu_et_al/4_mr_fat_depots/1_expo_outcome_df/gsatadj_PCOS (Venkatesh et al).txt")
liu_pcos_adj_age <- fread("output/2_replicating_liu_et_al/4_mr_fat_depots/1_expo_outcome_df/gsatadj_PCOS (adj age).txt")
liu_pcos_adj_age_bmi <- fread("output/2_replicating_liu_et_al/4_mr_fat_depots/1_expo_outcome_df/gsatadj_PCOS (adj age+BMI).txt")

##############################################################
#STEP 2: let's first rerun the analyses for the meta-analysis#
##############################################################

mr_df <- liu_pcos_meta_analysis[which(((as.numeric(liu_pcos_meta_analysis$beta.exposure)^2)/(as.numeric(liu_pcos_meta_analysis$se.exposure)^2)) > 10),] #19/19

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

mr_df$id.exposure <- "gsatadj"
mr_df$id.outcome <- "PCOS (meta-analysis)"
mr_df$exposure <- "gsatadj"
mr_df$outcome <- "PCOS (meta-analysis)"
mr_df$units.exposure <- "SD"
mr_df$units.outcome <- "LOR"

mr_steiger_df <- TwoSampleMR::steiger_filtering(mr_df)
mr_steiger_df <- mr_steiger_df[which(mr_steiger_df$steiger_dir == TRUE),] #19
mr_steiger_df$mr_keep <- TRUE

#Let's save the data:

dir.create("output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/")
dir.create("output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/gsatadj/")
dir.create("output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/gsatadj/meta_analysis/")
saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/gsatadj/meta_analysis/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure           id.outcome              outcome exposure                    method nsnp           b        se      pval
# 1     gsatadj PCOS (meta-analysis) PCOS (meta-analysis)  gsatadj                  MR Egger   19  0.04977140 0.3240623 0.8797447
# 2     gsatadj PCOS (meta-analysis) PCOS (meta-analysis)  gsatadj           Weighted median   19 -0.09852377 0.1365606 0.4706231
# 3     gsatadj PCOS (meta-analysis) PCOS (meta-analysis)  gsatadj Inverse variance weighted   19 -0.12452248 0.1009341 0.2173145
# 4     gsatadj PCOS (meta-analysis) PCOS (meta-analysis)  gsatadj               Simple mode   19 -0.09835536 0.2310956 0.6754413
# 5     gsatadj PCOS (meta-analysis) PCOS (meta-analysis)  gsatadj             Weighted mode   19 -0.09835536 0.2152474 0.6531797

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method    Estimate         SE      CI_low     CI_upp         P
# 1  Egger fixed effects -0.01636768 0.02359425 -0.06261156 0.02987621 0.4878612
# 2 Egger random effects -0.01636768 0.02359425 -0.06261156 0.02987621 0.7560694
# 
# [[1]]$Q
# Method        Q df         P
# 1   Q_ivw 11.63685 18 0.8654084
# 2 Q_egger 11.31650 17 0.8396459
# 3  Q_diff  0.32035  1 0.5713974
# 
# [[1]]$res
# [1] "A"

strict_het <- compute_strict_het(mr_steiger_df)

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0 [0.0000; 0.0826]; tau = 0 [0.0000; 0.2874]
# I^2 = 0.0% [0.0%; 48.9%]; H = 1.00 [1.00; 1.40]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 11.26   18  0.8831

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/gsatadj/meta_analysis/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/gsatadj/meta_analysis/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/gsatadj/meta_analysis//sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/gsatadj/meta_analysis//strict_het_test_before_outlier_extraction.RDS")

#No heterogeneity

#mr_post <- remove_outlier(mr_steiger_df, rucker) #no significant outliers!

#saveRDS(mr_post, "output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/gsatadj/meta_analysis/instruments_after_outlier_removal_df.RDS")

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
# ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/meta_analysis/plots_after_outlier_extraction.svg", height=16, width=16, units="in")
# 
# saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/meta_analysis/mr_res_after_outlier_extraction.RDS")
# saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/meta_analysis/sensitivity_test_after_outlier_extraction.RDS")
# saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/meta_analysis/strict_het_test_after_outlier_extraction.RDS")

###############################################################
#Let's first analyse the data for the whole set for gsatadj - PCOS#
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

dir.create("output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/")
dir.create("output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/doctor")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/doctor/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp          b        se       pval
# 1    exposure    outcome outcome exposure                  MR Egger   19 -0.2467544 0.3611379 0.50364336
# 2    exposure    outcome outcome exposure           Weighted median   19 -0.1811905 0.1457573 0.21383199
# 3    exposure    outcome outcome exposure Inverse variance weighted   19 -0.2188216 0.1107080 0.04809076
# 4    exposure    outcome outcome exposure               Simple mode   19 -0.2550089 0.2657584 0.34999376
# 5    exposure    outcome outcome exposure             Weighted mode   19 -0.1677994 0.2511683 0.51256151

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method    Estimate         SE      CI_low     CI_upp         P
# 1  Egger fixed effects 0.002648228 0.02777044 -0.05178083 0.05707728 0.9240278
# 2 Egger random effects 0.002648228 0.03249110 -0.06103315 0.06632961 0.4675197
# 
# [[1]]$Q
# Method            Q df         P
# 1   Q_ivw 23.279944406 18 0.1800604
# 2 Q_egger 23.270850611 17 0.1406390
# 3  Q_diff  0.009093795  1 0.9240278
# 
# [[1]]$res
# [1] "A"

strict_het <- compute_strict_het(mr_steiger_df)

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0.0404 [0.0000; 0.3488]; tau = 0.2009 [0.0000; 0.5906]
# I^2 = 17.6% [0.0%; 52.3%]; H = 1.10 [1.00; 1.45]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 21.86   18  0.2384

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/doctor/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/doctor/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/doctor/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/doctor/strict_het_test_before_outlier_extraction.RDS")

#Slight heterogeneity and no pleiotropy, but leave-one out reveals that removing some variants might affect the results.
#Funnel plot also highlights this.

mr_post <- remove_outlier(mr_steiger_df, rucker) #outliers detected!

saveRDS(mr_post, "output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/gsatadj/doctor/instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure                    method nsnp          b        se        pval
# 1    exposure    outcome outcome exposure                  MR Egger   18 -0.3971576 0.3326485 0.249908589
# 2    exposure    outcome outcome exposure           Weighted median   18 -0.2076181 0.1500709 0.166521818
# 3    exposure    outcome outcome exposure Inverse variance weighted   18 -0.2709330 0.1028075 0.008405264
# 4    exposure    outcome outcome exposure               Simple mode   18 -0.2603874 0.2533654 0.318485659
# 5    exposure    outcome outcome exposure             Weighted mode   18 -0.1817248 0.2346358 0.449275263

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method   Estimate         SE      CI_low     CI_upp         P
# 1  Egger fixed effects 0.01183864 0.02804747 -0.04313338 0.06681067 0.6729570
# 2 Egger random effects 0.01183864 0.02959024 -0.04615716 0.06983445 0.3445466
# 
# [[1]]$Q
# Method          Q df         P
# 1   Q_ivw 17.9867591 17 0.3896767
# 2 Q_egger 17.8085965 16 0.3352237
# 3  Q_diff  0.1781626  1 0.6729570
# 
# [[1]]$res
# [1] "A"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_post) #they disagree, but still

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0.0105 [0.0000; 0.2573]; tau = 0.1026 [0.0000; 0.5073]
# I^2 = 0.0% [0.0%; 50.0%]; H = 1.00 [1.00; 1.41]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 16.92   17  0.4597

#That is quite good:

mr_post$labels <- NA
mr_post$units.exposure <- "SD"
mr_post$units.outcome <- "LOR"
 
plotio=mr_plots(mr_post)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/doctor/plots_after_outlier_extraction.svg", height=16, width=16, units="in")
 
saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/doctor/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/doctor/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/doctor/strict_het_test_after_outlier_extraction.RDS")

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

dir.create("output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/broad")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/broad/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp          b         se        pval
# 1    exposure    outcome outcome exposure                  MR Egger   19 -0.1291731 0.27803844 0.648121744
# 2    exposure    outcome outcome exposure           Weighted median   19 -0.2275735 0.10565513 0.031245858
# 3    exposure    outcome outcome exposure Inverse variance weighted   19 -0.2700120 0.08548618 0.001585585
# 4    exposure    outcome outcome exposure               Simple mode   19 -0.2853009 0.21289070 0.196876121
# 5    exposure    outcome outcome exposure             Weighted mode   19 -0.1456507 0.20338986 0.483107836

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method    Estimate         SE      CI_low     CI_upp         P
# 1  Egger fixed effects -0.01332758 0.01906667 -0.05069757 0.02404240 0.4845526
# 2 Egger random effects -0.01332758 0.02498215 -0.06229170 0.03563653 0.7031508
# 
# [[1]]$Q
# Method          Q df          P
# 1   Q_ivw 29.6735396 18 0.04074387
# 2 Q_egger 29.1849400 17 0.03285842
# 3  Q_diff  0.4885997  1 0.48455262
# 
# [[1]]$res
# [1] "B"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_steiger_df) #they disagree, but still

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0.0478 [0.0000; 0.2077]; tau = 0.2185 [0.0000; 0.4557]
# I^2 = 33.6% [0.0%; 61.9%]; H = 1.23 [1.00; 1.62]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 27.12   18  0.0768

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/broad/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/broad/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/broad/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/broad/strict_het_test_before_outlier_extraction.RDS")

#Asymmetry and heterogeneity introduced!

mr_post <- remove_outlier(mr_steiger_df, rucker)
 
saveRDS(mr_post, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/broad/instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure                    method nsnp          b        se         pval
# 1    exposure    outcome outcome exposure                  MR Egger   15 -0.3663006 0.2254549 0.1282090223
# 2    exposure    outcome outcome exposure           Weighted median   15 -0.2364411 0.1027868 0.0214309330
# 3    exposure    outcome outcome exposure Inverse variance weighted   15 -0.2692453 0.0772035 0.0004876046
# 4    exposure    outcome outcome exposure               Simple mode   15 -0.2040899 0.2122624 0.3526153800
# 5    exposure    outcome outcome exposure             Weighted mode   15 -0.1395752 0.1578749 0.3915792252

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method    Estimate         SE      CI_low     CI_upp         P
# 1  Egger fixed effects 0.009309424 0.01752058 -0.02503028 0.04364913 0.5951816
# 2 Egger random effects 0.009309424 0.01752058 -0.02503028 0.04364913 0.2975908
# 
# [[1]]$Q
# Method         Q df         P
# 1   Q_ivw 9.8767194 14 0.7711400
# 2 Q_egger 9.6667832 13 0.7209291
# 3  Q_diff 0.2099362  1 0.6468174
# 
# [[1]]$res
# [1] "A"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_post) #they disagree, but still

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0 [0.0000; 0.0615]; tau = 0 [0.0000; 0.2479]
# I^2 = 0.0% [0.0%; 53.6%]; H = 1.00 [1.00; 1.47]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 9.45   14  0.8009

#That is quite good:

mr_post$labels <- NA
mr_post$units.exposure <- "SD"
mr_post$units.outcome <- "LOR"
 
plotio=mr_plots(mr_post)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/broad/plots_after_outlier_extraction.svg", height=16, width=16, units="in")
 
saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/broad/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/broad/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/broad/strict_het_test_after_outlier_extraction.RDS")

################################################################
#Let's first analyse the data for the whole set for gsatadj - PCOS# consortium definition
################################################################

#STEP 1: let's use only those that have F>10

mr_df <- liu_pcos_consortium[which(((liu_pcos_consortium$beta.exposure^2)/(liu_pcos_consortium$se.exposure^2)) > 10),] #19

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

dir.create("output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/consortium")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/consortium/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp          b         se         pval
# 1    exposure    outcome outcome exposure                  MR Egger   19 -0.1901460 0.11326508 0.1114819206
# 2    exposure    outcome outcome exposure           Weighted median   19 -0.1299165 0.03489492 0.0001968088
# 3    exposure    outcome outcome exposure Inverse variance weighted   19 -0.0783891 0.03527553 0.0262701512
# 4    exposure    outcome outcome exposure               Simple mode   19 -0.1447344 0.05947546 0.0255997290
# 5    exposure    outcome outcome exposure             Weighted mode   19 -0.1385567 0.05551984 0.0225117993

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method   Estimate          SE       CI_low     CI_upp         P
# 1  Egger fixed effects 0.01052965 0.006958159 -0.003108091 0.02416739 0.1302083
# 2 Egger random effects 0.01052965 0.010143326 -0.009350903 0.03041020 0.1496149
# 
# [[1]]$Q
# Method        Q df           P
# 1   Q_ivw 38.41611 18 0.003410828
# 2 Q_egger 36.12609 17 0.004412031
# 3  Q_diff  2.29002  1 0.130208271
# 
# [[1]]$res
# [1] "B"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_steiger_df) #they disagree, but still

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0.0105 [0.0011; 0.0451]; tau = 0.1026 [0.0329; 0.2124]
# I^2 = 47.6% [10.6%; 69.3%]; H = 1.38 [1.06; 1.80]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 34.33   18  0.0115

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/consortium/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/consortium/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/consortium/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/consortium/strict_het_test_before_outlier_extraction.RDS")

#There is pleiotropy.

#STEP 5: removing outliers!

mr_post <- remove_outlier(mr_steiger_df, rucker) #but not significant outliers!!
  
saveRDS(mr_post, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/consortium/instruments_after_outlier_removal_df.RDS")
  
#STEP 6: let's run the results for real now:
  
mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure                    method nsnp           b         se         pval
# 1    exposure    outcome outcome exposure                  MR Egger   17 -0.14328950 0.08244983 0.1027083675
# 2    exposure    outcome outcome exposure           Weighted median   17 -0.13485038 0.03563801 0.0001543949
# 3    exposure    outcome outcome exposure Inverse variance weighted   17 -0.08847733 0.02558887 0.0005449109
# 4    exposure    outcome outcome exposure               Simple mode   17 -0.14458503 0.06459219 0.0397612828
# 5    exposure    outcome outcome exposure             Weighted mode   17 -0.14160290 0.05607439 0.0224993239

# #Let's check the plot:
  
rucker <- TwoSampleMR::mr_rucker(mr_post)

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

strict_het <- compute_strict_het(mr_post) #they disagree, but still

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 < 0.0001 [0.0000; 0.0167]; tau = 0.0017 [0.0000; 0.1293]
# I^2 = 0.0% [0.0%; 51.1%]; H = 1.00 [1.00; 1.43]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 14.60   16  0.5542
#   
#That is quite good:
  
mr_post$labels <- NA
mr_post$units.exposure <- "SD"
mr_post$units.outcome <- "LOR"
  
plotio=mr_plots(mr_post)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/consortium/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/consortium/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/consortium/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/consortium/strict_het_test_after_outlier_extraction.RDS")

#I think we might have now more asymmetry, despite heterogeneity not popping up...

mr_post_2 <- remove_outlier(mr_post, rucker) #but not significant outliers!! - we are fine

###############################################################
#Let's first analyse the data for the whole set for gsatadj - PCOS# age adjusted
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

dir.create("output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/age_adjusted")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/age_adjusted/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp           b         se      pval
# 1    exposure    outcome outcome exposure                  MR Egger   18 -0.05492702 0.32435186 0.8676486
# 2    exposure    outcome outcome exposure           Weighted median   18  0.03300545 0.12208541 0.7868932
# 3    exposure    outcome outcome exposure Inverse variance weighted   18 -0.07097384 0.09525362 0.4562089
# 4    exposure    outcome outcome exposure               Simple mode   18  0.06728753 0.20783393 0.7500716
# 5    exposure    outcome outcome exposure             Weighted mode   18  0.07358572 0.18990852 0.7032094

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method     Estimate         SE      CI_low     CI_upp         P
# 1  Egger fixed effects -0.001522678 0.02468546 -0.04990529 0.04685994 0.9508151
# 2 Egger random effects -0.001522678 0.02933388 -0.05901603 0.05597068 0.5206992
# 
# [[1]]$Q
# Method            Q df         P
# 1   Q_ivw 22.596944972 17 0.1628378
# 2 Q_egger 22.593140154 16 0.1250512
# 3  Q_diff  0.003804817  1 0.9508151
# 
# [[1]]$res
# [1] "A"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_steiger_df) #they disagree, but still

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0.0200 [0.0000; 0.2992]; tau = 0.1413 [0.0000; 0.5470]
# I^2 = 18.7% [0.0%; 53.6%]; H = 1.11 [1.00; 1.47]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 20.91   17  0.2304

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/age_adjusted//plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/age_adjusted//mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/age_adjusted//sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/age_adjusted//strict_het_test_before_outlier_extraction.RDS")

#We have a weird outlier! It does not affect things much, but let's see...

mr_post <- remove_outlier(mr_steiger_df, rucker) #no significant outliers. 

saveRDS(mr_post, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/age_adjusted//instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure                    method nsnp             b         se      pval
# 1    exposure    outcome outcome exposure                  MR Egger   17 -0.2653283742 0.28216070 0.3619384
# 2    exposure    outcome outcome exposure           Weighted median   17 -0.0008582422 0.11885157 0.9942384
# 3    exposure    outcome outcome exposure Inverse variance weighted   17 -0.1159567239 0.08407172 0.1678147
# 4    exposure    outcome outcome exposure               Simple mode   17  0.0640096156 0.22443312 0.7791454
# 5    exposure    outcome outcome exposure             Weighted mode   17  0.0745336637 0.21030649 0.7276631

# #Let's check the plot:
  
rucker <- TwoSampleMR::mr_rucker(mr_post)

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

strict_het <- compute_strict_het(mr_post) 

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0 [0.0000; 0.1319]; tau = 0 [0.0000; 0.3631]
# I^2 = 0.0% [0.0%; 51.1%]; H = 1.00 [1.00; 1.43]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 13.69   16  0.6221

#That is quite good:

mr_post$labels <- NA
mr_post$units.exposure <- "SD"
mr_post$units.outcome <- "LOR"

plotio=mr_plots(mr_post)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/age_adjusted/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/age_adjusted/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/age_adjusted/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/age_adjusted/strict_het_test_after_outlier_extraction.RDS")

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

dir.create("output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/age_and_bmi_adjusted")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/age_and_bmi_adjusted/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp            b        se      pval
# 1    exposure    outcome outcome exposure                  MR Egger   18 -0.030237400 0.3606336 0.9342196
# 2    exposure    outcome outcome exposure           Weighted median   18  0.028934491 0.1449185 0.8417462
# 3    exposure    outcome outcome exposure Inverse variance weighted   18 -0.020172318 0.1056217 0.8485362
# 4    exposure    outcome outcome exposure               Simple mode   18  0.006336669 0.2591146 0.9807744
# 5    exposure    outcome outcome exposure             Weighted mode   18  0.012592872 0.2457655 0.9597320

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method     Estimate         SE      CI_low     CI_upp         P
# 1  Egger fixed effects 0.0009554199 0.02910966 -0.05609846 0.05800930 0.9738170
# 2 Egger random effects 0.0009554199 0.03263571 -0.06300939 0.06492023 0.4883225
# 
# [[1]]$Q
# Method            Q df         P
# 1   Q_ivw 20.111995239 17 0.2685410
# 2 Q_egger 20.110917994 16 0.2152664
# 3  Q_diff  0.001077245  1 0.9738170
# 
# [[1]]$res
# [1] "A"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_steiger_df) #they disagree, but still

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0.0466 [0.0000; 0.2447]; tau = 0.2159 [0.0000; 0.4947]
# I^2 = 11.3% [0.0%; 47.6%]; H = 1.06 [1.00; 1.38]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 19.16   17  0.3192

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/age_and_bmi_adjusted/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/age_and_bmi_adjusted/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/age_and_bmi_adjusted/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/age_and_bmi_adjusted/strict_het_test_before_outlier_extraction.RDS")

#We have a weird outlier

mr_post <- remove_outlier(mr_steiger_df, rucker) #no significant outliers. 

saveRDS(mr_post, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/age_and_bmi_adjusted//instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure                    method nsnp           b        se      pval
# 1    exposure    outcome outcome exposure                  MR Egger   16 0.111293664 0.3578487 0.7603774
# 2    exposure    outcome outcome exposure           Weighted median   16 0.041743732 0.1442462 0.7722813
# 3    exposure    outcome outcome exposure Inverse variance weighted   16 0.039357617 0.1063927 0.7114364
# 4    exposure    outcome outcome exposure               Simple mode   16 0.005694792 0.2304311 0.9806092
# 5    exposure    outcome outcome exposure             Weighted mode   16 0.001079627 0.2152593 0.9960643

# #Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method     Estimate         SE      CI_low     CI_upp         P
# 1  Egger fixed effects -0.006520702 0.02523305 -0.05597658 0.04293518 0.7960835
# 2 Egger random effects -0.006520702 0.02523305 -0.05597658 0.04293518 0.6019582
# 
# [[1]]$Q
# Method          Q df         P
# 1   Q_ivw 9.33754856 15 0.8592171
# 2 Q_egger 9.29321964 14 0.8118207
# 3  Q_diff 0.04432892  1 0.8332428
# 
# [[1]]$res
# [1] "A"

strict_het <- compute_strict_het(mr_post) 

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0 [0.0000; 0.0788]; tau = 0 [0.0000; 0.2806]
# I^2 = 0.0% [0.0%; 52.3%]; H = 1.00 [1.00; 1.45]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 9.05   15  0.8747

#That is quite good:

mr_post$labels <- NA
mr_post$units.exposure <- "SD"
mr_post$units.outcome <- "LOR"

plotio=mr_plots(mr_post)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/age_and_bmi_adjusted/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/age_and_bmi_adjusted/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/age_and_bmi_adjusted/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/age_and_bmi_adjusted/strict_het_test_after_outlier_extraction.RDS")

###############################################################
#Let's first analyse the data for the whole set for bmi - PCOS# Venkatesh
###############################################################

#STEP 1: let's use only those that have F>10

mr_df <- liu_pcos_venkatesh[which(((liu_pcos_venkatesh$beta.exposure^2)/(liu_pcos_venkatesh$se.exposure^2)) > 10),] #19

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

dir.create("output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/venkatesh")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/venkatesh/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp          b         se       pval
# 1    exposure    outcome outcome exposure                  MR Egger   21 -0.1939456 0.24968549 0.44686386
# 2    exposure    outcome outcome exposure           Weighted median   21 -0.1253383 0.11081455 0.25802824
# 3    exposure    outcome outcome exposure Inverse variance weighted   21 -0.1477136 0.08233897 0.07281813
# 4    exposure    outcome outcome exposure               Simple mode   21 -0.0481229 0.19067363 0.80332017
# 5    exposure    outcome outcome exposure             Weighted mode   21 -0.1273029 0.15432850 0.41916827

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method    Estimate         SE      CI_low    CI_upp         P
# 1  Egger fixed effects 0.004565305 0.02144085 -0.03745799 0.0465886 0.8313850
# 2 Egger random effects 0.004565305 0.02144085 -0.03745799 0.0465886 0.4156925
# 
# [[1]]$Q
# Method           Q df         P
# 1   Q_ivw 16.15966609 20 0.7066707
# 2 Q_egger 16.12119812 19 0.6491768
# 3  Q_diff  0.03846797  1 0.8445064
# 
# [[1]]$res
# [1] "A"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_steiger_df) #they disagree, but still

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0 [0.0000; 0.1180]; tau = 0 [0.0000; 0.3435]
# I^2 = 0.0% [0.0%; 47.0%]; H = 1.00 [1.00; 1.37]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 15.37   20  0.7547

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/venkatesh/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/venkatesh/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/venkatesh/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/venkatesh/strict_het_test_before_outlier_extraction.RDS")

#Asymmetry! And some outliers!

mr_post <- remove_outlier(mr_steiger_df, rucker) #but no significant outliers detected...

#saveRDS(mr_post, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/venkatesh/instruments_after_outlier_removal_df.RDS")

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
# ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/venkatesh/plots_after_outlier_extraction.svg", height=16, width=16, units="in")
# 
# saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/venkatesh/mr_res_after_outlier_extraction.RDS")
# saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/venkatesh/sensitivity_test_after_outlier_extraction.RDS")
# saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/gsatadj/venkatesh/strict_het_test_after_outlier_extraction.RDS")