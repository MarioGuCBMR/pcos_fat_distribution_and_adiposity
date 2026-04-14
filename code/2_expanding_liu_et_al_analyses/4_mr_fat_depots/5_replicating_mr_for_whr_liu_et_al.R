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

liu_pcos_doctor <- fread("output/2_replicating_liu_et_al/2_mr/1_expo_outcome_df/whr_liu_PCOS.txt")
liu_pcos_broad <- fread("output/2_replicating_liu_et_al/2_mr/1_expo_outcome_df/whr_liu_PCOS (Broad).txt")
liu_pcos_consortium <- fread("output/2_replicating_liu_et_al/2_mr/1_expo_outcome_df/whr_liu_PCOS (Consortium).txt")
liu_pcos_meta_analysis <- fread("output/2_replicating_liu_et_al/2_mr/1_expo_outcome_df/whr_liu_PCOS (Day et al).txt")
liu_pcos_adj_age <- fread("output/2_replicating_liu_et_al/2_mr/1_expo_outcome_df/whr_liu_PCOS (adj age).txt")
liu_pcos_adj_age_bmi <- fread("output/2_replicating_liu_et_al/2_mr/1_expo_outcome_df/whr_liu_PCOS (adj age+BMI).txt")
liu_pcos_venkatesh <- fread("output/2_replicating_liu_et_al/2_mr/1_expo_outcome_df/whr_liu_PCOS (Venkatesh et al).txt")

liu_amenorrhea<- fread("output/2_replicating_liu_et_al/2_mr/1_expo_outcome_df/whr_liu_Amenorrhea.txt")
liu_oligomenorrhea<- fread("output/2_replicating_liu_et_al/2_mr/1_expo_outcome_df/whr_liu_Oligomenorrhea.txt")
liu_menorrhagia<- fread("output/2_replicating_liu_et_al/2_mr/1_expo_outcome_df/whr_liu_Menorrhagia.txt")
liu_anovulation<- fread("output/2_replicating_liu_et_al/2_mr/1_expo_outcome_df/whr_liu_Anovulation.txt")
liu_hirsutism<- fread("output/2_replicating_liu_et_al/2_mr/1_expo_outcome_df/whr_liu_Hirsutism.txt")

##############################################################
#STEP 2: let's first rerun the analyses for the meta-analysis#
##############################################################

mr_df <- liu_pcos_meta_analysis[which(((as.numeric(liu_pcos_meta_analysis$beta.exposure)^2)/(as.numeric(liu_pcos_meta_analysis$se.exposure)^2)) > 10),] #250/250

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

mr_df$id.exposure <- "WHR"
mr_df$id.outcome <- "PCOS (meta-analysis)"
mr_df$exposure <- "WHR"
mr_df$outcome <- "PCOS (meta-analysis)"
mr_df$units.exposure <- "SD"
mr_df$units.outcome <- "LOR"

mr_steiger_df <- TwoSampleMR::steiger_filtering(mr_df)
mr_steiger_df <- mr_steiger_df[which(mr_steiger_df$steiger_dir == TRUE),] #259
mr_steiger_df$mr_keep <- TRUE

#Let's save the data:

dir.create("output/2_replicating_liu_et_al/2_mr/2_mr/")
dir.create("output/2_replicating_liu_et_al/2_mr/2_mr/whr/")
dir.create("output/2_replicating_liu_et_al/2_mr/2_mr/whr/meta_analysis/")
saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al/2_mr/2_mr/whr/meta_analysis/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure           id.outcome              outcome exposure                    method nsnp           b        se      pval
# 1         WHR PCOS (meta-analysis) PCOS (meta-analysis)      WHR                  MR Egger  186  0.15463856 0.2772434 0.5776781
# 2         WHR PCOS (meta-analysis) PCOS (meta-analysis)      WHR           Weighted median  186 -0.16280355 0.1845535 0.3776967
# 3         WHR PCOS (meta-analysis) PCOS (meta-analysis)      WHR Inverse variance weighted  186  0.08130322 0.1143126 0.4769380
# 4         WHR PCOS (meta-analysis) PCOS (meta-analysis)      WHR               Simple mode  186 -0.19636045 0.3931405 0.6180444
# 5         WHR PCOS (meta-analysis) PCOS (meta-analysis)      WHR             Weighted mode  186 -0.19636045 0.2653708 0.4602698

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method     Estimate          SE      CI_low     CI_upp         P
# 1  Egger fixed effects -0.001900031 0.006314635 -0.01427649 0.01047643 0.7634959
# 2 Egger random effects -0.001900031 0.006540706 -0.01471958 0.01091952 0.6142806
# 
# [[1]]$Q
# Method            Q  df         P
# 1   Q_ivw 197.50120080 185 0.2513348
# 2 Q_egger 197.41066402 184 0.2365828
# 3  Q_diff   0.09053678   1 0.7634959
# 
# [[1]]$res
# [1] "A"

strict_het <- compute_strict_het(mr_steiger_df)

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al/2_mr/2_mr/whr/meta_analysis/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al/2_mr/2_mr/whr/meta_analysis/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al/2_mr/2_mr/whr/meta_analysis//sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al/2_mr/2_mr/whr/meta_analysis//strict_het_test_before_outlier_extraction.RDS")

#Sensitivity plots show some outliers, though they are not strong enough to affect the results according to the rest of sensitivity analyses.
#Let's run the pipeline just in case.

mr_post <- remove_outlier(mr_steiger_df, rucker)

saveRDS(mr_post, "output/2_replicating_liu_et_al/2_mr/2_mr/whr/meta_analysis/instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure           id.outcome              outcome exposure                    method nsnp           b        se      pval
# 1         WHR PCOS (meta-analysis) PCOS (meta-analysis)      WHR                  MR Egger  176 -0.08816436 0.2717425 0.7459940
# 2         WHR PCOS (meta-analysis) PCOS (meta-analysis)      WHR           Weighted median  176 -0.22222222 0.1868190 0.2342409
# 3         WHR PCOS (meta-analysis) PCOS (meta-analysis)      WHR Inverse variance weighted  176 -0.02526873 0.1137197 0.8241566
# 4         WHR PCOS (meta-analysis) PCOS (meta-analysis)      WHR               Simple mode  176 -0.21395949 0.4165228 0.6081236
# 5         WHR PCOS (meta-analysis) PCOS (meta-analysis)      WHR             Weighted mode  176 -0.21395949 0.2672403 0.4244339

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method    Estimate          SE       CI_low     CI_upp         P
# 1  Egger fixed effects 0.001627436 0.005490137 -0.009133035 0.01238791 0.7669024
# 2 Egger random effects 0.001627436 0.005490137 -0.009133035 0.01238791 0.3834512
# 
# [[1]]$Q
# Method            Q  df         P
# 1   Q_ivw 128.66670863 175 0.9965356
# 2 Q_egger 128.60176460 174 0.9959537
# 3  Q_diff   0.06494403   1 0.7988458
# 
# [[1]]$res
#[1] "A"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_post) #they disagree, but still

#That is quite good:

mr_post$labels <- NA
mr_post$units.exposure <- "SD"
mr_post$units.outcome <- "LOR"

plotio=mr_plots(mr_post)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/whr/meta_analysis/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/meta_analysis/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/meta_analysis/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/meta_analysis/strict_het_test_after_outlier_extraction.RDS")

###############################################################
#Let's first analyse the data for the whole set for bmi - PCOS#
###############################################################

#STEP 1: let's use only those that have F>10

mr_df <- liu_pcos_doctor[which(((liu_pcos_doctor$beta.exposure^2)/(liu_pcos_doctor$se.exposure^2)) > 10),] 

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

dir.create("output/2_replicating_liu_et_al//2_mr/2_mr/whr/")
dir.create("output/2_replicating_liu_et_al//2_mr/2_mr/whr/doctor")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/doctor/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp         b        se         pval
# 1    exposure    outcome outcome exposure                  MR Egger  194 0.8586614 0.2583974 1.066248e-03
# 2    exposure    outcome outcome exposure           Weighted median  194 0.6251470 0.1789736 4.777038e-04
# 3    exposure    outcome outcome exposure Inverse variance weighted  194 0.5267479 0.1107185 1.959665e-06
# 4    exposure    outcome outcome exposure               Simple mode  194 0.6606958 0.4260416 1.225932e-01
# 5    exposure    outcome outcome exposure             Weighted mode  194 0.7789321 0.2529562 2.377232e-03

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method     Estimate          SE      CI_low      CI_upp         P
# 1  Egger fixed effects -0.008661279 0.005849150 -0.02012540 0.002802845 0.1386664
# 2 Egger random effects -0.008661279 0.006096126 -0.02060947 0.003286907 0.9223103
# 
# [[1]]$Q
# Method          Q  df         P
# 1   Q_ivw 210.749071 193 0.1812016
# 2 Q_egger 208.556374 192 0.1961240
# 3  Q_diff   2.192697   1 0.1386664
# 
# [[1]]$res
# [1] "A"

strict_het <- compute_strict_het(mr_steiger_df)

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/whr/doctor/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/doctor/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/doctor/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/doctor/strict_het_test_before_outlier_extraction.RDS")

#A bit of asymmetry due to one variant that is highlighted in the funnel plot. I was close to not remove it, but since both loo and funnel show it, I will do it in this case.

mr_post <- remove_outlier(mr_steiger_df, rucker)

saveRDS(mr_post, "output/2_replicating_liu_et_al/2_mr/2_mr/whr/doctor/instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure                    method nsnp         b        se         pval
# 1    exposure    outcome outcome exposure                  MR Egger  182 0.8733452 0.2539906 7.267527e-04
# 2    exposure    outcome outcome exposure           Weighted median  182 0.6449801 0.1761695 2.511026e-04
# 3    exposure    outcome outcome exposure Inverse variance weighted  182 0.5803432 0.1091846 1.065132e-07
# 4    exposure    outcome outcome exposure               Simple mode  182 0.6500320 0.3973327 1.035796e-01
# 5    exposure    outcome outcome exposure             Weighted mode  182 0.7834155 0.2389288 1.249887e-03
#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method     Estimate          SE      CI_low      CI_upp         P
# 1  Egger fixed effects -0.007685895 0.005277408 -0.01802943 0.002657635 0.1452885
# 2 Egger random effects -0.007685895 0.005277408 -0.01802943 0.002657635 0.9273558
# 
# [[1]]$Q
# Method          Q  df         P
# 1   Q_ivw 140.168657 181 0.9891886
# 2 Q_egger 138.536213 180 0.9904659
# 3  Q_diff   1.632444   1 0.2013653
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/whr/doctor/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/doctor/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/doctor/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/doctor/strict_het_test_after_outlier_extraction.RDS")

###############################################################
#Let's first analyse the data for the whole set for bmi - PCOS# broad definition - anovulation
###############################################################

#STEP 1: let's use only those that have F>10

mr_df <- liu_pcos_broad[which(((liu_pcos_broad$beta.exposure^2)/(liu_pcos_broad$se.exposure^2)) > 10),] 

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

dir.create("output/2_replicating_liu_et_al//2_mr/2_mr/whr/broad")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/broad/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp           b         se         pval
# 1    exposure    outcome outcome exposure                  MR Egger  194  0.77677140 0.19530004 9.867602e-05
# 2    exposure    outcome outcome exposure           Weighted median  194  0.44873326 0.13009311 5.619921e-04
# 3    exposure    outcome outcome exposure Inverse variance weighted  194  0.31789721 0.08425401 1.612440e-04
# 4    exposure    outcome outcome exposure               Simple mode  194 -0.06459667 0.33961288 8.493470e-01
# 5    exposure    outcome outcome exposure             Weighted mode  194  0.58139940 0.18860739 2.352457e-03

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method    Estimate          SE      CI_low       CI_upp           P
# 1  Egger fixed effects -0.01194524 0.004015236 -0.01981495 -0.004075520 0.002930098
# 2 Egger random effects -0.01194524 0.004601661 -0.02096433 -0.002926148 0.995282185
# 
# [[1]]$Q
# Method          Q  df            P
# 1   Q_ivw 261.029155 193 0.0007990933
# 2 Q_egger 252.178662 192 0.0023153101
# 3  Q_diff   8.850493   1 0.0029300979
# 
# [[1]]$res
# [1] "D"

#Pleitropy!

strict_het <- compute_strict_het(mr_steiger_df) #they disagree, but still

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/whr/broad/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/broad/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/broad/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/broad/strict_het_test_before_outlier_extraction.RDS")

#Pleitropy detected

mr_post <- remove_outlier(mr_steiger_df, rucker)

saveRDS(mr_post, "output/2_replicating_liu_et_al/2_mr/2_mr/whr/broad/instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure                    method nsnp          b         se         pval
# 1    exposure    outcome outcome exposure                  MR Egger  179  0.7408406 0.17292253 3.004147e-05
# 2    exposure    outcome outcome exposure           Weighted median  179  0.4488271 0.13032243 5.732302e-04
# 3    exposure    outcome outcome exposure Inverse variance weighted  179  0.3117668 0.07460186 2.926960e-05
# 4    exposure    outcome outcome exposure               Simple mode  179 -0.1743263 0.35607888 6.250395e-01
# 5    exposure    outcome outcome exposure             Weighted mode  179  0.5903759 0.20476042 4.421270e-03

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method    Estimate          SE      CI_low       CI_upp           P
# 1  Egger fixed effects -0.01130263 0.003742004 -0.01863682 -0.003968435 0.002523789
# 2 Egger random effects -0.01130263 0.003742004 -0.01863682 -0.003968435 0.998738106
# 
# [[1]]$Q
# Method          Q  df           P
# 1   Q_ivw 154.330411 178 0.899568550
# 2 Q_egger 146.765545 177 0.952819408
# 3  Q_diff   7.564866   1 0.005951709
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/whr/broad/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/broad/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/broad/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/broad/strict_het_test_after_outlier_extraction.RDS")

#Still biased, though due to some heterogeneity.
#I don't want to cherry-pick, so I am going to try reporting it.

###############################################################
#Let's first analyse the data for the whole set for bmi - PCOS# consortium definition - anovulation
###############################################################

#STEP 1: let's use only those that have F>10

mr_df <- liu_pcos_consortium[which(((liu_pcos_consortium$beta.exposure^2)/(liu_pcos_consortium$se.exposure^2)) > 10),] #250

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

dir.create("output/2_replicating_liu_et_al//2_mr/2_mr/whr/consortium")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/consortium/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp         b         se         pval
# 1    exposure    outcome outcome exposure                  MR Egger  194 0.2487087 0.07998893 2.160731e-03
# 2    exposure    outcome outcome exposure           Weighted median  194 0.1695541 0.05220394 1.162528e-03
# 3    exposure    outcome outcome exposure Inverse variance weighted  194 0.1425803 0.03367269 2.292629e-05
# 4    exposure    outcome outcome exposure               Simple mode  194 0.1101866 0.11287037 3.301764e-01
# 5    exposure    outcome outcome exposure             Weighted mode  194 0.1740928 0.06164578 5.239355e-03

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method     Estimate          SE       CI_low       CI_upp          P
# 1  Egger fixed effects -0.002744691 0.001468956 -0.005623792 0.0001344098 0.06169746
# 2 Egger random effects -0.002744691 0.001877627 -0.006424772 0.0009353902 0.92810023
# 
# [[1]]$Q
# Method          Q  df            P
# 1   Q_ivw 317.182329 193 4.339317e-08
# 2 Q_egger 313.691172 192 6.805327e-08
# 3  Q_diff   3.491157   1 6.169746e-02
# 
# [[1]]$res
# [1] "B"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_steiger_df) #they disagree, but still

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/whr/consortium/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/consortium/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/consortium/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/consortium/strict_het_test_before_outlier_extraction.RDS")

#Very small heterogeneity, overall. Small outlier in loo and funnel, so let's check it out.

#STEP 5: removing outliers!

mr_post <- remove_outlier(mr_steiger_df, rucker)
  
saveRDS(mr_post, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/consortium/instruments_after_outlier_removal_df.RDS")
  
#STEP 6: let's run the results for real now:
  
mr_res <- TwoSampleMR::mr(mr_post)
  
# id.exposure id.outcome outcome exposure                    method nsnp         b         se         pval
# 1    exposure    outcome outcome exposure                  MR Egger  173 0.2134965 0.06606434 1.476598e-03
# 2    exposure    outcome outcome exposure           Weighted median  173 0.1694391 0.05354817 1.554905e-03
# 3    exposure    outcome outcome exposure Inverse variance weighted  173 0.1379075 0.02784589 7.325688e-07
# 4    exposure    outcome outcome exposure               Simple mode  173 0.1159140 0.10477475 2.701347e-01
# 5    exposure    outcome outcome exposure             Weighted mode  173 0.1787304 0.06082623 3.752241e-03

# #Let's check the plot:
#  
rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method     Estimate          SE       CI_low      CI_upp        P
# 1  Egger fixed effects -0.001959496 0.001518234 -0.004935181 0.001016188 0.196828
# 2 Egger random effects -0.001959496 0.001518234 -0.004935181 0.001016188 0.901586
# 
# [[1]]$Q
# Method          Q  df         P
# 1   Q_ivw 165.016456 172 0.6352047
# 2 Q_egger 163.424495 171 0.6478984
# 3  Q_diff   1.591961   1 0.2070462
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/whr/consortium/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/consortium/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/consortium/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/consortium/strict_het_test_after_outlier_extraction.RDS")

###############################################################
#Let's first analyse the data for the whole set for bmi - PCOS# age adjusted
###############################################################

#STEP 1: let's use only those that have F>10

mr_df <- liu_pcos_adj_age[which(((liu_pcos_adj_age$beta.exposure^2)/(liu_pcos_adj_age$se.exposure^2)) > 10),] #194
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

dir.create("output/2_replicating_liu_et_al//2_mr/2_mr/whr/age_adjusted")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/age_adjusted/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp         b         se         pval
# 1    exposure    outcome outcome exposure                  MR Egger  194 0.2718276 0.22782531 0.2342867128
# 2    exposure    outcome outcome exposure           Weighted median  194 0.3645531 0.14976937 0.0149290118
# 3    exposure    outcome outcome exposure Inverse variance weighted  194 0.3399151 0.09429699 0.0003124787
# 4    exposure    outcome outcome exposure               Simple mode  194 0.3838887 0.34071318 0.2612594059
# 5    exposure    outcome outcome exposure             Weighted mode  194 0.3838887 0.20444964 0.0619349035

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method    Estimate          SE       CI_low     CI_upp         P
# 1  Egger fixed effects 0.001757082 0.004956336 -0.007957158 0.01147132 0.7229550
# 2 Egger random effects 0.001757082 0.005349498 -0.008727740 0.01224191 0.3712829
# 
# [[1]]$Q
# Method          Q  df          P
# 1   Q_ivw 223.794628 193 0.06368651
# 2 Q_egger 223.668949 192 0.05838147
# 3  Q_diff   0.125679   1 0.72295497
# 
# [[1]]$res
# [1] "A"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_steiger_df) #they disagree, but still

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/whr/age_adjusted/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/age_adjusted/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/age_adjusted/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/age_adjusted/strict_het_test_before_outlier_extraction.RDS")

#A single bizarre weak outlier pops up. 
#Removing outliers will make it go.

mr_post <- remove_outlier(mr_steiger_df, rucker)

saveRDS(mr_post, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/age_adjusted/instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure                    method nsnp         b         se         pval
# 1    exposure    outcome outcome exposure                  MR Egger  178 0.2957167 0.21511194 1.709689e-01
# 2    exposure    outcome outcome exposure           Weighted median  178 0.3661591 0.14754456 1.307628e-02
# 3    exposure    outcome outcome exposure Inverse variance weighted  178 0.3596927 0.09055834 7.128827e-05
# 4    exposure    outcome outcome exposure               Simple mode  178 0.3768050 0.32426352 2.467864e-01
# 5    exposure    outcome outcome exposure             Weighted mode  178 0.3932193 0.22482605 8.202555e-02

# #Let's check the plot:
#  
rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method   Estimate          SE       CI_low     CI_upp         P
# 1  Egger fixed effects 0.00166314 0.004292205 -0.006749427 0.01007571 0.6984015
# 2 Egger random effects 0.00166314 0.004292205 -0.006749427 0.01007571 0.3492008
# 
# [[1]]$Q
# Method           Q  df         P
# 1   Q_ivw 126.1280667 177 0.9985637
# 2 Q_egger 126.0205625 176 0.9983157
# 3  Q_diff   0.1075042   1 0.7430037
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/whr/age_adjusted/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/age_adjusted/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/age_adjusted/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/age_adjusted/strict_het_test_after_outlier_extraction.RDS")

###############################################################
#Let's first analyse the data for the whole set for bmi - PCOS# age and BMI adjusted
###############################################################

#STEP 1: let's use only those that have F>10

mr_df <- liu_pcos_adj_age_bmi[which(((liu_pcos_adj_age_bmi$beta.exposure^2)/(liu_pcos_adj_age_bmi$se.exposure^2)) > 10),] #250
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

dir.create("output/2_replicating_liu_et_al//2_mr/2_mr/whr/age_and_bmi_adjusted")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/age_and_bmi_adjusted/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp         b        se       pval
# 1    exposure    outcome outcome exposure                  MR Egger  194 0.1222136 0.2747535 0.65695669
# 2    exposure    outcome outcome exposure           Weighted median  194 0.2649234 0.1918063 0.16721657
# 3    exposure    outcome outcome exposure Inverse variance weighted  194 0.2241805 0.1135356 0.04832042
# 4    exposure    outcome outcome exposure               Simple mode  194 0.2644583 0.4098644 0.51954200
# 5    exposure    outcome outcome exposure             Weighted mode  194 0.2277979 0.2516759 0.36652814

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method    Estimate          SE       CI_low     CI_upp         P
# 1  Egger fixed effects 0.002629326 0.005829213 -0.008795721 0.01405437 0.6519461
# 2 Egger random effects 0.002629326 0.006448751 -0.010009994 0.01526865 0.3417373
# 
# [[1]]$Q
# Method           Q  df          P
# 1   Q_ivw 235.1844260 193 0.02060548
# 2 Q_egger 234.9809706 192 0.01868691
# 3  Q_diff   0.2034554   1 0.65194612
# 
# [[1]]$res
# [1] "B"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_steiger_df) #they disagree, but still

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/whr/age_and_bmi_adjusted/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/age_and_bmi_adjusted/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/age_and_bmi_adjusted/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/age_and_bmi_adjusted/strict_het_test_before_outlier_extraction.RDS")

#again, some small outliers...
#Let's go one more time with feeling

mr_post <- remove_outlier(mr_steiger_df, rucker)

saveRDS(mr_post, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/age_and_bmi_adjusted/instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure                    method nsnp         b        se        pval
# 1    exposure    outcome outcome exposure                  MR Egger  176 0.3151475 0.2610679 0.229013749
# 2    exposure    outcome outcome exposure           Weighted median  176 0.2942148 0.1920616 0.125552788
# 3    exposure    outcome outcome exposure Inverse variance weighted  176 0.3018660 0.1085042 0.005401401
# 4    exposure    outcome outcome exposure               Simple mode  176 0.2998709 0.4407171 0.497140361
# 5    exposure    outcome outcome exposure             Weighted mode  176 0.2571809 0.2614962 0.326720234

# #Let's check the plot:
#  
rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method      Estimate          SE      CI_low     CI_upp         P
# 1  Egger fixed effects -0.0003405564 0.005302258 -0.01073279 0.01005168 0.9487882
# 2 Egger random effects -0.0003405564 0.005302258 -0.01073279 0.01005168 0.5256059
# 
# [[1]]$Q
# Method            Q  df         P
# 1   Q_ivw 1.319614e+02 175 0.9935302
# 2 Q_egger 1.319582e+02 174 0.9924406
# 3  Q_diff 3.128553e-03   1 0.9553948
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/whr/age_and_bmi_adjusted/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/age_and_bmi_adjusted/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/age_and_bmi_adjusted/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/age_and_bmi_adjusted/strict_het_test_after_outlier_extraction.RDS")

###############################################################
#Let's first analyse the data for the whole set for bmi - PCOS# Venkatesh
###############################################################

#STEP 1: let's use only those that have F>10

mr_df <- liu_pcos_venkatesh[which(((liu_pcos_venkatesh$beta.exposure^2)/(liu_pcos_venkatesh$se.exposure^2)) > 10),] #194

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

dir.create("output/2_replicating_liu_et_al//2_mr/2_mr/whr/venkatesh")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/venkatesh/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp         b        se         pval
# 1    exposure    outcome outcome exposure                  MR Egger  196 0.5128677 0.2456733 0.0381398582
# 2    exposure    outcome outcome exposure           Weighted median  196 0.4849042 0.1754168 0.0057045480
# 3    exposure    outcome outcome exposure Inverse variance weighted  196 0.4140557 0.1062803 0.0000978415
# 4    exposure    outcome outcome exposure               Simple mode  196 0.4625879 0.3110727 0.1386105937
# 5    exposure    outcome outcome exposure             Weighted mode  196 0.4267917 0.2112394 0.0447080602

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method     Estimate          SE      CI_low      CI_upp         P
# 1  Egger fixed effects -0.002707809 0.005255069 -0.01300755 0.007591937 0.6063605
# 2 Egger random effects -0.002707809 0.006066880 -0.01459867 0.009183057 0.6723192
# 
# [[1]]$Q
# Method           Q  df           P
# 1   Q_ivw 258.8340760 195 0.001500521
# 2 Q_egger 258.5685671 194 0.001328602
# 3  Q_diff   0.2655089   1 0.606360462
# 
# [[1]]$res
# [1] "B"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_steiger_df) #they disagree, but still

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/whr/venkatesh/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/venkatesh/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/venkatesh/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/venkatesh/strict_het_test_before_outlier_extraction.RDS")

#Quite the heterogeneity!

mr_post <- remove_outlier(mr_steiger_df, rucker)

saveRDS(mr_post, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/venkatesh/instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure                    method nsnp         b         se        pval
# 1    exposure    outcome outcome exposure                  MR Egger  182 0.3850080 0.21755072 0.078463714
# 2    exposure    outcome outcome exposure           Weighted median  182 0.4687462 0.17700812 0.008093063
# 3    exposure    outcome outcome exposure Inverse variance weighted  182 0.3182588 0.09863094 0.001251985
# 4    exposure    outcome outcome exposure               Simple mode  182 0.4537847 0.30790687 0.142279873
# 5    exposure    outcome outcome exposure             Weighted mode  182 0.4537847 0.19299455 0.019784703

# #Let's check the plot:
#  
rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method     Estimate         SE      CI_low      CI_upp         P
# 1  Egger fixed effects -0.001876149 0.00481659 -0.01131649 0.007564194 0.6968929
# 2 Egger random effects -0.001876149 0.00481659 -0.01131649 0.007564194 0.6515535
# 
# [[1]]$Q
# Method           Q  df         P
# 1   Q_ivw 140.6969867 181 0.9882127
# 2 Q_egger 140.5784913 180 0.9866701
# 3  Q_diff   0.1184955   1 0.7306720
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/whr/venkatesh/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/venkatesh/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/venkatesh/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/venkatesh/strict_het_test_after_outlier_extraction.RDS")

######################
#Let's try amenorrhea#
######################

#STEP 1: let's use only those that have F>10

mr_df <- liu_amenorrhea[which(((liu_amenorrhea$beta.exposure^2)/(liu_amenorrhea$se.exposure^2)) > 10),] #278

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

dir.create("output/2_replicating_liu_et_al//2_mr/2_mr/whr/amenorrhea")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/amenorrhea/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp          b         se       pval
# 1    exposure    outcome outcome exposure                  MR Egger  194 0.49474723 0.22104839 0.02635764
# 2    exposure    outcome outcome exposure           Weighted median  194 0.10522096 0.17027758 0.53661631
# 3    exposure    outcome outcome exposure Inverse variance weighted  194 0.01203904 0.09400941 0.89809988
# 4    exposure    outcome outcome exposure               Simple mode  194 0.47459352 0.36984035 0.20094712
# 5    exposure    outcome outcome exposure             Weighted mode  194 0.31869063 0.21720401 0.14393815

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method    Estimate          SE      CI_low       CI_upp          P
# 1  Egger fixed effects -0.01248858 0.005046045 -0.02237864 -0.002598510 0.01332646
# 2 Egger random effects -0.01248858 0.005189794 -0.02266039 -0.002316767 0.99194408
# 
# [[1]]$Q
# Method          Q  df          P
# 1   Q_ivw 209.220249 193 0.20122835
# 2 Q_egger 203.095002 192 0.27776119
# 3  Q_diff   6.125248   1 0.01332646
# 
# [[1]]$res
# [1] "A"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_steiger_df) #they disagree, but still

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/whr/amenorrhea/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/amenorrhea/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/amenorrhea/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/amenorrhea/strict_het_test_before_outlier_extraction.RDS")

#Pleiotropy! Let's remove the outliers
 
mr_post <- remove_outlier(mr_steiger_df, rucker)
 
saveRDS(mr_post, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/amenorrhea/instruments_after_outlier_removal_df.RDS")
 
#STEP 6: let's run the results for real now:
 
mr_res <- TwoSampleMR::mr(mr_post)
 
# id.exposure id.outcome outcome exposure                    method nsnp            b         se      pval
# 1    exposure    outcome outcome exposure                  MR Egger  180  0.347208145 0.22291614 0.1211107
# 2    exposure    outcome outcome exposure           Weighted median  180  0.059938027 0.17200999 0.7274973
# 3    exposure    outcome outcome exposure Inverse variance weighted  180 -0.006302846 0.09408723 0.9465902
# 4    exposure    outcome outcome exposure               Simple mode  180  0.483134024 0.38585786 0.2121657
# 5    exposure    outcome outcome exposure             Weighted mode  180  0.300950228 0.23553442 0.2029980

#Let's check the plot:
  
rucker <- TwoSampleMR::mr_rucker(mr_post)
 
# [[1]]$intercept
# Method     Estimate          SE      CI_low        CI_upp          P
# 1  Egger fixed effects -0.009111712 0.004391308 -0.01771852 -0.0005049066 0.03799185
# 2 Egger random effects -0.009111712 0.004391308 -0.01771852 -0.0005049066 0.98100407
# 
# [[1]]$Q
# Method          Q  df          P
# 1   Q_ivw 129.573427 179 0.99791333
# 2 Q_egger 126.513376 178 0.99870014
# 3  Q_diff   3.060051   1 0.08023916
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/whr/amenorrhea/plots_after_outlier_extraction.svg", height=16, width=16, units="in")
 
saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/amenorrhea/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/amenorrhea/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/amenorrhea/strict_het_test_after_outlier_extraction.RDS")

##########################
#Let's try oligomenorrhea#
##########################

#STEP 1: let's use only those that have F>10

mr_df <- liu_oligomenorrhea[which(((liu_oligomenorrhea$beta.exposure^2)/(liu_oligomenorrhea$se.exposure^2)) > 10),] #276

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

dir.create("output/2_replicating_liu_et_al//2_mr/2_mr/whr/oligomenorrhea")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/oligomenorrhea/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp         b        se        pval
# 1    exposure    outcome outcome exposure                  MR Egger  194 0.7973695 0.3133257 0.011718591
# 2    exposure    outcome outcome exposure           Weighted median  194 0.6032181 0.2096101 0.004004468
# 3    exposure    outcome outcome exposure Inverse variance weighted  194 0.3622943 0.1314876 0.005862943
# 4    exposure    outcome outcome exposure               Simple mode  194 0.7235030 0.5024013 0.151462165
# 5    exposure    outcome outcome exposure             Weighted mode  194 0.7823771 0.3054404 0.011186589

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method    Estimate          SE      CI_low     CI_upp         P
# 1  Egger fixed effects -0.01125042 0.007035623 -0.02503999 0.00253915 0.1098062
# 2 Egger random effects -0.01125042 0.007035623 -0.02503999 0.00253915 0.9450969
# 
# [[1]]$Q
# Method          Q  df         P
# 1   Q_ivw 178.065923 193 0.7723176
# 2 Q_egger 175.725652 192 0.7941539
# 3  Q_diff   2.340271   1 0.1260676
# 
# [[1]]$res
# [1] "A"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_steiger_df) #they disagree, but still

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/whr/oligomenorrhea/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/oligomenorrhea/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/oligomenorrhea/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/oligomenorrhea/strict_het_test_before_outlier_extraction.RDS")

#Very controlled..., let's run it just in case - I can see some stuff...

mr_post <- remove_outlier(mr_steiger_df, rucker) #removing outliers removes variants and makes it pleiotropic!!! We are not doing it.
#mr_post <- remove_outlier(mr_post, rucker) #Run twice
#mr_post <- remove_outlier(mr_post, rucker) #does not run. the best we can do honestly

# saveRDS(mr_post, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/oligomenorrhea/instruments_after_outlier_removal_df.RDS")
# 
# #STEP 6: let's run the results for real now:
# 
# mr_res <- TwoSampleMR::mr(mr_post)
# 
# # id.exposure id.outcome outcome exposure                    method nsnp         b        se        pval
# # 1    exposure    outcome outcome exposure                  MR Egger  183 0.9303201 0.3188638 0.003974389
# # 2    exposure    outcome outcome exposure           Weighted median  183 0.6267944 0.2088986 0.002695614
# # 3    exposure    outcome outcome exposure Inverse variance weighted  183 0.4170915 0.1356865 0.002112553
# # 4    exposure    outcome outcome exposure               Simple mode  183 0.7441177 0.4834139 0.125468525
# # 5    exposure    outcome outcome exposure             Weighted mode  183 0.8156315 0.3287275 0.014002305
# 
# #Let's check the plot:
# #  
# rucker <- TwoSampleMR::mr_rucker(mr_post)
# 
# # [[1]]$intercept
# # Method    Estimate          SE      CI_low        CI_upp          P
# # 1  Egger fixed effects -0.01248619 0.006093084 -0.02442842 -0.0005439683 0.04043862
# # 2 Egger random effects -0.01248619 0.006093084 -0.02442842 -0.0005439683 0.97978069
# # 
# # [[1]]$Q
# # Method          Q  df          P
# # 1   Q_ivw 174.837171 261 0.99998981
# # 2 Q_egger 172.058177 260 0.99999391
# # 3  Q_diff   2.778995   1 0.09550808
# # 
# # [[1]]$res
# # [1] "A"
# 
# # #Here there is heterogeneity and pleitropy
# 
# strict_het <- compute_strict_het(mr_post) #they disagree, but still
# 
# #That is quite good:
# 
# mr_post$labels <- NA
# mr_post$units.exposure <- "SD"
# mr_post$units.outcome <- "LOR"
# 
# plotio=mr_plots(mr_post)
# ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/whr/oligomenorrhea/plots_after_outlier_extraction.svg", height=16, width=16, units="in")
# 
# saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/oligomenorrhea/mr_res_after_outlier_extraction.RDS")
# saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/oligomenorrhea/sensitivity_test_after_outlier_extraction.RDS")
# saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/oligomenorrhea/strict_het_test_after_outlier_extraction.RDS")

#######################
#Let's try menorrhagia#
#######################

#STEP 1: let's use only those that have F>10

mr_df <- liu_menorrhagia[which(((liu_menorrhagia$beta.exposure^2)/(liu_menorrhagia$se.exposure^2)) > 10),] #276

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

dir.create("output/2_replicating_liu_et_al//2_mr/2_mr/whr/menorrhagia")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/menorrhagia/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp           b         se         pval
# 1    exposure    outcome outcome exposure                  MR Egger  194  0.39154347 0.10598001 2.872245e-04
# 2    exposure    outcome outcome exposure           Weighted median  194  0.22722877 0.05707733 6.860616e-05
# 3    exposure    outcome outcome exposure Inverse variance weighted  194  0.18031093 0.04508883 6.360644e-05
# 4    exposure    outcome outcome exposure               Simple mode  194 -0.02108491 0.19121609 9.123118e-01
# 5    exposure    outcome outcome exposure             Weighted mode  194  0.32781843 0.09265327 5.048018e-04

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method     Estimate          SE       CI_low        CI_upp           P
# 1  Egger fixed effects -0.005472767 0.001855173 -0.009108839 -0.0018366952 0.003177698
# 2 Egger random effects -0.005472767 0.002490248 -0.010353563 -0.0005919712 0.986014033
# 
# [[1]]$Q
# Method          Q  df            P
# 1   Q_ivw 354.655932 193 1.177653e-11
# 2 Q_egger 345.953408 192 6.479301e-11
# 3  Q_diff   8.702524   1 3.177698e-03
# 
# [[1]]$res
# [1] "D"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_steiger_df) #they disagree, but still

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/whr/menorrhagia/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/menorrhagia/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/menorrhagia/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/menorrhagia/strict_het_test_before_outlier_extraction.RDS")

#Pleiotropy! Let's remove the outliers

mr_post <- remove_outlier(mr_steiger_df, rucker) 
#mr_post <- remove_outlier(mr_post, rucker) #run twice with updated rucker ofc, but no outliers are found...

saveRDS(mr_post, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/menorrhagia/instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure                    method nsnp           b         se         pval
# 1    exposure    outcome outcome exposure                  MR Egger  170  0.31550372 0.08371305 2.267294e-04
# 2    exposure    outcome outcome exposure           Weighted median  170  0.22365875 0.06018399 2.022013e-04
# 3    exposure    outcome outcome exposure Inverse variance weighted  170  0.17431993 0.03569570 1.042197e-06
# 4    exposure    outcome outcome exposure               Simple mode  170 -0.04465857 0.18377914 8.082991e-01
# 5    exposure    outcome outcome exposure             Weighted mode  170  0.33719486 0.09314363 3.888169e-04

#Let's check the plot:
#  
rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method     Estimate          SE       CI_low       CI_upp         P
# 1  Egger fixed effects -0.003641323 0.001882245 -0.007330456 4.781019e-05 0.0530439
# 2 Egger random effects -0.003641323 0.001882245 -0.007330456 4.781019e-05 0.9734781
# 
# [[1]]$Q
# Method          Q  df          P
# 1   Q_ivw 159.531816 169 0.68740277
# 2 Q_egger 156.055371 168 0.73602886
# 3  Q_diff   3.476445   1 0.06224833
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/whr/menorrhagia/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/menorrhagia/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/menorrhagia/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/menorrhagia/strict_het_test_after_outlier_extraction.RDS")

#######################
#Let's try anovulation#
#######################

#STEP 1: let's use only those that have F>10

mr_df <- liu_anovulation[which(((liu_anovulation$beta.exposure^2)/(liu_anovulation$se.exposure^2)) > 10),] #276

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

dir.create("output/2_replicating_liu_et_al//2_mr/2_mr/whr/anovulation")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/anovulation/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp         b         se         pval
# 1    exposure    outcome outcome exposure                  MR Egger  194 0.8425925 0.23329583 0.0003883017
# 2    exposure    outcome outcome exposure           Weighted median  194 0.3784274 0.14749777 0.0102983176
# 3    exposure    outcome outcome exposure Inverse variance weighted  194 0.2889640 0.09949459 0.0036805278
# 4    exposure    outcome outcome exposure               Simple mode  194 0.1548898 0.43941312 0.7248536961
# 5    exposure    outcome outcome exposure             Weighted mode  194 0.5361899 0.21515585 0.0135414804

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method    Estimate          SE      CI_low       CI_upp           P
# 1  Egger fixed effects -0.01432261 0.004713596 -0.02356108 -0.005084129 0.002377011
# 2 Egger random effects -0.01432261 0.005476866 -0.02505707 -0.003588147 0.995540065
# 
# [[1]]$Q
# Method          Q  df            P
# 1   Q_ivw 268.448306 193 0.0002672082
# 2 Q_egger 259.215377 192 0.0008772518
# 3  Q_diff   9.232929   1 0.0023770112
# 
# [[1]]$res
# [1] "D"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_steiger_df) #they disagree, but still

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/whr/anovulation/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/anovulation/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/anovulation/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/anovulation/strict_het_test_before_outlier_extraction.RDS")

#BMI is not linked to anovulation

mr_post <- remove_outlier(mr_steiger_df, rucker) 
#mr_post <- remove_outlier(mr_post, rucker) #run twice with updated rucker ofc, but no outliers are found...
 
saveRDS(mr_post, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/anovulation/instruments_after_outlier_removal_df.RDS")
 
#STEP 6: let's run the results for real now:
 
mr_res <- TwoSampleMR::mr(mr_post)
 
# id.exposure id.outcome outcome exposure                    method nsnp         b        se         pval
# 1    exposure    outcome outcome exposure                  MR Egger  175 0.7352752 0.2089461 0.0005536403
# 2    exposure    outcome outcome exposure           Weighted median  175 0.3300277 0.1445203 0.0223945873
# 3    exposure    outcome outcome exposure Inverse variance weighted  175 0.2029107 0.0891567 0.0228526749
# 4    exposure    outcome outcome exposure               Simple mode  175 0.1703723 0.3738374 0.6491459145
# 5    exposure    outcome outcome exposure             Weighted mode  175 0.5352923 0.2226570 0.0172617652

#Let's check the plot:
  
rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method   Estimate          SE      CI_low       CI_upp           P
# 1  Egger fixed effects -0.0138082 0.004424322 -0.02247971 -0.005136686 0.001802534
# 2 Egger random effects -0.0138082 0.004424322 -0.02247971 -0.005136686 0.999098733
# 
# [[1]]$Q
# Method          Q  df           P
# 1   Q_ivw 148.897696 174 0.916237011
# 2 Q_egger 140.961110 173 0.964503559
# 3  Q_diff   7.936586   1 0.004844514
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/whr/anovulation/plots_after_outlier_extraction.svg", height=16, width=16, units="in")
 
saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/anovulation/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/anovulation/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/anovulation/strict_het_test_after_outlier_extraction.RDS")

#####################
#Let's try hirsutism#
#####################

#STEP 1: let's use only those that have F>10

mr_df <- liu_hirsutism[which(((liu_hirsutism$beta.exposure^2)/(liu_hirsutism$se.exposure^2)) > 10),] #276

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

dir.create("output/2_replicating_liu_et_al//2_mr/2_mr/whr/hirsutism")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/hirsutism/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp         b        se         pval
# 1    exposure    outcome outcome exposure                  MR Egger  194 0.4680765 0.4727086 0.3233220861
# 2    exposure    outcome outcome exposure           Weighted median  194 0.7484299 0.3758769 0.0464636497
# 3    exposure    outcome outcome exposure Inverse variance weighted  194 0.7430471 0.1979824 0.0001746631
# 4    exposure    outcome outcome exposure               Simple mode  194 0.3326228 0.7435195 0.6551149624
# 5    exposure    outcome outcome exposure             Weighted mode  194 0.6777709 0.5224158 0.1960496939

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method   Estimate         SE      CI_low     CI_upp         P
# 1  Egger fixed effects 0.00710801 0.01081485 -0.01408871 0.02830473 0.5110232
# 2 Egger random effects 0.00710801 0.01109255 -0.01463300 0.02884902 0.2608292
# 
# [[1]]$Q
# Method           Q  df         P
# 1   Q_ivw 202.4188197 193 0.3065595
# 2 Q_egger 201.9868484 192 0.2963436
# 3  Q_diff   0.4319713   1 0.5110232
# 
# [[1]]$res
# [1] "A"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_steiger_df) #they disagree, but still

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/whr/hirsutism/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/hirsutism/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/hirsutism/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/hirsutism/strict_het_test_before_outlier_extraction.RDS")

#Quite good actually, we do not need anything else

# mr_post <- remove_outlier(mr_steiger_df, rucker) 
# #mr_post <- remove_outlier(mr_post, rucker) #run twice with updated rucker ofc, but no outliers are found...
# 
# saveRDS(mr_post, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/hirsutism/instruments_after_outlier_removal_df.RDS")
# 
# #STEP 6: let's run the results for real now:
# 
# mr_res <- TwoSampleMR::mr(mr_post)
# 
# # id.exposure id.outcome outcome exposure                    method nsnp         b         se         pval
# # 1    exposure    outcome outcome exposure                  MR Egger  244 0.0195787 0.09640408 8.392354e-01
# # 2    exposure    outcome outcome exposure           Weighted median  244 0.1780444 0.05638549 1.590584e-03
# # 3    exposure    outcome outcome exposure Inverse variance weighted  244 0.2231649 0.03436266 8.336465e-11
# # 4    exposure    outcome outcome exposure               Simple mode  244 0.3047140 0.16828207 7.141711e-02
# # 5    exposure    outcome outcome exposure             Weighted mode  244 0.1927115 0.10166373 5.920252e-02
# 
# #Let's check the plot:
# #  
# rucker <- TwoSampleMR::mr_rucker(mr_post)
# 
# # [[1]]$intercept
# # Method    Estimate          SE       CI_low     CI_upp           P
# # 1  Egger fixed effects 0.004384173 0.001844563 0.0007688961 0.00799945 0.017463149
# # 2 Egger random effects 0.004384173 0.001844563 0.0007688961 0.00799945 0.008731575
# # 
# # [[1]]$Q
# # Method          Q  df          P
# # 1   Q_ivw 223.957769 243 0.80427150
# # 2 Q_egger 218.848984 242 0.85481831
# # 3  Q_diff   5.108784   1 0.02380499
# # 
# # [[1]]$res
# # [1] "A"
# 
# # #Here there is heterogeneity and pleitropy
# 
# strict_het <- compute_strict_het(mr_post) #they disagree, but still
# 
# #That is quite good:
# 
# mr_post$labels <- NA
# mr_post$units.exposure <- "SD"
# mr_post$units.outcome <- "LOR"
# 
# plotio=mr_plots(mr_post)
# ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/whr/hirsutism/plots_after_outlier_extraction.svg", height=16, width=16, units="in")
# 
# saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/hirsutism/mr_res_after_outlier_extraction.RDS")
# saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/hirsutism/sensitivity_test_after_outlier_extraction.RDS")
# saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/whr/hirsutism/strict_het_test_after_outlier_extraction.RDS")
# 
