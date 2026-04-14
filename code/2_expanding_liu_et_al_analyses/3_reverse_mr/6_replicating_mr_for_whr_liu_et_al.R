##############
#INTRODUCTION#
##############

#This code performs MR analyses between PCOS and BMI - replicating the results from Liu et al.

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

liu_pcos_doctor <- fread("output/2_replicating_liu_et_al/3_reverse_mr/1_expo_outcome_df/pcos_finngen_WHR.txt")
liu_pcos_broad <- fread("output/2_replicating_liu_et_al/3_reverse_mr/1_expo_outcome_df/pcos_broad_WHR.txt")
liu_pcos_consortium <- fread("output/2_replicating_liu_et_al/3_reverse_mr/1_expo_outcome_df/pcos_consortium_WHR.txt")
liu_pcos_meta_analysis <- fread("output/2_replicating_liu_et_al/3_reverse_mr/1_expo_outcome_df/pcos_day_WHR.txt")
liu_pcos_venkatesh <- fread("output/2_replicating_liu_et_al/3_reverse_mr/1_expo_outcome_df/pcos_venkatesh_WHR.txt")
liu_pcos_adj_age <- fread("output/2_replicating_liu_et_al/3_reverse_mr/1_expo_outcome_df/pcos_tyrmi_WHR.txt")
liu_pcos_adj_age_bmi <- fread("output/2_replicating_liu_et_al/3_reverse_mr/1_expo_outcome_df/pcosadjbmi_tyrmi_WHR.txt")

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

mr_df$r.outcome <- TwoSampleMR::get_r_from_pn(
  p = as.numeric(mr_df$pval.outcome),n = as.numeric(mr_df$samplesize.outcome))

mr_df$r.exposure <- TwoSampleMR::get_r_from_lor(
  lor = as.numeric(mr_df$beta.exposure),
  af = as.numeric(mr_df$eaf.exposure),
  ncase = as.numeric(mr_df$ncase.exposure),
  ncontrol = as.numeric(mr_df$ncontrol.exposure),
  prevalence = mr_df$prevalence.exposure[1] 
)

#Adding additional columns to make steiger run...

mr_df$id.outcome <- "WHR"
mr_df$id.exposure <- "PCOS (meta-analysis)"
mr_df$outcome <- "WHR"
mr_df$exposure <- "PCOS (meta-analysis)"
mr_df$units.outcome <- "SD"
mr_df$units.exposure <- "LOR"

mr_steiger_df <- TwoSampleMR::steiger_filtering(mr_df)
mr_steiger_df <- mr_steiger_df[which(mr_steiger_df$steiger_dir == TRUE),] 
mr_steiger_df$mr_keep <- TRUE

#Let's save the data:

dir.create("output/2_replicating_liu_et_al/3_reverse_mr/2_mr/")
dir.create("output/2_replicating_liu_et_al/3_reverse_mr/2_mr/whr/")
dir.create("output/2_replicating_liu_et_al/3_reverse_mr/2_mr/whr/meta_analysis/")
saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al/3_reverse_mr/2_mr/whr/meta_analysis/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome             exposure                    method nsnp           b         se      pval
# 1 PCOS (meta-analysis)        WHR     WHR PCOS (meta-analysis)                  MR Egger   14 0.061491982 0.06359853 0.3526884
# 2 PCOS (meta-analysis)        WHR     WHR PCOS (meta-analysis)           Weighted median   14 0.007204944 0.01111335 0.5167816
# 3 PCOS (meta-analysis)        WHR     WHR PCOS (meta-analysis) Inverse variance weighted   14 0.009450339 0.01357175 0.4862259
# 4 PCOS (meta-analysis)        WHR     WHR PCOS (meta-analysis)               Simple mode   14 0.014986262 0.02493037 0.5580934
# 5 PCOS (meta-analysis)        WHR     WHR PCOS (meta-analysis)             Weighted mode   14 0.008495380 0.02272898 0.7146013

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method     Estimate          SE      CI_low       CI_upp          P
# 1  Egger fixed effects -0.006708668 0.003698838 -0.01395826 0.0005409218 0.06972038
# 2 Egger random effects -0.006708668 0.008005136 -0.02239845 0.0089811095 0.79899744
# 
# [[1]]$Q
# Method         Q df            P
# 1   Q_ivw 59.496270 13 6.463192e-08
# 2 Q_egger 56.206678 12 1.098705e-07
# 3  Q_diff  3.289591  1 6.972038e-02
# 
# [[1]]$res
# [1] "B"

strict_het <- compute_strict_het(mr_steiger_df)

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "LOR"
mr_steiger_df$units.outcome <- "SD"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al/3_reverse_mr/2_mr/whr/meta_analysis/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al/3_reverse_mr/2_mr/whr/meta_analysis/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al/3_reverse_mr/2_mr/whr/meta_analysis//sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al/3_reverse_mr/2_mr/whr/meta_analysis//strict_het_test_before_outlier_extraction.RDS")

#A lot of heterogeneity

mr_post <- remove_outlier(mr_steiger_df, rucker)

saveRDS(mr_post, "output/2_replicating_liu_et_al/3_reverse_mr/2_mr/whr/meta_analysis/instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome             exposure                    method nsnp           b         se      pval
# 1 PCOS (meta-analysis)        WHR     WHR PCOS (meta-analysis)                  MR Egger    8 0.087024385 0.05069733 0.1368757
# 2 PCOS (meta-analysis)        WHR     WHR PCOS (meta-analysis)           Weighted median    8 0.008582837 0.01131901 0.4482910
# 3 PCOS (meta-analysis)        WHR     WHR PCOS (meta-analysis) Inverse variance weighted    8 0.013806676 0.00911914 0.1300175
# 4 PCOS (meta-analysis)        WHR     WHR PCOS (meta-analysis)               Simple mode    8 0.008158114 0.01933879 0.6857862
# 5 PCOS (meta-analysis)        WHR     WHR PCOS (meta-analysis)             Weighted mode    8 0.006978473 0.01786181 0.7076424

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method     Estimate          SE      CI_low      CI_upp         P
# 1  Egger fixed effects -0.008876212 0.005919802 -0.02047881 0.002726388 0.1337673
# 2 Egger random effects -0.008876212 0.005919802 -0.02047881 0.002726388 0.9331164
# 
# [[1]]$Q
# Method        Q df         P
# 1   Q_ivw 7.878608  7 0.3434199
# 2 Q_egger 5.731126  6 0.4539734
# 3  Q_diff 2.147482  1 0.1428039
# 
# [[1]]$res
# [1] "A"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_post) #heterogeneity removed

#That is quite good:

mr_post$labels <- NA
mr_post$units.exposure <- "LOR"
mr_post$units.outcome <- "SD"

plotio=mr_plots(mr_post)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/meta_analysis/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/meta_analysis/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/meta_analysis/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/meta_analysis/strict_het_test_after_outlier_extraction.RDS")

###############################################################
#Let's first analyse the data for the whole set for bmi - PCOS#
###############################################################

#STEP 1: let's use only those that have F>10

mr_df <- liu_pcos_doctor[which(((liu_pcos_doctor$beta.exposure^2)/(liu_pcos_doctor$se.exposure^2)) > 10),] 

#STEP 2: let's compute MR steiger:

mr_df$beta.exposure <- as.numeric(mr_df$beta.exposure)
mr_df$se.exposure <- as.numeric(mr_df$se.exposure)
mr_df$eaf.exposure <- as.numeric(mr_df$eaf.exposure)

mr_df$beta.outcome <- as.numeric(mr_df$beta.outcome)
mr_df$se.outcome <- as.numeric(mr_df$se.outcome)
mr_df$eaf.outcome <- as.numeric(mr_df$eaf.outcome)

#STEP 2: let's compute MR steiger:

mr_df$r.outcome <- TwoSampleMR::get_r_from_pn(
  p = as.numeric(mr_df$pval.outcome),n = as.numeric(mr_df$samplesize.outcome))

mr_df$r.exposure <- TwoSampleMR::get_r_from_lor(
  lor = as.numeric(mr_df$beta.exposure),
  af = as.numeric(mr_df$eaf.exposure),
  ncase = as.numeric(mr_df$ncase.exposure),
  ncontrol = as.numeric(mr_df$ncontrol.exposure),
  prevalence = mr_df$prevalence.exposure[1] 
)

mr_steiger_df <- TwoSampleMR::steiger_filtering(mr_df)
mr_steiger_df <- mr_steiger_df[which(mr_steiger_df$steiger_dir == TRUE),]
mr_steiger_df$mr_keep <- TRUE

#Let's save the data:

dir.create("output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/")
dir.create("output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/doctor")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/doctor/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp           b          se      pval
# 1    exposure    outcome outcome exposure                  MR Egger    6 0.041432367 0.026316888 0.1905193
# 2    exposure    outcome outcome exposure           Weighted median    6 0.006380350 0.008592208 0.4577397
# 3    exposure    outcome outcome exposure Inverse variance weighted    6 0.001672174 0.012305279 0.8919076
# 4    exposure    outcome outcome exposure               Simple mode    6 0.002769548 0.015979809 0.8692014
# 5    exposure    outcome outcome exposure             Weighted mode    6 0.008520029 0.015220481 0.5997821

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method    Estimate          SE      CI_low       CI_upp            P
# 1  Egger fixed effects -0.01120697 0.003116601 -0.01731540 -0.005098549 0.0003232761
# 2 Egger random effects -0.01120697 0.006788240 -0.02451168  0.002097732 0.9506245434
# 
# [[1]]$Q
# Method        Q df            P
# 1   Q_ivw 31.90680  5 6.199155e-06
# 2 Q_egger 18.97632  4 7.944075e-04
# 3  Q_diff 12.93047  1 3.232761e-04
# 
# [[1]]$res
# [1] "D"

strict_het <- compute_strict_het(mr_steiger_df)

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "LOR"
mr_steiger_df$units.outcome <- "SD"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/doctor/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/doctor/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/doctor/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/doctor/strict_het_test_before_outlier_extraction.RDS")

#Biased by pleiotropy  (which is in part due to small number of instruments)

mr_post <- remove_outlier(mr_steiger_df, rucker) #no significant outliers...

saveRDS(mr_post, "output/2_replicating_liu_et_al/3_reverse_mr/2_mr/whr/doctor/instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure                    method nsnp           b          se      pval
# 1    exposure    outcome outcome exposure                  MR Egger    6 0.041432367 0.026316888 0.1905193
# 2    exposure    outcome outcome exposure           Weighted median    6 0.006380350 0.008424903 0.4488580
# 3    exposure    outcome outcome exposure Inverse variance weighted    6 0.001672174 0.012305279 0.8919076
# 4    exposure    outcome outcome exposure               Simple mode    6 0.002769548 0.015514136 0.8653243
# 5    exposure    outcome outcome exposure             Weighted mode    6 0.008520029 0.014159130 0.5735923

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method    Estimate          SE      CI_low       CI_upp            P
# 1  Egger fixed effects -0.01120697 0.003116601 -0.01731540 -0.005098549 0.0003232761
# 2 Egger random effects -0.01120697 0.006788240 -0.02451168  0.002097732 0.9506245434
# 
# [[1]]$Q
# Method        Q df            P
# 1   Q_ivw 31.90680  5 6.199155e-06
# 2 Q_egger 18.97632  4 7.944075e-04
# 3  Q_diff 12.93047  1 3.232761e-04
# 
# [[1]]$res
# [1] "D"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_post) #heterogeneity removed

#That is quite good:

mr_post$labels <- NA
mr_post$units.exposure <- "LOR"
mr_post$units.outcome <- "SD"

plotio=mr_plots(mr_post)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/doctor/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/doctor/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/doctor/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/doctor/strict_het_test_after_outlier_extraction.RDS")

###############################################################
#Let's first analyse the data for the whole set for bmi - PCOS# broad definition - anovulation
###############################################################

#STEP 1: let's use only those that have F>10

mr_df <- liu_pcos_broad[which(((liu_pcos_broad$beta.exposure^2)/(liu_pcos_broad$se.exposure^2)) > 10),] 

#STEP 2: let's compute MR steiger:

mr_df$beta.exposure <- as.numeric(mr_df$beta.exposure)
mr_df$se.exposure <- as.numeric(mr_df$se.exposure)
mr_df$eaf.exposure <- as.numeric(mr_df$eaf.exposure)

mr_df$beta.outcome <- as.numeric(mr_df$beta.outcome)
mr_df$se.outcome <- as.numeric(mr_df$se.outcome)
mr_df$eaf.outcome <- as.numeric(mr_df$eaf.outcome)

#STEP 2: let's compute MR steiger:

mr_df$r.outcome <- TwoSampleMR::get_r_from_pn(
  p = as.numeric(mr_df$pval.outcome),n = as.numeric(mr_df$samplesize.outcome))

mr_df$r.exposure <- TwoSampleMR::get_r_from_lor(
  lor = as.numeric(mr_df$beta.exposure),
  af = as.numeric(mr_df$eaf.exposure),
  ncase = as.numeric(mr_df$ncase.exposure),
  ncontrol = as.numeric(mr_df$ncontrol.exposure),
  prevalence = mr_df$prevalence.exposure[1] 
)


mr_steiger_df <- TwoSampleMR::steiger_filtering(mr_df)
mr_steiger_df <- mr_steiger_df[which(mr_steiger_df$steiger_dir == TRUE),]
mr_steiger_df$mr_keep <- TRUE

#Let's save the data:

dir.create("output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/broad")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/broad/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp           b         se       pval
# 1    exposure    outcome outcome exposure                  MR Egger   12 0.050205355 0.03054271 0.13125058
# 2    exposure    outcome outcome exposure           Weighted median   12 0.012856238 0.01015509 0.20551696
# 3    exposure    outcome outcome exposure Inverse variance weighted   12 0.024699078 0.01133877 0.02938484
# 4    exposure    outcome outcome exposure               Simple mode   12 0.011508198 0.01394873 0.42688131
# 5    exposure    outcome outcome exposure             Weighted mode   12 0.009832489 0.01365834 0.48660800

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method     Estimate          SE       CI_low      CI_upp          P
# 1  Egger fixed effects -0.004673723 0.002662037 -0.009891221 0.000543774 0.07914065
# 2 Egger random effects -0.004673723 0.005189373 -0.014844707 0.005497260 0.81610841
# 
# [[1]]$Q
# Method         Q df            P
# 1   Q_ivw 41.084023 11 2.331108e-05
# 2 Q_egger 38.001561 10 3.792796e-05
# 3  Q_diff  3.082462  1 7.914065e-02
# 
# [[1]]$res
# [1] "B"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_steiger_df) #they disagree, but still

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "LOR"
mr_steiger_df$units.outcome <- "SD"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/broad/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/broad/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/broad/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/broad/strict_het_test_before_outlier_extraction.RDS")

#One variant is an outlier removing the association!

mr_post <- remove_outlier(mr_steiger_df, rucker)

saveRDS(mr_post, "output/2_replicating_liu_et_al/3_reverse_mr/2_mr/whr/broad/instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure                    method nsnp          b          se       pval
# 1    exposure    outcome outcome exposure                  MR Egger    8 0.01197918 0.024778891 0.64591712
# 2    exposure    outcome outcome exposure           Weighted median    8 0.01101984 0.009698502 0.25585534
# 3    exposure    outcome outcome exposure Inverse variance weighted    8 0.01490963 0.007903227 0.05922435
# 4    exposure    outcome outcome exposure               Simple mode    8 0.01127119 0.012709860 0.40461922
# 5    exposure    outcome outcome exposure             Weighted mode    8 0.00948932 0.011946528 0.45309426

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method     Estimate          SE       CI_low      CI_upp         P
# 1  Egger fixed effects 0.0005042404 0.003502918 -0.006361353 0.007369834 0.8855410
# 2 Egger random effects 0.0005042404 0.004003389 -0.007342258 0.008350739 0.4498844
# 
# [[1]]$Q
# Method          Q df         P
# 1   Q_ivw 7.85766756  7 0.3453119
# 2 Q_egger 7.83694634  6 0.2502941
# 3  Q_diff 0.02072122  1 0.8855410
# 
# [[1]]$res
# [1] "A"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_post) #heterogeneity removed

#That is quite good:

mr_post$labels <- NA
mr_post$units.exposure <- "LOR"
mr_post$units.outcome <- "SD"

plotio=mr_plots(mr_post)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/broad/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/broad/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/broad/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/broad/strict_het_test_after_outlier_extraction.RDS")

###############################################################
#Let's first analyse the data for the whole set for bmi - PCOS# consortium definition - anovulation
###############################################################

#STEP 1: let's use only those that have F>10

mr_df <- liu_pcos_consortium[which(((liu_pcos_consortium$beta.exposure^2)/(liu_pcos_consortium$se.exposure^2)) > 10),] 

#STEP 2: let's compute MR steiger:

mr_df$beta.exposure <- as.numeric(mr_df$beta.exposure)
mr_df$se.exposure <- as.numeric(mr_df$se.exposure)
mr_df$eaf.exposure <- as.numeric(mr_df$eaf.exposure)

mr_df$beta.outcome <- as.numeric(mr_df$beta.outcome)
mr_df$se.outcome <- as.numeric(mr_df$se.outcome)
mr_df$eaf.outcome <- as.numeric(mr_df$eaf.outcome)

#STEP 2: let's compute MR steiger:

mr_df$r.outcome <- TwoSampleMR::get_r_from_pn(
  p = as.numeric(mr_df$pval.outcome),n = as.numeric(mr_df$samplesize.outcome))

mr_df$r.exposure <- TwoSampleMR::get_r_from_lor(
  lor = as.numeric(mr_df$beta.exposure),
  af = as.numeric(mr_df$eaf.exposure),
  ncase = as.numeric(mr_df$ncase.exposure),
  ncontrol = as.numeric(mr_df$ncontrol.exposure),
  prevalence = mr_df$prevalence.exposure[1] 
)

mr_steiger_df <- TwoSampleMR::steiger_filtering(mr_df)
mr_steiger_df <- mr_steiger_df[which(mr_steiger_df$steiger_dir == TRUE),]
mr_steiger_df$mr_keep <- TRUE

#Let's save the data:

dir.create("output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/consortium")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/consortium/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp          b         se        pval
# 1    exposure    outcome outcome exposure                  MR Egger    5 0.14008674 0.05430996 0.081822883
# 2    exposure    outcome outcome exposure           Weighted median    5 0.09596419 0.03084172 0.001861358
# 3    exposure    outcome outcome exposure Inverse variance weighted    5 0.07528791 0.02330738 0.001236961
# 4    exposure    outcome outcome exposure               Simple mode    5 0.10200330 0.04114911 0.068293784
# 5    exposure    outcome outcome exposure             Weighted mode    5 0.10200330 0.03648070 0.049006570

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method     Estimate          SE      CI_low       CI_upp         P
# 1  Egger fixed effects -0.004565237 0.002813365 -0.01007933 0.0009488575 0.1046544
# 2 Egger random effects -0.004565237 0.002813365 -0.01007933 0.0009488575 0.9476728
# 
# [[1]]$Q
# Method        Q df         P
# 1   Q_ivw 3.732967  4 0.4433493
# 2 Q_egger 1.988038  3 0.5748932
# 3  Q_diff 1.744929  1 0.1865155
# 
# [[1]]$res
# [1] "A"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_steiger_df) #they disagree, but still

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "LOR"
mr_steiger_df$units.outcome <- "SD"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/consortium/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/consortium/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/consortium/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/consortium/strict_het_test_before_outlier_extraction.RDS")

#Asymmetry due to low variants. Let's see...

#STEP 5: removing outliers!

mr_post <- remove_outlier(mr_steiger_df, rucker) #no significant outliers
  
# saveRDS(mr_post, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/consortium/instruments_after_outlier_removal_df.RDS")
#   
# #STEP 6: let's run the results for real now:
#   
# mr_res <- TwoSampleMR::mr(mr_post)
#   
# # id.exposure id.outcome outcome exposure                    method nsnp          b         se         pval
# # 1    exposure    outcome outcome exposure                  MR Egger    5 0.14008674 0.05430996 0.0818228833
# # 2    exposure    outcome outcome exposure           Weighted median    5 0.09596419 0.02858347 0.0007869863
# # 3    exposure    outcome outcome exposure Inverse variance weighted    5 0.07528791 0.02330738 0.0012369613
# # 4    exposure    outcome outcome exposure               Simple mode    5 0.10200330 0.04112303 0.0681788019
# # 5    exposure    outcome outcome exposure             Weighted mode    5 0.10200330 0.03731347 0.0522441440
# 
# # #Let's check the plot:
# #  
# rucker <- TwoSampleMR::mr_rucker(mr_post)
# 
# # [[1]]$intercept
# # Method     Estimate          SE      CI_low       CI_upp         P
# # 1  Egger fixed effects -0.004565237 0.002813365 -0.01007933 0.0009488575 0.1046544
# # 2 Egger random effects -0.004565237 0.002813365 -0.01007933 0.0009488575 0.9476728
# # 
# # [[1]]$Q
# # Method        Q df         P
# # 1   Q_ivw 3.732967  4 0.4433493
# # 2 Q_egger 1.988038  3 0.5748932
# # 3  Q_diff 1.744929  1 0.1865155
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
# mr_post$units.exposure <- "LOR"
# mr_post$units.outcome <- "SD"
#   
# plotio=mr_plots(mr_post)
# ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/consortium/plots_after_outlier_extraction.svg", height=16, width=16, units="in")
# 
# saveRDS(mr_res, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/consortium/mr_res_after_outlier_extraction.RDS")
# saveRDS(rucker, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/consortium/sensitivity_test_after_outlier_extraction.RDS")
# saveRDS(strict_het, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/consortium/strict_het_test_after_outlier_extraction.RDS")

###############################################################
#Let's first analyse the data for the whole set for bmi - PCOS# age adjusted
###############################################################

#STEP 1: let's use only those that have F>10

mr_df <- liu_pcos_adj_age[which(((liu_pcos_adj_age$beta.exposure^2)/(liu_pcos_adj_age$se.exposure^2)) > 10),] 

#Here we are just gonna add the eaf from outcome to exposure.

mr_df$eaf.exposure=mr_df$eaf.outcome #a trick to run the analyses. The frequencies match. I checked them all manually

#STEP 2: let's compute MR steiger:

mr_df$beta.exposure <- as.numeric(mr_df$beta.exposure)
mr_df$se.exposure <- as.numeric(mr_df$se.exposure)
mr_df$eaf.exposure <- as.numeric(mr_df$eaf.exposure)

mr_df$beta.outcome <- as.numeric(mr_df$beta.outcome)
mr_df$se.outcome <- as.numeric(mr_df$se.outcome)
mr_df$eaf.outcome <- as.numeric(mr_df$eaf.outcome)

#STEP 2: let's compute MR steiger:

mr_df$r.outcome <- TwoSampleMR::get_r_from_pn(
  p = as.numeric(mr_df$pval.outcome),n = as.numeric(mr_df$samplesize.outcome))

mr_df$r.exposure <- TwoSampleMR::get_r_from_lor(
  lor = as.numeric(mr_df$beta.exposure),
  af = as.numeric(mr_df$eaf.exposure),
  ncase = as.numeric(mr_df$ncase.exposure),
  ncontrol = as.numeric(mr_df$ncontrol.exposure),
  prevalence = mr_df$prevalence.exposure[1] 
)

mr_steiger_df <- TwoSampleMR::steiger_filtering(mr_df)
mr_steiger_df <- mr_steiger_df[which(mr_steiger_df$steiger_dir == TRUE),]
mr_steiger_df$mr_keep <- TRUE

#Let's save the data:

dir.create("output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/age_adjusted")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/age_adjusted/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp            b         se      pval
# 1    exposure    outcome outcome exposure                  MR Egger    6  0.100817996 0.06285990 0.1840123
# 2    exposure    outcome outcome exposure           Weighted median    6 -0.001151857 0.01350149 0.9320123
# 3    exposure    outcome outcome exposure Inverse variance weighted    6  0.003902897 0.01890335 0.8364270
# 4    exposure    outcome outcome exposure               Simple mode    6 -0.012296286 0.01966573 0.5592231
# 5    exposure    outcome outcome exposure             Weighted mode    6 -0.016682674 0.01436970 0.2980674

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method    Estimate          SE      CI_low       CI_upp            P
# 1  Egger fixed effects -0.02011457 0.005346363 -0.03059325 -0.009635888 0.0001683648
# 2 Egger random effects -0.02011457 0.012588323 -0.04478723  0.004558094 0.9449645947
# 
# [[1]]$Q
# Method        Q df            P
# 1   Q_ivw 36.33060  5 8.156110e-07
# 2 Q_egger 22.17578  4 1.849018e-04
# 3  Q_diff 14.15482  1 1.683648e-04
# 
# [[1]]$res
# [1] "D"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_steiger_df) #they disagree, but still

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "LOR"
mr_steiger_df$units.outcome <- "SD"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/age_adjusted/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/age_adjusted/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/age_adjusted/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/age_adjusted/strict_het_test_before_outlier_extraction.RDS")

#There is a bit of asymmetry. I am almost certain that Wmed and Wmod correct for this, but we are gonna run the analyses again just to be on the safe side.

mr_post <- remove_outlier(mr_steiger_df, rucker)

saveRDS(mr_post, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/age_adjusted/instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure                    method nsnp             b         se      pval
# 1    exposure    outcome outcome exposure                  MR Egger    4  0.0940814806 0.11128457 0.4868949
# 2    exposure    outcome outcome exposure           Weighted median    4 -0.0011747006 0.01401070 0.9331812
# 3    exposure    outcome outcome exposure Inverse variance weighted    4 -0.0019318895 0.01993723 0.9228069
# 4    exposure    outcome outcome exposure               Simple mode    4 -0.0008279725 0.01693139 0.9640712
# 5    exposure    outcome outcome exposure             Weighted mode    4 -0.0056785341 0.01489340 0.7283971

# #Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method   Estimate          SE      CI_low        CI_upp         P
# 1  Egger fixed effects -0.0173382 0.008607785 -0.03420915 -0.0004672493 0.0439837
# 2 Egger random effects -0.0173382 0.019743666 -0.05603507  0.0213586757 0.8100729
# 
# [[1]]$Q
# Method         Q df           P
# 1   Q_ivw 14.579295  3 0.002213863
# 2 Q_egger 10.522107  2 0.005189833
# 3  Q_diff  4.057188  1 0.043983695
# 
# [[1]]$res
# [1] "D"

# #Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_post) #they disagree, but still

#That is quite good:

mr_post$labels <- NA
mr_post$units.exposure <- "LOR"
mr_post$units.outcome <- "SD"

plotio=mr_plots(mr_post)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/age_adjusted/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/age_adjusted/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/age_adjusted/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/age_adjusted/strict_het_test_after_outlier_extraction.RDS")

###############################################################
#Let's first analyse the data for the whole set for bmi - PCOS# age and BMI adjusted
###############################################################

#STEP 1: let's use only those that have F>10

mr_df <- liu_pcos_adj_age_bmi[which(((liu_pcos_adj_age_bmi$beta.exposure^2)/(liu_pcos_adj_age_bmi$se.exposure^2)) > 10),] #250
mr_df$eaf.exposure=mr_df$eaf.outcome #all correctly aligned. They make sense

#STEP 2: let's compute MR steiger:

mr_df$beta.exposure <- as.numeric(mr_df$beta.exposure)
mr_df$se.exposure <- as.numeric(mr_df$se.exposure)
mr_df$eaf.exposure <- as.numeric(mr_df$eaf.exposure)

mr_df$beta.outcome <- as.numeric(mr_df$beta.outcome)
mr_df$se.outcome <- as.numeric(mr_df$se.outcome)
mr_df$eaf.outcome <- as.numeric(mr_df$eaf.outcome)

#STEP 2: let's compute MR steiger:

mr_df$r.outcome <- TwoSampleMR::get_r_from_pn(
  p = as.numeric(mr_df$pval.outcome),n = as.numeric(mr_df$samplesize.outcome))

mr_df$r.exposure <- TwoSampleMR::get_r_from_lor(
  lor = as.numeric(mr_df$beta.exposure),
  af = as.numeric(mr_df$eaf.exposure),
  ncase = as.numeric(mr_df$ncase.exposure),
  ncontrol = as.numeric(mr_df$ncontrol.exposure),
  prevalence = mr_df$prevalence.exposure[1] 
)


mr_steiger_df <- TwoSampleMR::steiger_filtering(mr_df)
mr_steiger_df <- mr_steiger_df[which(mr_steiger_df$steiger_dir == TRUE),]
mr_steiger_df$mr_keep <- TRUE

#Let's save the data:

dir.create("output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/age_and_bmi_adjusted")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/age_and_bmi_adjusted/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp           b         se       pval
# 1    exposure    outcome outcome exposure                  MR Egger    3  0.17830835 0.03228667 0.11403848
# 2    exposure    outcome outcome exposure           Weighted median    3  0.01842493 0.01443683 0.20186914
# 3    exposure    outcome outcome exposure Inverse variance weighted    3  0.01915482 0.03462596 0.58013162
# 4    exposure    outcome outcome exposure               Simple mode    3  0.07025584 0.01984255 0.07133775
# 5    exposure    outcome outcome exposure             Weighted mode    3 -0.03352081 0.02052935 0.24410459

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method    Estimate          SE     CI_low      CI_upp            P
# 1  Egger fixed effects -0.03788505 0.002176956 -0.0421518 -0.03361829 7.862028e-68
# 2 Egger random effects -0.03788505 0.002176956 -0.0421518 -0.03361829 1.000000e+00
# 
# [[1]]$Q
# Method           Q df            P
# 1   Q_ivw 26.67941216  2 1.609308e-06
# 2 Q_egger  0.08780279  1 7.669893e-01
# 3  Q_diff 26.59160937  1 2.513329e-07
# 
# [[1]]$res
# [1] "C"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_steiger_df) #they disagree, but still

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "LOR"
mr_steiger_df$units.outcome <- "SD"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/age_and_bmi_adjusted/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/age_and_bmi_adjusted/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/age_and_bmi_adjusted/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/age_and_bmi_adjusted/strict_het_test_before_outlier_extraction.RDS")

#Let's remove outliers:

mr_post <- remove_outlier(mr_steiger_df, rucker)

saveRDS(mr_post, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/age_and_bmi_adjusted/instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure     method nsnp          b         se        pval
# 1    exposure    outcome outcome exposure Wald ratio    1 0.08060942 0.02825485 0.004331664

# #Let's check the plot:
#  
#rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method   Estimate          SE       CI_low     CI_upp         P
# 1  Egger fixed effects 0.00571065 0.006517913 -0.007064224 0.01848552 0.3809501
# 2 Egger random effects 0.00571065 0.006517913 -0.007064224 0.01848552 0.1904750
# 
# [[1]]$Q
# Method         Q df         P
# 1   Q_ivw 1.6801867  3 0.6413477
# 2 Q_egger 1.2141684  2 0.5449375
# 3  Q_diff 0.4660183  1 0.4948247
# 
# [[1]]$res
# [1] "A"

# #Here there is heterogeneity and pleitropy

#strict_het <- compute_strict_het(mr_post) #they disagree, but still

#That is quite good:

# mr_post$labels <- NA
# mr_post$units.exposure <- "LOR"
# mr_post$units.outcome <- "SD"
# 
# plotio=mr_plots(mr_post)
# ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/age_and_bmi_adjusted/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/age_and_bmi_adjusted/mr_res_after_outlier_extraction.RDS")
# saveRDS(rucker, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/age_and_bmi_adjusted/sensitivity_test_after_outlier_extraction.RDS")
# saveRDS(strict_het, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/age_and_bmi_adjusted/strict_het_test_after_outlier_extraction.RDS")

##########################################################
#STEP 2: let's first rerun the analyses for the Venkatesh#
##########################################################

mr_df <- liu_pcos_venkatesh[which(((as.numeric(liu_pcos_venkatesh$beta.exposure)^2)/(as.numeric(liu_pcos_venkatesh$se.exposure)^2)) > 10),] #17/17

#Let's make all numeric or it will go crazy

mr_df$beta.exposure <- as.numeric(mr_df$beta.exposure)
mr_df$se.exposure <- as.numeric(mr_df$se.exposure)
mr_df$eaf.exposure <- as.numeric(mr_df$eaf.exposure)

mr_df$beta.outcome <- as.numeric(mr_df$beta.outcome)
mr_df$se.outcome <- as.numeric(mr_df$se.outcome)
mr_df$eaf.outcome <- as.numeric(mr_df$eaf.outcome)

#STEP 2: let's compute MR steiger:

mr_df$r.outcome <- TwoSampleMR::get_r_from_pn(
  p = as.numeric(mr_df$pval.outcome),n = as.numeric(mr_df$samplesize.outcome))

mr_df$r.exposure <- TwoSampleMR::get_r_from_lor(
  lor = as.numeric(mr_df$beta.exposure),
  af = as.numeric(mr_df$eaf.exposure),
  ncase = as.numeric(mr_df$ncase.exposure),
  ncontrol = as.numeric(mr_df$ncontrol.exposure),
  prevalence = mr_df$prevalence.exposure[1] 
)

#Adding additional columns to make steiger run...

mr_df$id.outcome <- "WHR"
mr_df$id.exposure <- "PCOS (Venkatesh)"
mr_df$outcome <- "WHR"
mr_df$exposure <- "PCOS (Venkatesh)"
mr_df$units.outcome <- "SD"
mr_df$units.exposure <- "LOR"

mr_steiger_df <- TwoSampleMR::steiger_filtering(mr_df)
mr_steiger_df <- mr_steiger_df[which(mr_steiger_df$steiger_dir == TRUE),] 
mr_steiger_df$mr_keep <- TRUE

#Let's save the data:

dir.create("output/2_replicating_liu_et_al/3_reverse_mr/2_mr/")
dir.create("output/2_replicating_liu_et_al/3_reverse_mr/2_mr/whr/")
dir.create("output/2_replicating_liu_et_al/3_reverse_mr/2_mr/whr/venkatesh/")
saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al/3_reverse_mr/2_mr/whr/venkatesh/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome         exposure                    method nsnp          b         se         pval
# 1 PCOS (Venkatesh)        BMI     BMI PCOS (Venkatesh)                  MR Egger   16 0.03952724 0.03795348 0.3153116340
# 2 PCOS (Venkatesh)        BMI     BMI PCOS (Venkatesh)           Weighted median   16 0.03978611 0.01204030 0.0009517624
# 3 PCOS (Venkatesh)        BMI     BMI PCOS (Venkatesh) Inverse variance weighted   16 0.02706482 0.01420938 0.0568171668
# 4 PCOS (Venkatesh)        BMI     BMI PCOS (Venkatesh)               Simple mode   16 0.05729115 0.01964903 0.0106495186
# 5 PCOS (Venkatesh)        BMI     BMI PCOS (Venkatesh)             Weighted mode   16 0.05597081 0.02438383 0.0365472322

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method     Estimate          SE       CI_low      CI_upp         P
# 1  Egger fixed effects -0.001640664 0.002094909 -0.005746611 0.002465282 0.4335289
# 2 Egger random effects -0.001640664 0.004609740 -0.010675588 0.007394260 0.6390470
# 
# [[1]]$Q
# Method         Q df            P
# 1   Q_ivw 68.400977 15 8.593297e-09
# 2 Q_egger 67.787626 14 4.842000e-09
# 3  Q_diff  0.613351  1 4.335289e-01
# 
# [[1]]$res
# [1] "B"

strict_het <- compute_strict_het(mr_steiger_df)

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "LOR"
mr_steiger_df$units.outcome <- "SD"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al/3_reverse_mr/2_mr/whr/venkatesh/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al/3_reverse_mr/2_mr/whr/venkatesh/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al/3_reverse_mr/2_mr/whr/venkatesh//sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al/3_reverse_mr/2_mr/whr/venkatesh//strict_het_test_before_outlier_extraction.RDS")

#A lot of heterogeneity

mr_post <- remove_outlier(mr_steiger_df, rucker)

saveRDS(mr_post, "output/2_replicating_liu_et_al/3_reverse_mr/2_mr/whr/venkatesh/instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome         exposure                    method nsnp          b          se         pval
# 1 PCOS (Venkatesh)        BMI     BMI PCOS (Venkatesh)                  MR Egger    9 0.01307595 0.034112359 7.128582e-01
# 2 PCOS (Venkatesh)        BMI     BMI PCOS (Venkatesh)           Weighted median    9 0.03876824 0.012676862 2.226791e-03
# 3 PCOS (Venkatesh)        BMI     BMI PCOS (Venkatesh) Inverse variance weighted    9 0.03619503 0.009225006 8.724356e-05
# 4 PCOS (Venkatesh)        BMI     BMI PCOS (Venkatesh)               Simple mode    9 0.04502708 0.022773025 8.340789e-02
# 5 PCOS (Venkatesh)        BMI     BMI PCOS (Venkatesh)             Weighted mode    9 0.04297446 0.019129407 5.486504e-02

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method    Estimate          SE       CI_low      CI_upp         P
# 1  Egger fixed effects 0.002478536 0.003283147 -0.003956314 0.008913387 0.4502929
# 2 Egger random effects 0.002478536 0.003283147 -0.003956314 0.008913387 0.2251464
# 
# [[1]]$Q
# Method         Q df         P
# 1   Q_ivw 6.5823550  8 0.5822880
# 2 Q_egger 6.0867906  7 0.5296527
# 3  Q_diff 0.4955643  1 0.4814556
# 
# [[1]]$res
# [1] "A"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_post) #heterogeneity removed

#That is quite good:

mr_post$labels <- NA
mr_post$units.exposure <- "LOR"
mr_post$units.outcome <- "SD"

plotio=mr_plots(mr_post)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/venkatesh/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/venkatesh/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/venkatesh/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whr/venkatesh/strict_het_test_after_outlier_extraction.RDS")
