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

liu_pcos_doctor <- fread("output/2_replicating_liu_et_al/3_reverse_mr/1_expo_outcome_df/pcos_finngen_BMI.txt")
liu_pcos_broad <- fread("output/2_replicating_liu_et_al/3_reverse_mr/1_expo_outcome_df/pcos_broad_BMI.txt")
liu_pcos_consortium <- fread("output/2_replicating_liu_et_al/3_reverse_mr/1_expo_outcome_df/pcos_consortium_BMI.txt")
liu_pcos_meta_analysis <- fread("output/2_replicating_liu_et_al/3_reverse_mr/1_expo_outcome_df/pcos_day_BMI.txt")
liu_pcos_venkatesh <- fread("output/2_replicating_liu_et_al/3_reverse_mr/1_expo_outcome_df/pcos_venkatesh_BMI.txt")
liu_pcos_adj_age <- fread("output/2_replicating_liu_et_al/3_reverse_mr/1_expo_outcome_df/pcos_tyrmi_BMI.txt")
liu_pcos_adj_age_bmi <- fread("output/2_replicating_liu_et_al/3_reverse_mr/1_expo_outcome_df/pcosadjbmi_tyrmi_BMI.txt")

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

mr_df$id.outcome <- "BMI"
mr_df$id.exposure <- "PCOS (meta-analysis)"
mr_df$outcome <- "BMI"
mr_df$exposure <- "PCOS (meta-analysis)"
mr_df$units.outcome <- "SD"
mr_df$units.exposure <- "LOR"

mr_steiger_df <- TwoSampleMR::steiger_filtering(mr_df)
mr_steiger_df <- mr_steiger_df[which(mr_steiger_df$steiger_dir == TRUE),] 
mr_steiger_df$mr_keep <- TRUE

#Let's save the data:

dir.create("output/2_replicating_liu_et_al/3_reverse_mr/2_mr/")
dir.create("output/2_replicating_liu_et_al/3_reverse_mr/2_mr/bmi/")
dir.create("output/2_replicating_liu_et_al/3_reverse_mr/2_mr/bmi/meta_analysis/")
saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al/3_reverse_mr/2_mr/bmi/meta_analysis/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome             exposure                    method nsnp             b          se      pval
# 1 PCOS (meta-analysis)        BMI     BMI PCOS (meta-analysis)                  MR Egger   14  0.0388024839 0.061481375 0.5397936
# 2 PCOS (meta-analysis)        BMI     BMI PCOS (meta-analysis)           Weighted median   14  0.0000000000 0.009680109 1.0000000
# 3 PCOS (meta-analysis)        BMI     BMI PCOS (meta-analysis) Inverse variance weighted   14 -0.0064073274 0.012986528 0.6217422
# 4 PCOS (meta-analysis)        BMI     BMI PCOS (meta-analysis)               Simple mode   14  0.0027493989 0.014608394 0.8536207
# 5 PCOS (meta-analysis)        BMI     BMI PCOS (meta-analysis)             Weighted mode   14  0.0007678908 0.014816960 0.9594559

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method     Estimate          SE      CI_low      CI_upp         P
# 1  Egger fixed effects -0.005831021 0.003616782 -0.01291978 0.001257740 0.1069157
# 2 Egger random effects -0.005831021 0.007744520 -0.02101000 0.009347959 0.7742517
# 
# [[1]]$Q
# Method        Q df            P
# 1   Q_ivw 57.61988 13 1.394081e-07
# 2 Q_egger 55.02065 12 1.794560e-07
# 3  Q_diff  2.59923  1 1.069157e-01
# 
# [[1]]$res
# [1] "B"

strict_het <- compute_strict_het(mr_steiger_df)

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "LOR"
mr_steiger_df$units.outcome <- "SD"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al/3_reverse_mr/2_mr/bmi/meta_analysis/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al/3_reverse_mr/2_mr/bmi/meta_analysis/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al/3_reverse_mr/2_mr/bmi/meta_analysis//sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al/3_reverse_mr/2_mr/bmi/meta_analysis//strict_het_test_before_outlier_extraction.RDS")

#A lot of heterogeneity

mr_post <- remove_outlier(mr_steiger_df, rucker)

saveRDS(mr_post, "output/2_replicating_liu_et_al/3_reverse_mr/2_mr/bmi/meta_analysis/instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome             exposure                    method nsnp             b          se      pval
# 1 PCOS (meta-analysis)        BMI     BMI PCOS (meta-analysis)                  MR Egger   11 -0.0663870478 0.042505809 0.1527617
# 2 PCOS (meta-analysis)        BMI     BMI PCOS (meta-analysis)           Weighted median   11  0.0000000000 0.009485380 1.0000000
# 3 PCOS (meta-analysis)        BMI     BMI PCOS (meta-analysis) Inverse variance weighted   11 -0.0003096461 0.007765933 0.9681949
# 4 PCOS (meta-analysis)        BMI     BMI PCOS (meta-analysis)               Simple mode   11 -0.0001420672 0.015860213 0.9930293
# 5 PCOS (meta-analysis)        BMI     BMI PCOS (meta-analysis)             Weighted mode   11  0.0004866622 0.014540585 0.9739590

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method    Estimate          SE       CI_low     CI_upp          P
# 1  Egger fixed effects 0.008196544 0.005110261 -0.001819384 0.01821247 0.10872766
# 2 Egger random effects 0.008196544 0.005195455 -0.001986360 0.01837945 0.05732446
# 
# [[1]]$Q
# Method         Q df         P
# 1   Q_ivw 11.875198 10 0.2934959
# 2 Q_egger  9.302580  9 0.4098253
# 3  Q_diff  2.572618  1 0.1087277
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/meta_analysis/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/meta_analysis/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/meta_analysis/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/meta_analysis/strict_het_test_after_outlier_extraction.RDS")

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

dir.create("output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/")
dir.create("output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/doctor")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/doctor/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp           b          se      pval
# 1    exposure    outcome outcome exposure                  MR Egger    6 0.031739719 0.028405917 0.3264193
# 2    exposure    outcome outcome exposure           Weighted median    6 0.003428967 0.007442574 0.6449971
# 3    exposure    outcome outcome exposure Inverse variance weighted    6 0.002740462 0.011663700 0.8142425
# 4    exposure    outcome outcome exposure               Simple mode    6 0.001636645 0.009097219 0.8642903
# 5    exposure    outcome outcome exposure             Weighted mode    6 0.002102924 0.010827275 0.8536446

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method     Estimate          SE      CI_low       CI_upp           P
# 1  Egger fixed effects -0.008113028 0.003036424 -0.01406431 -0.002161747 0.007542258
# 2 Egger random effects -0.008113028 0.007280029 -0.02238162  0.006155567 0.867450997
# 
# [[1]]$Q
# Method         Q df            P
# 1   Q_ivw 30.132372  5 1.388961e-05
# 2 Q_egger 22.993309  4 1.270165e-04
# 3  Q_diff  7.139062  1 7.542258e-03
# 
# [[1]]$res
# [1] "D"

strict_het <- compute_strict_het(mr_steiger_df)

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "LOR"
mr_steiger_df$units.outcome <- "SD"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/doctor/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/doctor/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/doctor/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/doctor/strict_het_test_before_outlier_extraction.RDS")

#Biased by pleiotropy  (which is in part due to small number of instruments)

mr_post <- remove_outlier(mr_steiger_df, rucker)

saveRDS(mr_post, "output/2_replicating_liu_et_al/3_reverse_mr/2_mr/bmi/doctor/instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure                    method nsnp           b          se      pval
# 1    exposure    outcome outcome exposure                  MR Egger    5 0.029435408 0.018784391 0.2150998
# 2    exposure    outcome outcome exposure           Weighted median    5 0.003679437 0.007001111 0.5992006
# 3    exposure    outcome outcome exposure Inverse variance weighted    5 0.011434944 0.008267327 0.1666192
# 4    exposure    outcome outcome exposure               Simple mode    5 0.002701910 0.008526859 0.7671908
# 5    exposure    outcome outcome exposure             Weighted mode    5 0.002825789 0.010057436 0.7926714

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method     Estimate          SE      CI_low       CI_upp         P
# 1  Egger fixed effects -0.005256211 0.003122085 -0.01137539 0.0008629627 0.0922671
# 2 Egger random effects -0.005256211 0.004943943 -0.01494616 0.0044337389 0.8561457
# 
# [[1]]$Q
# Method         Q df          P
# 1   Q_ivw 10.357154  4 0.03482258
# 2 Q_egger  7.522787  3 0.05697585
# 3  Q_diff  2.834368  1 0.09226710
# 
# [[1]]$res
# [1] "B"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_post) #heterogeneity removed

#That is quite good:

mr_post$labels <- NA
mr_post$units.exposure <- "LOR"
mr_post$units.outcome <- "SD"

plotio=mr_plots(mr_post)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/doctor/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/doctor/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/doctor/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/doctor/strict_het_test_after_outlier_extraction.RDS")

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

dir.create("output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/broad")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/broad/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp            b          se       pval
# 1    exposure    outcome outcome exposure                  MR Egger   12  0.060232823 0.029171900 0.06586033
# 2    exposure    outcome outcome exposure           Weighted median   12  0.009977684 0.009648502 0.30108128
# 3    exposure    outcome outcome exposure Inverse variance weighted   12  0.011404138 0.011952469 0.34002084
# 4    exposure    outcome outcome exposure               Simple mode   12 -0.009223191 0.017984328 0.61820064
# 5    exposure    outcome outcome exposure             Weighted mode   12  0.005317060 0.016776108 0.75722136

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method     Estimate          SE      CI_low        CI_upp            P
# 1  Egger fixed effects -0.008926693 0.002611376 -0.01404490 -0.0038084902 0.0006299353
# 2 Egger random effects -0.008926693 0.004947694 -0.01862399  0.0007706087 0.9644010306
# 
# [[1]]$Q
# Method        Q df            P
# 1   Q_ivw 47.58308 11 1.693938e-06
# 2 Q_egger 35.89772 10 8.765159e-05
# 3  Q_diff 11.68537  1 6.299353e-04
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/broad/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/broad/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/broad/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/broad/strict_het_test_before_outlier_extraction.RDS")

#Biased by pleiotropy  (which is in part due to small number of instruments)

mr_post <- remove_outlier(mr_steiger_df, rucker)

saveRDS(mr_post, "output/2_replicating_liu_et_al/3_reverse_mr/2_mr/bmi/broad/instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure                    method nsnp            b          se      pval
# 1    exposure    outcome outcome exposure                  MR Egger    9 -0.003028638 0.032135261 0.9275543
# 2    exposure    outcome outcome exposure           Weighted median    9  0.006351046 0.009907238 0.5214894
# 3    exposure    outcome outcome exposure Inverse variance weighted    9  0.002615779 0.009170402 0.7754592
# 4    exposure    outcome outcome exposure               Simple mode    9 -0.008259913 0.015290171 0.6037608
# 5    exposure    outcome outcome exposure             Weighted mode    9  0.003095574 0.012974281 0.8174193

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method     Estimate          SE       CI_low      CI_upp         P
# 1  Egger fixed effects 0.0009207779 0.003524073 -0.005986278 0.007827834 0.7938748
# 2 Egger random effects 0.0009207779 0.004993585 -0.008866469 0.010708025 0.4268529
# 
# [[1]]$Q
# Method           Q df          P
# 1   Q_ivw 14.12334328  8 0.07860663
# 2 Q_egger 14.05507484  7 0.05021042
# 3  Q_diff  0.06826844  1 0.79387482
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/broad/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/broad/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/broad/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/broad/strict_het_test_after_outlier_extraction.RDS")

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

dir.create("output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/consortium")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/consortium/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp           b         se      pval
# 1    exposure    outcome outcome exposure                  MR Egger    5 -0.01602405 0.18072218 0.9349342
# 2    exposure    outcome outcome exposure           Weighted median    5  0.01522455 0.03267190 0.6412276
# 3    exposure    outcome outcome exposure Inverse variance weighted    5  0.03326845 0.06916669 0.6305240
# 4    exposure    outcome outcome exposure               Simple mode    5  0.05976396 0.05361307 0.3274231
# 5    exposure    outcome outcome exposure             Weighted mode    5  0.01944331 0.03202716 0.5765660

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method    Estimate          SE      CI_low     CI_upp         P
# 1  Egger fixed effects 0.003519555 0.003324064 -0.00299549 0.01003460 0.2896859
# 2 Egger random effects 0.003519555 0.011616987 -0.01924932 0.02628843 0.3809578
# 
# [[1]]$Q
# Method         Q df            P
# 1   Q_ivw 37.762271  4 1.254493e-07
# 2 Q_egger 36.641190  3 5.480157e-08
# 3  Q_diff  1.121081  1 2.896859e-01
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/consortium/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/consortium/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/consortium/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/consortium/strict_het_test_before_outlier_extraction.RDS")

#Very small heterogeneity, overall. However, the interecept is the intercept

#STEP 5: removing outliers!

mr_post <- remove_outlier(mr_steiger_df, rucker)
  
saveRDS(mr_post, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/consortium/instruments_after_outlier_removal_df.RDS")
  
#STEP 6: let's run the results for real now:
  
mr_res <- TwoSampleMR::mr(mr_post)
  
# id.exposure id.outcome outcome exposure                    method nsnp           b         se      pval
# 1    exposure    outcome outcome exposure                  MR Egger    3 -0.07272813 0.07010874 0.4883267
# 2    exposure    outcome outcome exposure           Weighted median    3  0.01540990 0.03056754 0.6141727
# 3    exposure    outcome outcome exposure Inverse variance weighted    3  0.03005371 0.03782842 0.4269194
# 4    exposure    outcome outcome exposure               Simple mode    3  0.00834072 0.04011125 0.8545285
# 5    exposure    outcome outcome exposure             Weighted mode    3  0.00834072 0.03352585 0.8267431

# #Let's check the plot:
#  
rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method    Estimate          SE       CI_low     CI_upp          P
# 1  Egger fixed effects 0.009870945 0.005828482 -0.001552669 0.02129456 0.09034690
# 2 Egger random effects 0.009870945 0.006159305 -0.002201071 0.02194296 0.05451073
# 
# [[1]]$Q
# Method        Q df         P
# 1   Q_ivw 3.984923  2 0.1363594
# 2 Q_egger 1.116741  1 0.2906212
# 3  Q_diff 2.868182  1 0.0903469
# 
# [[1]]$res
# [1] "A"

# #Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_post) #they disagree, but still
  
#That is quite good:
  
mr_post$labels <- NA
mr_post$units.exposure <- "LOR"
mr_post$units.outcome <- "SD"
  
plotio=mr_plots(mr_post)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/consortium/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/consortium/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/consortium/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/consortium/strict_het_test_after_outlier_extraction.RDS")

###############################################################
#Let's first analyse the data for the whole set for bmi - PCOS# age adjusted
###############################################################

#STEP 1: let's use only those that have F>10

mr_df <- liu_pcos_adj_age[which(((liu_pcos_adj_age$beta.exposure^2)/(liu_pcos_adj_age$se.exposure^2)) > 10),] #194

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

dir.create("output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/age_adjusted")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/age_adjusted/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp             b         se      pval
# 1    exposure    outcome outcome exposure                  MR Egger    6  0.0892515433 0.05543207 0.1826620
# 2    exposure    outcome outcome exposure           Weighted median    6  0.0001966058 0.01106880 0.9858286
# 3    exposure    outcome outcome exposure Inverse variance weighted    6 -0.0022846542 0.01706995 0.8935285
# 4    exposure    outcome outcome exposure               Simple mode    6 -0.0020610916 0.01551033 0.8994667
# 5    exposure    outcome outcome exposure             Weighted mode    6  0.0002030604 0.01562073 0.9901310

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[[1]]$intercept
# Method    Estimate          SE      CI_low       CI_upp            P
# 1  Egger fixed effects -0.01899707 0.005264378 -0.02931506 -0.008679076 0.0003078463
# 2 Egger random effects -0.01899707 0.011103467 -0.04075946  0.002765330 0.9564514084
# 
# [[1]]$Q
# Method        Q df            P
# 1   Q_ivw 30.81643  5 0.0000101823
# 2 Q_egger 17.79439  4 0.0013536577
# 3  Q_diff 13.02204  1 0.0003078463
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/age_adjusted/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/age_adjusted/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/age_adjusted/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/age_adjusted/strict_het_test_before_outlier_extraction.RDS")

#There is a bit of asymmetry. I am almost certain that Wmed and Wmod correct for this, but we are gonna run the analyses again just to be on the safe side.

mr_post <- remove_outlier(mr_steiger_df, rucker)

saveRDS(mr_post, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/age_adjusted/instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure                    method nsnp             b         se      pval
# 1    exposure    outcome outcome exposure                  MR Egger    4 -0.0319672411 0.04352783 0.5391323
# 2    exposure    outcome outcome exposure           Weighted median    4  0.0003005909 0.01078215 0.9777590
# 3    exposure    outcome outcome exposure Inverse variance weighted    4 -0.0029343102 0.00926933 0.7515766
# 4    exposure    outcome outcome exposure               Simple mode    4  0.0004870692 0.01281205 0.9720629
# 5    exposure    outcome outcome exposure             Weighted mode    4  0.0005995679 0.01180445 0.9626841

# #Let's check the plot:
#  
rucker <- TwoSampleMR::mr_rucker(mr_post)

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

strict_het <- compute_strict_het(mr_post) #they disagree, but still

#That is quite good:

mr_post$labels <- NA
mr_post$units.exposure <- "LOR"
mr_post$units.outcome <- "SD"

plotio=mr_plots(mr_post)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/age_adjusted/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/age_adjusted/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/age_adjusted/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/age_adjusted/strict_het_test_after_outlier_extraction.RDS")

###############################################################
#Let's first analyse the data for the whole set for bmi - PCOS# age and BMI adjusted
###############################################################

#STEP 1: let's use only those that have F>10

mr_df <- liu_pcos_adj_age_bmi[which(((liu_pcos_adj_age_bmi$beta.exposure^2)/(liu_pcos_adj_age_bmi$se.exposure^2)) > 10),] 
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

dir.create("output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/age_and_bmi_adjusted")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/age_and_bmi_adjusted/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp             b         se       pval
# 1    exposure    outcome outcome exposure                  MR Egger    3  0.1454273445 0.11012871 0.41261991
# 2    exposure    outcome outcome exposure           Weighted median    3 -0.0387331827 0.01886354 0.04004000
# 3    exposure    outcome outcome exposure Inverse variance weighted    3 -0.0003185621 0.03894617 0.99347374
# 4    exposure    outcome outcome exposure               Simple mode    3 -0.0439507740 0.01936062 0.15122853
# 5    exposure    outcome outcome exposure             Weighted mode    3 -0.0479462844 0.01589206 0.09454063

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method   Estimate          SE      CI_low      CI_upp            P
# 1  Egger fixed effects -0.0346346 0.007076341 -0.04850397 -0.02076523 9.859515e-07
# 2 Egger random effects -0.0346346 0.025022968 -0.08367872  0.01440952 9.168380e-01
# 
# [[1]]$Q
# Method        Q df            P
# 1   Q_ivw 36.45969  2 1.210263e-08
# 2 Q_egger 12.50432  1 4.060120e-04
# 3  Q_diff 23.95537  1 9.859515e-07
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/age_and_bmi_adjusted/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/age_and_bmi_adjusted/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/age_and_bmi_adjusted/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/age_and_bmi_adjusted/strict_het_test_before_outlier_extraction.RDS")

#Let's remove outliers:

mr_post <- remove_outlier(mr_steiger_df, rucker)

saveRDS(mr_post, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/age_and_bmi_adjusted/instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure     method nsnp           b         se      pval
# 1    exposure    outcome outcome exposure Wald ratio    1 -0.03684211 0.02797784 0.1878951

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
# ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/age_and_bmi_adjusted/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/age_and_bmi_adjusted/mr_res_after_outlier_extraction.RDS")
# saveRDS(rucker, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/age_and_bmi_adjusted/sensitivity_test_after_outlier_extraction.RDS")
# saveRDS(strict_het, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/age_and_bmi_adjusted/strict_het_test_after_outlier_extraction.RDS")

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

mr_df$id.outcome <- "BMI"
mr_df$id.exposure <- "PCOS (Venkatesh)"
mr_df$outcome <- "BMI"
mr_df$exposure <- "PCOS (Venkatesh)"
mr_df$units.outcome <- "SD"
mr_df$units.exposure <- "LOR"

mr_steiger_df <- TwoSampleMR::steiger_filtering(mr_df)
mr_steiger_df <- mr_steiger_df[which(mr_steiger_df$steiger_dir == TRUE),] 
mr_steiger_df$mr_keep <- TRUE

#Let's save the data:

dir.create("output/2_replicating_liu_et_al/3_reverse_mr/2_mr/")
dir.create("output/2_replicating_liu_et_al/3_reverse_mr/2_mr/bmi/")
dir.create("output/2_replicating_liu_et_al/3_reverse_mr/2_mr/bmi/venkatesh/")
saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al/3_reverse_mr/2_mr/bmi/venkatesh/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome         exposure                    method nsnp             b         se       pval
# 1 PCOS (Venkatesh)        BMI     BMI PCOS (Venkatesh)                  MR Egger   16  0.0682012482 0.03354796 0.06147269
# 2 PCOS (Venkatesh)        BMI     BMI PCOS (Venkatesh)           Weighted median   16  0.0035779473 0.01118075 0.74896104
# 3 PCOS (Venkatesh)        BMI     BMI PCOS (Venkatesh) Inverse variance weighted   16  0.0048879618 0.01427776 0.73208912
# 4 PCOS (Venkatesh)        BMI     BMI PCOS (Venkatesh)               Simple mode   16 -0.0003988849 0.01999619 0.98434775
# 5 PCOS (Venkatesh)        BMI     BMI PCOS (Venkatesh)             Weighted mode   16  0.0060466836 0.01532134 0.69865038

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method     Estimate          SE      CI_low        CI_upp            P
# 1  Egger fixed effects -0.008349065 0.002044412 -0.01235604 -0.0043420918 4.429619e-05
# 2 Egger random effects -0.008349065 0.004080141 -0.01634599 -0.0003521361 9.796350e-01
# 
# [[1]]$Q
# Method        Q df            P
# 1   Q_ivw 72.44022 15 1.636030e-09
# 2 Q_egger 55.76242 14 6.400905e-07
# 3  Q_diff 16.67781  1 4.429619e-05
# 
# [[1]]$res
# [1] "D"

strict_het <- compute_strict_het(mr_steiger_df)

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "LOR"
mr_steiger_df$units.outcome <- "SD"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al/3_reverse_mr/2_mr/bmi/venkatesh/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al/3_reverse_mr/2_mr/bmi/venkatesh/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al/3_reverse_mr/2_mr/bmi/venkatesh//sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al/3_reverse_mr/2_mr/bmi/venkatesh//strict_het_test_before_outlier_extraction.RDS")

#A lot of heterogeneity

mr_post <- remove_outlier(mr_steiger_df, rucker)

saveRDS(mr_post, "output/2_replicating_liu_et_al/3_reverse_mr/2_mr/bmi/venkatesh/instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome         exposure                    method nsnp            b          se      pval
# 1 PCOS (Venkatesh)        BMI     BMI PCOS (Venkatesh)                  MR Egger   13  0.001374263 0.029138127 0.9632281
# 2 PCOS (Venkatesh)        BMI     BMI PCOS (Venkatesh)           Weighted median   13  0.002265390 0.010755838 0.8331841
# 3 PCOS (Venkatesh)        BMI     BMI PCOS (Venkatesh) Inverse variance weighted   13 -0.002019875 0.009593903 0.8332483
# 4 PCOS (Venkatesh)        BMI     BMI PCOS (Venkatesh)               Simple mode   13  0.004674108 0.019537308 0.8149567
# 5 PCOS (Venkatesh)        BMI     BMI PCOS (Venkatesh)             Weighted mode   13  0.004674108 0.015148956 0.7629597

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method      Estimate          SE       CI_low      CI_upp         P
# 1  Egger fixed effects -0.0004162858 0.002498082 -0.005312437 0.004479865 0.8676516
# 2 Egger random effects -0.0004162858 0.003356089 -0.006994099 0.006161528 0.5493578
# 
# [[1]]$Q
# Method          Q df          P
# 1   Q_ivw 19.8816842 12 0.06935746
# 2 Q_egger 19.8539146 11 0.04738481
# 3  Q_diff  0.0277696  1 0.86765163
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/venkatesh/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/venkatesh/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/venkatesh/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/bmi/venkatesh/strict_het_test_after_outlier_extraction.RDS")


