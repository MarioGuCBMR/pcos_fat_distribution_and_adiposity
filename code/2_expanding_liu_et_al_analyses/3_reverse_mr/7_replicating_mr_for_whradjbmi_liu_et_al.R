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

liu_pcos_doctor <- fread("output/2_replicating_liu_et_al/3_reverse_mr/1_expo_outcome_df/pcos_finngen_WHRadjBMI.txt")
liu_pcos_broad <- fread("output/2_replicating_liu_et_al/3_reverse_mr/1_expo_outcome_df/pcos_broad_WHRadjBMI.txt")
liu_pcos_consortium <- fread("output/2_replicating_liu_et_al/3_reverse_mr/1_expo_outcome_df/pcos_consortium_WHRadjBMI.txt")
liu_pcos_meta_analysis <- fread("output/2_replicating_liu_et_al/3_reverse_mr/1_expo_outcome_df/pcos_day_WHRadjBMI.txt")
liu_pcos_venkatesh <- fread("output/2_replicating_liu_et_al/3_reverse_mr/1_expo_outcome_df/pcos_venkatesh_WHRadjBMI.txt")
liu_pcos_adj_age <- fread("output/2_replicating_liu_et_al/3_reverse_mr/1_expo_outcome_df/pcos_tyrmi_WHRadjBMI.txt")
liu_pcos_adj_age_bmi <- fread("output/2_replicating_liu_et_al/3_reverse_mr/1_expo_outcome_df/pcosadjbmi_tyrmi_WHRadjBMI.txt")

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

mr_df$id.outcome <- "WHRadjBMI"
mr_df$id.exposure <- "PCOS (meta-analysis)"
mr_df$outcome <- "WHRadjBMI"
mr_df$exposure <- "PCOS (meta-analysis)"
mr_df$units.outcome <- "SD"
mr_df$units.exposure <- "LOR"

mr_steiger_df <- TwoSampleMR::steiger_filtering(mr_df)
mr_steiger_df <- mr_steiger_df[which(mr_steiger_df$steiger_dir == TRUE),] 
mr_steiger_df$mr_keep <- TRUE

#Let's save the data:

dir.create("output/2_replicating_liu_et_al/3_reverse_mr/2_mr/")
dir.create("output/2_replicating_liu_et_al/3_reverse_mr/2_mr/whradjbmi/")
dir.create("output/2_replicating_liu_et_al/3_reverse_mr/2_mr/whradjbmi/meta_analysis/")
saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al/3_reverse_mr/2_mr/whradjbmi/meta_analysis/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome   outcome             exposure                    method nsnp          b         se      pval
# 1 PCOS (meta-analysis)  WHRadjBMI WHRadjBMI PCOS (meta-analysis)                  MR Egger   14 0.06411302 0.05538912 0.2695880
# 2 PCOS (meta-analysis)  WHRadjBMI WHRadjBMI PCOS (meta-analysis)           Weighted median   14 0.01301332 0.01073755 0.2255335
# 3 PCOS (meta-analysis)  WHRadjBMI WHRadjBMI PCOS (meta-analysis) Inverse variance weighted   14 0.01596743 0.01185029 0.1778421
# 4 PCOS (meta-analysis)  WHRadjBMI WHRadjBMI PCOS (meta-analysis)               Simple mode   14 0.01416071 0.02065678 0.5050553
# 5 PCOS (meta-analysis)  WHRadjBMI WHRadjBMI PCOS (meta-analysis)             Weighted mode   14 0.01300276 0.01735662 0.4671052

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method     Estimate          SE      CI_low      CI_upp          P
# 1  Egger fixed effects -0.006214701 0.003715168 -0.01349630 0.001066895 0.09436838
# 2 Egger random effects -0.006214701 0.006981440 -0.01989807 0.007468669 0.81331397
# 
# [[1]]$Q
# Method         Q df            P
# 1   Q_ivw 45.173674 13 1.957636e-05
# 2 Q_egger 42.375443 12 2.879163e-05
# 3  Q_diff  2.798231  1 9.436838e-02
# 
# [[1]]$res
# [1] "B"

strict_het <- compute_strict_het(mr_steiger_df)

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "LOR"
mr_steiger_df$units.outcome <- "SD"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al/3_reverse_mr/2_mr/whradjbmi/meta_analysis/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al/3_reverse_mr/2_mr/whradjbmi/meta_analysis/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al/3_reverse_mr/2_mr/whradjbmi/meta_analysis//sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al/3_reverse_mr/2_mr/whradjbmi/meta_analysis//strict_het_test_before_outlier_extraction.RDS")

#A lot of heterogeneity

mr_post <- remove_outlier(mr_steiger_df, rucker)

saveRDS(mr_post, "output/2_replicating_liu_et_al/3_reverse_mr/2_mr/whradjbmi/meta_analysis/instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome   outcome             exposure                    method nsnp          b          se       pval
# 1 PCOS (meta-analysis)  WHRadjBMI WHRadjBMI PCOS (meta-analysis)                  MR Egger   10 0.03184661 0.059945784 0.60966993
# 2 PCOS (meta-analysis)  WHRadjBMI WHRadjBMI PCOS (meta-analysis)           Weighted median   10 0.01325474 0.011266259 0.23939575
# 3 PCOS (meta-analysis)  WHRadjBMI WHRadjBMI PCOS (meta-analysis) Inverse variance weighted   10 0.02064184 0.009210753 0.02502228
# 4 PCOS (meta-analysis)  WHRadjBMI WHRadjBMI PCOS (meta-analysis)               Simple mode   10 0.01074087 0.017430652 0.55302266
# 5 PCOS (meta-analysis)  WHRadjBMI WHRadjBMI PCOS (meta-analysis)             Weighted mode   10 0.01015068 0.016033792 0.54243090

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method     Estimate          SE      CI_low      CI_upp         P
# 1  Egger fixed effects -0.001370354 0.005599029 -0.01234425 0.009603541 0.8066511
# 2 Egger random effects -0.001370354 0.007233850 -0.01554844 0.012807731 0.5751246
# 
# [[1]]$Q
# Method           Q df         P
# 1   Q_ivw 13.41366197  9 0.1447643
# 2 Q_egger 13.35376011  8 0.1002437
# 3  Q_diff  0.05990186  1 0.8066511
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/meta_analysis/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/meta_analysis/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/meta_analysis/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/meta_analysis/strict_het_test_after_outlier_extraction.RDS")

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

dir.create("output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/")
dir.create("output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/doctor")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/doctor/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp           b          se      pval
# 1    exposure    outcome outcome exposure                  MR Egger    6 0.030536385 0.016218966 0.1328572
# 2    exposure    outcome outcome exposure           Weighted median    6 0.007230627 0.007458245 0.3323053
# 3    exposure    outcome outcome exposure Inverse variance weighted    6 0.002268867 0.008075979 0.7787563
# 4    exposure    outcome outcome exposure               Simple mode    6 0.006319403 0.013556014 0.6606974
# 5    exposure    outcome outcome exposure             Weighted mode    6 0.010355703 0.011416686 0.4059626

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method     Estimate          SE      CI_low        CI_upp         P
# 1  Egger fixed effects -0.007967601 0.003116601 -0.01407603 -0.0018591756 0.0105730
# 2 Egger random effects -0.007967601 0.004183558 -0.01616722  0.0002320226 0.9715777
# 
# [[1]]$Q
# Method         Q df          P
# 1   Q_ivw 13.743287  5 0.01732501
# 2 Q_egger  7.207577  4 0.12531699
# 3  Q_diff  6.535710  1 0.01057300
# 
# [[1]]$res
# [1] "C"

strict_het <- compute_strict_het(mr_steiger_df)

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "LOR"
mr_steiger_df$units.outcome <- "SD"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/doctor/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/doctor/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/doctor/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/doctor/strict_het_test_before_outlier_extraction.RDS")

#Biased by pleiotropy  (which is in part due to small number of instruments)

mr_post <- remove_outlier(mr_steiger_df, rucker)

# saveRDS(mr_post, "output/2_replicating_liu_et_al/3_reverse_mr/2_mr/whradjbmi/doctor/instruments_after_outlier_removal_df.RDS")
# 
# #STEP 6: let's run the results for real now:
# 
# mr_res <- TwoSampleMR::mr(mr_post)
# 
# # id.exposure id.outcome outcome exposure                    method nsnp           b          se      pval
# # 1    exposure    outcome outcome exposure                  MR Egger    5 0.029435408 0.018784391 0.2150998
# # 2    exposure    outcome outcome exposure           Weighted median    5 0.003679437 0.007001111 0.5992006
# # 3    exposure    outcome outcome exposure Inverse variance weighted    5 0.011434944 0.008267327 0.1666192
# # 4    exposure    outcome outcome exposure               Simple mode    5 0.002701910 0.008526859 0.7671908
# # 5    exposure    outcome outcome exposure             Weighted mode    5 0.002825789 0.010057436 0.7926714
# 
# #Let's check the plot:
# 
# rucker <- TwoSampleMR::mr_rucker(mr_post)
# 
# # [[1]]$intercept
# # Method     Estimate          SE      CI_low       CI_upp         P
# # 1  Egger fixed effects -0.005256211 0.003122085 -0.01137539 0.0008629627 0.0922671
# # 2 Egger random effects -0.005256211 0.004943943 -0.01494616 0.0044337389 0.8561457
# # 
# # [[1]]$Q
# # Method         Q df          P
# # 1   Q_ivw 10.357154  4 0.03482258
# # 2 Q_egger  7.522787  3 0.05697585
# # 3  Q_diff  2.834368  1 0.09226710
# # 
# # [[1]]$res
# # [1] "B"
# 
# #Here there is heterogeneity and pleitropy
# 
# strict_het <- compute_strict_het(mr_post) #heterogeneity removed
# 
# #That is quite good:
# 
# mr_post$labels <- NA
# mr_post$units.exposure <- "LOR"
# mr_post$units.outcome <- "SD"
# 
# plotio=mr_plots(mr_post)
# ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/doctor/plots_after_outlier_extraction.svg", height=16, width=16, units="in")
# 
# saveRDS(mr_res, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/doctor/mr_res_after_outlier_extraction.RDS")
# saveRDS(rucker, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/doctor/sensitivity_test_after_outlier_extraction.RDS")
# saveRDS(strict_het, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/doctor/strict_het_test_after_outlier_extraction.RDS")

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

dir.create("output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/broad")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/broad/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp          b          se       pval
# 1    exposure    outcome outcome exposure                  MR Egger   12 0.02861134 0.032996611 0.40621307
# 2    exposure    outcome outcome exposure           Weighted median   12 0.01026801 0.009392994 0.27432542
# 3    exposure    outcome outcome exposure Inverse variance weighted   12 0.02284465 0.011802177 0.05291297
# 4    exposure    outcome outcome exposure               Simple mode   12 0.01362893 0.015679934 0.40330728
# 5    exposure    outcome outcome exposure             Weighted mode   12 0.01430709 0.012334482 0.27063589

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method     Estimate          SE       CI_low      CI_upp         P
# 1  Egger fixed effects -0.001056677 0.002662037 -0.006274174 0.004160821 0.6914096
# 2 Egger random effects -0.001056677 0.005606304 -0.012044831 0.009931477 0.5747498
# 
# [[1]]$Q
# Method          Q df            P
# 1   Q_ivw 44.5107664 11 5.913801e-06
# 2 Q_egger 44.3532028 10 2.843974e-06
# 3  Q_diff  0.1575636  1 6.914096e-01
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/broad/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/broad/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/broad/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/broad/strict_het_test_before_outlier_extraction.RDS")

#Biased by pleiotropy  (which is in part due to small number of instruments)

mr_post <- remove_outlier(mr_steiger_df, rucker)

saveRDS(mr_post, "output/2_replicating_liu_et_al/3_reverse_mr/2_mr/whradjbmi/broad/instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure                    method nsnp           b          se        pval
# 1    exposure    outcome outcome exposure                  MR Egger    9 0.018021051 0.020499696 0.408514525
# 2    exposure    outcome outcome exposure           Weighted median    9 0.010690188 0.009897735 0.280113655
# 3    exposure    outcome outcome exposure Inverse variance weighted    9 0.023019684 0.007734691 0.002918818
# 4    exposure    outcome outcome exposure               Simple mode    9 0.005310836 0.015491301 0.740562705
# 5    exposure    outcome outcome exposure             Weighted mode    9 0.007510294 0.014874338 0.627230037

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method     Estimate          SE       CI_low      CI_upp         P
# 1  Egger fixed effects 0.0009785374 0.002904909 -0.004714979 0.006672053 0.7362251
# 2 Egger random effects 0.0009785374 0.003675675 -0.006225654 0.008182729 0.3950350
# 
# [[1]]$Q
# Method          Q df         P
# 1   Q_ivw 11.3209377  8 0.1841690
# 2 Q_egger 11.2074654  7 0.1298220
# 3  Q_diff  0.1134723  1 0.7362251
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/broad/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/broad/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/broad/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/broad/strict_het_test_after_outlier_extraction.RDS")

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

dir.create("output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/consortium")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/consortium/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp          b         se        pval
# 1    exposure    outcome outcome exposure                  MR Egger    5 0.16518282 0.06246500 0.077363346
# 2    exposure    outcome outcome exposure           Weighted median    5 0.10151805 0.03085586 0.001001612
# 3    exposure    outcome outcome exposure Inverse variance weighted    5 0.07470022 0.03163906 0.018225086
# 4    exposure    outcome outcome exposure               Simple mode    5 0.10186500 0.03916802 0.060002753
# 5    exposure    outcome outcome exposure             Weighted mode    5 0.10802302 0.03527188 0.037563721

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method     Estimate          SE      CI_low       CI_upp         P
# 1  Egger fixed effects -0.006374721 0.003456005 -0.01314837 0.0003989253 0.0651053
# 2 Egger random effects -0.006374721 0.003974950 -0.01416548 0.0014160374 0.9456125
# 
# [[1]]$Q
# Method        Q df         P
# 1   Q_ivw 7.370892  4 0.1175390
# 2 Q_egger 3.968585  3 0.2648764
# 3  Q_diff 3.402307  1 0.0651053
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/consortium/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/consortium/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/consortium/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/consortium/strict_het_test_before_outlier_extraction.RDS")

#Very small heterogeneity, overall. However, the interecept is the intercept

#STEP 5: removing outliers!

mr_post <- remove_outlier(mr_steiger_df, rucker)
  
saveRDS(mr_post, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/consortium/instruments_after_outlier_removal_df.RDS")
  
#STEP 6: let's run the results for real now:
  
mr_res <- TwoSampleMR::mr(mr_post)
  
# id.exposure id.outcome outcome exposure                    method nsnp         b         se         pval
# 1    exposure    outcome outcome exposure                  MR Egger    4 0.1297249 0.05750180 1.527121e-01
# 2    exposure    outcome outcome exposure           Weighted median    4 0.1092953 0.03225352 7.024219e-04
# 3    exposure    outcome outcome exposure Inverse variance weighted    4 0.1008112 0.02542370 7.332009e-05
# 4    exposure    outcome outcome exposure               Simple mode    4 0.1094561 0.04049217 7.359002e-02
# 5    exposure    outcome outcome exposure             Weighted mode    4 0.1115777 0.03967383 6.716053e-02

# #Let's check the plot:
#  
rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method     Estimate         SE       CI_low     CI_upp         P
# 1  Egger fixed effects -0.002289636 0.00192808 -0.006068603 0.00148933 0.2350219
# 2 Egger random effects -0.002289636 0.00192808 -0.006068603 0.00148933 0.8824891
# 
# [[1]]$Q
# Method         Q df         P
# 1   Q_ivw 0.7599916  3 0.8590106
# 2 Q_egger 0.4457157  2 0.8002286
# 3  Q_diff 0.3142759  1 0.5750681
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/consortium/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/consortium/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/consortium/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/consortium/strict_het_test_after_outlier_extraction.RDS")

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

dir.create("output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/age_adjusted")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/age_adjusted/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp            b         se      pval
# 1    exposure    outcome outcome exposure                  MR Egger    6 0.0661536915 0.06469753 0.3643464
# 2    exposure    outcome outcome exposure           Weighted median    6 0.0009234219 0.01102643 0.9332582
# 3    exposure    outcome outcome exposure Inverse variance weighted    6 0.0071368273 0.01687239 0.6723036
# 4    exposure    outcome outcome exposure               Simple mode    6 0.0078448821 0.01468276 0.6160322
# 5    exposure    outcome outcome exposure             Weighted mode    6 0.0050094695 0.01360558 0.7278038

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method    Estimate         SE      CI_low       CI_upp          P
# 1  Egger fixed effects -0.01223849 0.00535247 -0.02272913 -0.001747838 0.02222433
# 2 Egger random effects -0.01223849 0.01294202 -0.03760437  0.013127401 0.82783382
# 
# [[1]]$Q
# Method         Q df            P
# 1   Q_ivw 28.614113  5 2.759984e-05
# 2 Q_egger 23.385976  4 1.060141e-04
# 3  Q_diff  5.228137  1 2.222433e-02
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/age_adjusted/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/age_adjusted/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/age_adjusted/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/age_adjusted/strict_het_test_before_outlier_extraction.RDS")

#There is a bit of asymmetry. I am almost certain that Wmed and Wmod correct for this, but we are gonna run the analyses again just to be on the safe side.

mr_post <- remove_outlier(mr_steiger_df, rucker)

saveRDS(mr_post, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/age_adjusted/instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure                    method nsnp             b         se      pval
# 1    exposure    outcome outcome exposure                  MR Egger    4  0.1151800851 0.09795721 0.3606795
# 2    exposure    outcome outcome exposure           Weighted median    4  0.0001514778 0.01125955 0.9892662
# 3    exposure    outcome outcome exposure Inverse variance weighted    4  0.0120171283 0.01861620 0.5185902
# 4    exposure    outcome outcome exposure               Simple mode    4  0.0015111349 0.01149289 0.9037146
# 5    exposure    outcome outcome exposure             Weighted mode    4 -0.0007556944 0.01131452 0.9509511

# #Let's check the plot:
#  
rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method    Estimate          SE      CI_low       CI_upp          P
# 1  Egger fixed effects -0.01861458 0.008651671 -0.03557155 -0.001657621 0.03143207
# 2 Egger random effects -0.01861458 0.017368484 -0.05265619  0.015427019 0.85808269
# 
# [[1]]$Q
# Method         Q df           P
# 1   Q_ivw 12.689557  3 0.005358398
# 2 Q_egger  8.060349  2 0.017771228
# 3  Q_diff  4.629208  1 0.031432069
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/age_adjusted/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/age_adjusted/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/age_adjusted/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/age_adjusted/strict_het_test_after_outlier_extraction.RDS")

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

dir.create("output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/age_and_bmi_adjusted")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/age_and_bmi_adjusted/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp          b         se      pval
# 1    exposure    outcome outcome exposure                  MR Egger    3 0.12681270 0.05684047 0.2682562
# 2    exposure    outcome outcome exposure           Weighted median    3 0.01881476 0.01211009 0.1202696
# 3    exposure    outcome outcome exposure Inverse variance weighted    3 0.02600165 0.02485034 0.2954087
# 4    exposure    outcome outcome exposure               Simple mode    3 0.02868414 0.02430103 0.3592212
# 5    exposure    outcome outcome exposure             Weighted mode    3 0.01363807 0.01248033 0.3885652

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method    Estimate          SE      CI_low       CI_upp           P
# 1  Egger fixed effects -0.02397934 0.007356621 -0.03839806 -0.009560632 0.001115854
# 2 Egger random effects -0.02397934 0.012925422 -0.04931271  0.001354018 0.968216718
# 
# [[1]]$Q
# Method         Q df           P
# 1   Q_ivw 13.711700  2 0.001053276
# 2 Q_egger  3.086972  1 0.078921535
# 3  Q_diff 10.624728  1 0.001115854
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/age_and_bmi_adjusted/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/age_and_bmi_adjusted/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/age_and_bmi_adjusted/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/age_and_bmi_adjusted/strict_het_test_before_outlier_extraction.RDS")

#Let's remove outliers:

mr_post <- remove_outlier(mr_steiger_df, rucker)

# saveRDS(mr_post, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/age_and_bmi_adjusted/instruments_after_outlier_removal_df.RDS")
# 
# #STEP 6: let's run the results for real now:
# 
# mr_res <- TwoSampleMR::mr(mr_post)
# 
# # id.exposure id.outcome outcome exposure     method nsnp           b         se      pval
# # 1    exposure    outcome outcome exposure Wald ratio    1 -0.03684211 0.02797784 0.1878951
# 
# # #Let's check the plot:
# #  
# #rucker <- TwoSampleMR::mr_rucker(mr_post)
# 
# # [[1]]$intercept
# # Method   Estimate          SE       CI_low     CI_upp         P
# # 1  Egger fixed effects 0.00571065 0.006517913 -0.007064224 0.01848552 0.3809501
# # 2 Egger random effects 0.00571065 0.006517913 -0.007064224 0.01848552 0.1904750
# # 
# # [[1]]$Q
# # Method         Q df         P
# # 1   Q_ivw 1.6801867  3 0.6413477
# # 2 Q_egger 1.2141684  2 0.5449375
# # 3  Q_diff 0.4660183  1 0.4948247
# # 
# # [[1]]$res
# # [1] "A"
# 
# # #Here there is heterogeneity and pleitropy
# 
# #strict_het <- compute_strict_het(mr_post) #they disagree, but still
# 
# #That is quite good:
# 
# # mr_post$labels <- NA
# # mr_post$units.exposure <- "LOR"
# # mr_post$units.outcome <- "SD"
# # 
# # plotio=mr_plots(mr_post)
# # ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/age_and_bmi_adjusted/plots_after_outlier_extraction.svg", height=16, width=16, units="in")
# 
# saveRDS(mr_res, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/age_and_bmi_adjusted/mr_res_after_outlier_extraction.RDS")
# # saveRDS(rucker, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/age_and_bmi_adjusted/sensitivity_test_after_outlier_extraction.RDS")
# # saveRDS(strict_het, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/age_and_bmi_adjusted/strict_het_test_after_outlier_extraction.RDS")

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
dir.create("output/2_replicating_liu_et_al/3_reverse_mr/2_mr/whradjbmi/")
dir.create("output/2_replicating_liu_et_al/3_reverse_mr/2_mr/whradjbmi/venkatesh/")
saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al/3_reverse_mr/2_mr/whradjbmi/venkatesh/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome         exposure                    method nsnp          b         se       pval
# 1 PCOS (Venkatesh)        BMI     BMI PCOS (Venkatesh)                  MR Egger   17 0.01427323 0.03487295 0.68810973
# 2 PCOS (Venkatesh)        BMI     BMI PCOS (Venkatesh)           Weighted median   17 0.02533037 0.01148414 0.02740627
# 3 PCOS (Venkatesh)        BMI     BMI PCOS (Venkatesh) Inverse variance weighted   17 0.02865222 0.01317074 0.02959673
# 4 PCOS (Venkatesh)        BMI     BMI PCOS (Venkatesh)               Simple mode   17 0.05858128 0.02078062 0.01234662
# 5 PCOS (Venkatesh)        BMI     BMI PCOS (Venkatesh)             Weighted mode   17 0.03951055 0.01872974 0.05100514

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method    Estimate          SE       CI_low      CI_upp         P
# 1  Egger fixed effects 0.001853264 0.002011245 -0.002088704 0.005795231 0.3568149
# 2 Egger random effects 0.001853264 0.004143520 -0.006267885 0.009974413 0.3273408
# 
# [[1]]$Q
# Method          Q df            P
# 1   Q_ivw 64.5139854 16 8.929252e-08
# 2 Q_egger 63.6649132 15 5.848597e-08
# 3  Q_diff  0.8490722  1 3.568149e-01
# 
# [[1]]$res
# [1] "B"

strict_het <- compute_strict_het(mr_steiger_df)

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "LOR"
mr_steiger_df$units.outcome <- "SD"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al/3_reverse_mr/2_mr/whradjbmi/venkatesh/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al/3_reverse_mr/2_mr/whradjbmi/venkatesh/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al/3_reverse_mr/2_mr/whradjbmi/venkatesh//sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al/3_reverse_mr/2_mr/whradjbmi/venkatesh//strict_het_test_before_outlier_extraction.RDS")

#A lot of heterogeneity

mr_post <- remove_outlier(mr_steiger_df, rucker)

saveRDS(mr_post, "output/2_replicating_liu_et_al/3_reverse_mr/2_mr/whradjbmi/venkatesh/instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome         exposure                    method nsnp          b         se         pval
# 1 PCOS (Venkatesh)        BMI     BMI PCOS (Venkatesh)                  MR Egger   11 0.03911003 0.02249238 1.160676e-01
# 2 PCOS (Venkatesh)        BMI     BMI PCOS (Venkatesh)           Weighted median   11 0.04494303 0.01116153 5.658899e-05
# 3 PCOS (Venkatesh)        BMI     BMI PCOS (Venkatesh) Inverse variance weighted   11 0.04168546 0.00807618 2.449369e-07
# 4 PCOS (Venkatesh)        BMI     BMI PCOS (Venkatesh)               Simple mode   11 0.05892587 0.02201603 2.323193e-02
# 5 PCOS (Venkatesh)        BMI     BMI PCOS (Venkatesh)             Weighted mode   11 0.01609010 0.01945153 4.274275e-01

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method     Estimate          SE       CI_low      CI_upp         P
# 1  Egger fixed effects 0.0003168633 0.002458682 -0.004502064 0.005135791 0.8974564
# 2 Egger random effects 0.0003168633 0.002565703 -0.004711822 0.005345549 0.4508558
# 
# [[1]]$Q
# Method          Q df         P
# 1   Q_ivw 9.81716241 10 0.4566782
# 2 Q_egger 9.80055358  9 0.3668717
# 3  Q_diff 0.01660883  1 0.8974564
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/venkatesh/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/venkatesh/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/venkatesh/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//3_reverse_mr/2_mr/whradjbmi/venkatesh/strict_het_test_after_outlier_extraction.RDS")


