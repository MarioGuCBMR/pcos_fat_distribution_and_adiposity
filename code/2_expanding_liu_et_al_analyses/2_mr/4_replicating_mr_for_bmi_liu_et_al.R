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

path_2_input <- "your_path"

setwd(path_2_input)

#Let's load the data for all WHRadjBMI mathches

liu_pcos_doctor <- fread("output/2_replicating_liu_et_al/2_mr/1_expo_outcome_df/bmi_liu_PCOS.txt")
liu_pcos_broad <- fread("output/2_replicating_liu_et_al/2_mr/1_expo_outcome_df/bmi_liu_PCOS (Broad).txt")
liu_pcos_consortium <- fread("output/2_replicating_liu_et_al/2_mr/1_expo_outcome_df/bmi_liu_PCOS (Consortium).txt")
liu_pcos_meta_analysis <- fread("output/2_replicating_liu_et_al/2_mr/1_expo_outcome_df/bmi_liu_PCOS (Day et al).txt")
liu_pcos_adj_age <- fread("output/2_replicating_liu_et_al/2_mr/1_expo_outcome_df/bmi_liu_PCOS (adj age).txt")
liu_pcos_adj_age_bmi <- fread("output/2_replicating_liu_et_al/2_mr/1_expo_outcome_df/bmi_liu_PCOS (adj age+BMI).txt")

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

mr_df$id.exposure <- "BMI"
mr_df$id.outcome <- "PCOS (meta-analysis)"
mr_df$exposure <- "BMI"
mr_df$outcome <- "PCOS (meta-analysis)"
mr_df$units.exposure <- "SD"
mr_df$units.outcome <- "LOR"

mr_steiger_df <- TwoSampleMR::steiger_filtering(mr_df)
mr_steiger_df <- mr_steiger_df[which(mr_steiger_df$steiger_dir == TRUE),] #259
mr_steiger_df$mr_keep <- TRUE

#Let's save the data:

dir.create("output/2_replicating_liu_et_al/2_mr/2_mr/")
dir.create("output/2_replicating_liu_et_al/2_mr/2_mr/bmi/")
dir.create("output/2_replicating_liu_et_al/2_mr/2_mr/bmi/meta_analysis/")
saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al/2_mr/2_mr/bmi/meta_analysis/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

#id.exposure           id.outcome              outcome exposure                    method nsnp         b        se         pval
#1         BMI PCOS (meta-analysis) PCOS (meta-analysis)      BMI                  MR Egger  259 1.7163951 0.3063766 5.436103e-08
#2         BMI PCOS (meta-analysis) PCOS (meta-analysis)      BMI           Weighted median  259 1.0258410 0.1751019 4.669668e-09
#3         BMI PCOS (meta-analysis) PCOS (meta-analysis)      BMI Inverse variance weighted  259 0.9066126 0.1083403 5.850745e-17
#4         BMI PCOS (meta-analysis) PCOS (meta-analysis)      BMI               Simple mode  259 1.9182835 0.5868733 1.227439e-03
#5         BMI PCOS (meta-analysis) PCOS (meta-analysis)      BMI             Weighted mode  259 1.7971054 0.3705776 2.142012e-06

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method    Estimate          SE      CI_low       CI_upp           P
# 1  Egger fixed effects -0.01753452 0.006169648 -0.02962681 -0.005442237 0.004482272
# 2 Egger random effects -0.01753452 0.006169648 -0.02962681 -0.005442237 0.997758864
# 
# [[1]]$Q
# Method          Q  df           P
# 1   Q_ivw 261.523109 258 0.427117966
# 2 Q_egger 253.554092 257 0.549038816
# 3  Q_diff   7.969016   1 0.004758479

strict_het <- compute_strict_het(mr_steiger_df)

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al/2_mr/2_mr/bmi/meta_analysis/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al/2_mr/2_mr/bmi/meta_analysis/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al/2_mr/2_mr/bmi/meta_analysis//sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al/2_mr/2_mr/bmi/meta_analysis//strict_het_test_before_outlier_extraction.RDS")

#No heterogeneity, no pleiotropy... so no need for any tough stuff

mr_post <- remove_outlier(mr_steiger_df, rucker)

saveRDS(mr_post, "output/2_replicating_liu_et_al/2_mr/2_mr/bmi/meta_analysis/instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

#id.exposure           id.outcome              outcome exposure                    method nsnp         b        se         pval
#1         BMI PCOS (meta-analysis) PCOS (meta-analysis)      BMI                  MR Egger  250 1.1690588 0.3562058 1.179077e-03
#2         BMI PCOS (meta-analysis) PCOS (meta-analysis)      BMI           Weighted median  250 0.9942966 0.1853537 8.125514e-08
#3         BMI PCOS (meta-analysis) PCOS (meta-analysis)      BMI Inverse variance weighted  250 0.8285690 0.1130266 2.288804e-13
#4         BMI PCOS (meta-analysis) PCOS (meta-analysis)      BMI               Simple mode  250 1.9090451 0.5654594 8.528345e-04
#5         BMI PCOS (meta-analysis) PCOS (meta-analysis)      BMI             Weighted mode  250 1.4743710 0.3961635 2.447376e-04

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method     Estimate          SE      CI_low      CI_upp        P
# 1  Egger fixed effects -0.007062322 0.006488549 -0.01977964 0.005655001 0.276406
# 2 Egger random effects -0.007062322 0.006488549 -0.01977964 0.005655001 0.861797
# 
# [[1]]$Q
# Method        Q  df         P
# 1   Q_ivw 213.7053 249 0.9487585
# 2 Q_egger 212.6893 248 0.9492242
# 3  Q_diff   1.0160   1 0.3134697
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/bmi/meta_analysis/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/bmi/meta_analysis/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/bmi/meta_analysis/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/bmi/meta_analysis/strict_het_test_after_outlier_extraction.RDS")

###############################################################
#Let's first analyse the data for the whole set for bmi - PCOS#
###############################################################

#STEP 1: let's use only those that have F>10

mr_df <- liu_pcos_doctor[which(((liu_pcos_doctor$beta.exposure^2)/(liu_pcos_doctor$se.exposure^2)) > 10),] #250

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

dir.create("output/2_replicating_liu_et_al//2_mr/2_mr/bmi/")
dir.create("output/2_replicating_liu_et_al//2_mr/2_mr/bmi/doctor")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//2_mr/2_mr/bmi/doctor/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp         b        se         pval
# 1    exposure    outcome outcome exposure                  MR Egger  275 0.7334426 0.3344885 2.917083e-02
# 2    exposure    outcome outcome exposure           Weighted median  275 0.5790306 0.1867581 1.932388e-03
# 3    exposure    outcome outcome exposure Inverse variance weighted  275 0.6726508 0.1156346 5.989907e-09
# 4    exposure    outcome outcome exposure               Simple mode  275 0.6115877 0.4873679 2.105921e-01
# 5    exposure    outcome outcome exposure             Weighted mode  275 0.6115877 0.3186109 5.595387e-02

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method     Estimate          SE      CI_low     CI_upp         P
# 1  Egger fixed effects -0.001299226 0.005999493 -0.01305802 0.01045956 0.8285544
# 2 Egger random effects -0.001299226 0.006706231 -0.01444320 0.01184474 0.5768080
# 
# [[1]]$Q
# Method           Q  df           P
# 1   Q_ivw 341.1537711 274 0.003562577
# 2 Q_egger 341.1068746 273 0.003162593
# 3  Q_diff   0.0468965   1 0.828554384
# 
# [[1]]$res
# [1] "B"

strict_het <- compute_strict_het(mr_steiger_df)

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/bmi/doctor/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/bmi/doctor/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/bmi/doctor/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/bmi/doctor/strict_het_test_before_outlier_extraction.RDS")

#Heterogeneity comes from variants with high SEs, which is controlled by the methods. All agree on the directions and have very similar effects.
#No need to remove outliers.

###############################################################
#Let's first analyse the data for the whole set for bmi - PCOS# broad definition - anovulation
###############################################################

#STEP 1: let's use only those that have F>10

mr_df <- liu_pcos_broad[which(((liu_pcos_broad$beta.exposure^2)/(liu_pcos_broad$se.exposure^2)) > 10),] #250

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

dir.create("output/2_replicating_liu_et_al//2_mr/2_mr/bmi/broad")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//2_mr/2_mr/bmi/broad/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp         b         se         pval
# 1    exposure    outcome outcome exposure                  MR Egger  275 0.6236240 0.21300606 3.702226e-03
# 2    exposure    outcome outcome exposure           Weighted median  275 0.3937326 0.12309957 1.381497e-03
# 3    exposure    outcome outcome exposure Inverse variance weighted  275 0.3551419 0.07385828 1.521250e-06
# 4    exposure    outcome outcome exposure               Simple mode  275 0.5429965 0.32753404 9.849573e-02
# 5    exposure    outcome outcome exposure             Weighted mode  275 0.5214703 0.21726128 1.705358e-02

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method    Estimate          SE      CI_low      CI_upp         P
# 1  Egger fixed effects -0.00573712 0.004098573 -0.01377018 0.002295937 0.1615779
# 2 Egger random effects -0.00573712 0.004270138 -0.01410644 0.002632196 0.9104521
# 
# [[1]]$Q
# Method          Q  df         P
# 1   Q_ivw 298.293014 274 0.1499028
# 2 Q_egger 296.333617 273 0.1587003
# 3  Q_diff   1.959397   1 0.1615779
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/bmi/broad/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/bmi/broad/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/bmi/broad/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/bmi/broad/strict_het_test_before_outlier_extraction.RDS")

#Very slight heterogeneity - Funnel plot looks OK. No need to remove outliers.

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

dir.create("output/2_replicating_liu_et_al//2_mr/2_mr/bmi/consortium")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//2_mr/2_mr/bmi/consortium/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp          b         se         pval
# 1    exposure    outcome outcome exposure                  MR Egger  276 0.06046123 0.08669526 4.861442e-01
# 2    exposure    outcome outcome exposure           Weighted median  276 0.12367840 0.04523388 6.253230e-03
# 3    exposure    outcome outcome exposure Inverse variance weighted  276 0.19693362 0.03010932 6.126407e-11
# 4    exposure    outcome outcome exposure               Simple mode  276 0.02243765 0.12882363 8.618574e-01
# 5    exposure    outcome outcome exposure             Weighted mode  276 0.07920336 0.08005254 3.233410e-01

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method    Estimate          SE        CI_low      CI_upp          P
# 1  Egger fixed effects 0.002914392 0.001477659  1.823389e-05 0.005810551 0.04857495
# 2 Egger random effects 0.002914392 0.001736931 -4.899303e-04 0.006318715 0.04668357
# 
# [[1]]$Q
# Method          Q  df            P
# 1   Q_ivw 382.478422 275 1.886452e-05
# 2 Q_egger 378.588441 274 2.819125e-05
# 3  Q_diff   3.889982   1 4.857495e-02
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/bmi/consortium/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/bmi/consortium/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/bmi/consortium/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/bmi/consortium/strict_het_test_before_outlier_extraction.RDS")

#Very small heterogeneity, overall. However, the interecept is the intercept

#STEP 5: removing outliers!

mr_post <- remove_outlier(mr_steiger_df, rucker)
  
saveRDS(mr_post, "output/2_replicating_liu_et_al//2_mr/2_mr/bmi/consortium/instruments_after_outlier_removal_df.RDS")
  
#STEP 6: let's run the results for real now:
  
mr_res <- TwoSampleMR::mr(mr_post)
  
# id.exposure id.outcome outcome exposure                    method nsnp           b         se         pval
# 1    exposure    outcome outcome exposure                  MR Egger  247 0.056831770 0.07730106 4.629206e-01
# 2    exposure    outcome outcome exposure           Weighted median  247 0.119593073 0.04557912 8.694042e-03
# 3    exposure    outcome outcome exposure Inverse variance weighted  247 0.183151678 0.02704946 1.279014e-11
# 4    exposure    outcome outcome exposure               Simple mode  247 0.004427853 0.11991397 9.705746e-01
# 5    exposure    outcome outcome exposure             Weighted mode  247 0.077706189 0.07475054 2.995743e-01

# #Let's check the plot:
#  
rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method    Estimate          SE        CI_low      CI_upp          P
# 1  Egger fixed effects 0.002687618 0.001434164 -0.0001232921 0.005498529 0.06093096
# 2 Egger random effects 0.002687618 0.001434164 -0.0001232921 0.005498529 0.03046548
# 
# [[1]]$Q
# Method         Q  df          P
# 1   Q_ivw 215.33211 246 0.92135761
# 2 Q_egger 212.28913 245 0.93557736
# 3  Q_diff   3.04298   1 0.08108696
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/bmi/consortium/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/bmi/consortium/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/bmi/consortium/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/bmi/consortium/strict_het_test_after_outlier_extraction.RDS")

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

dir.create("output/2_replicating_liu_et_al//2_mr/2_mr/bmi/age_adjusted")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//2_mr/2_mr/bmi/age_adjusted/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp          b         se         pval
# 1    exposure    outcome outcome exposure                  MR Egger  275  0.7408916 0.24858905 0.0031382375
# 2    exposure    outcome outcome exposure           Weighted median  275  0.4938321 0.14577068 0.0007047288
# 3    exposure    outcome outcome exposure Inverse variance weighted  275  0.3256303 0.08687274 0.0001779949
# 4    exposure    outcome outcome exposure               Simple mode  275 -0.1474837 0.44801663 0.7422616229
# 5    exposure    outcome outcome exposure             Weighted mode  275  0.6081913 0.28151921 0.0316094490

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method     Estimate          SE      CI_low       CI_upp          P
# 1  Egger fixed effects -0.008899382 0.004891475 -0.01848650 0.0006877338 0.06885568
# 2 Egger random effects -0.008899382 0.004994300 -0.01868803 0.0008892666 0.96261786
# 
# [[1]]$Q
# Method          Q  df          P
# 1   Q_ivw 287.908310 274 0.27002632
# 2 Q_egger 284.598219 273 0.30225212
# 3  Q_diff   3.310091   1 0.06885568
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/bmi/age_adjusted/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/bmi/age_adjusted/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/bmi/age_adjusted/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/bmi/age_adjusted/strict_het_test_before_outlier_extraction.RDS")

#There is a bit of asymmetry. I am almost certain that Wmed and Wmod correct for this, but we are gonna run the analyses again just to be on the safe side.

mr_post <- remove_outlier(mr_steiger_df, rucker)

saveRDS(mr_post, "output/2_replicating_liu_et_al//2_mr/2_mr/bmi/age_adjusted/instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure                    method nsnp          b         se         pval
# 1    exposure    outcome outcome exposure                  MR Egger  259  0.7969216 0.25138553 1.708474e-03
# 2    exposure    outcome outcome exposure           Weighted median  259  0.5673460 0.15306923 2.101707e-04
# 3    exposure    outcome outcome exposure Inverse variance weighted  259  0.4244093 0.08774916 1.320660e-06
# 4    exposure    outcome outcome exposure               Simple mode  259 -0.2429821 0.43697517 5.786554e-01
# 5    exposure    outcome outcome exposure             Weighted mode  259  0.6212669 0.28328387 2.919379e-02

# #Let's check the plot:
#  
rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method     Estimate          SE      CI_low       CI_upp          P
# 1  Egger fixed effects -0.007952474 0.004355944 -0.01648997 0.0005850196 0.06790146
# 2 Egger random effects -0.007952474 0.004355944 -0.01648997 0.0005850196 0.96604927
# 
# [[1]]$Q
# Method          Q  df         P
# 1   Q_ivw 195.307331 258 0.9986129
# 2 Q_egger 192.806818 257 0.9989610
# 3  Q_diff   2.500513   1 0.1138092
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/bmi/age_adjusted/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/bmi/age_adjusted/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/bmi/age_adjusted/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/bmi/age_adjusted/strict_het_test_after_outlier_extraction.RDS")

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

dir.create("output/2_replicating_liu_et_al//2_mr/2_mr/bmi/age_and_bmi_adjusted")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//2_mr/2_mr/bmi/age_and_bmi_adjusted/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp          b         se      pval
# 1    exposure    outcome outcome exposure                  MR Egger  275 0.27690169 0.28543683 0.3328563
# 2    exposure    outcome outcome exposure           Weighted median  275 0.14369737 0.16964581 0.3969709
# 3    exposure    outcome outcome exposure Inverse variance weighted  275 0.07926089 0.09944844 0.4254482
# 4    exposure    outcome outcome exposure               Simple mode  275 0.43823771 0.48117838 0.3632227
# 5    exposure    outcome outcome exposure             Weighted mode  275 0.30370600 0.31575495 0.3369776

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method    Estimate          SE     CI_low      CI_upp         P
# 1  Egger fixed effects -0.00423703 0.005558708 -0.0151319 0.006657837 0.4459210
# 2 Egger random effects -0.00423703 0.005558708 -0.0151319 0.006657837 0.7770395
# 
# [[1]]$Q
# Method           Q  df         P
# 1   Q_ivw 256.9487876 274 0.7628315
# 2 Q_egger 256.4031101 273 0.7569500
# 3  Q_diff   0.5456775   1 0.4600893
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/bmi/age_and_bmi_adjusted/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/bmi/age_and_bmi_adjusted/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/bmi/age_and_bmi_adjusted/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/bmi/age_and_bmi_adjusted/strict_het_test_before_outlier_extraction.RDS")


