##############
#INTRODUCTION#
##############

#This code performs MR analyses between vat_asat_ratio and PCOS - replicating the results from Liu et al.

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

liu_pcos_doctor <- fread("output/2_replicating_liu_et_al/4_mr_fat_depots/1_expo_outcome_df/vat_asat_ratio_PCOS.txt")
liu_pcos_broad <- fread("output/2_replicating_liu_et_al/4_mr_fat_depots/1_expo_outcome_df/vat_asat_ratio_PCOS (Broad).txt")
liu_pcos_consortium <- fread("output/2_replicating_liu_et_al/4_mr_fat_depots/1_expo_outcome_df/vat_asat_ratio_PCOS (Consortium).txt")
liu_pcos_meta_analysis <- fread("output/2_replicating_liu_et_al/4_mr_fat_depots/1_expo_outcome_df/vat_asat_ratio_PCOS (Day et al).txt")
liu_pcos_venkatesh <- fread("output/2_replicating_liu_et_al/4_mr_fat_depots/1_expo_outcome_df/vat_asat_ratio_PCOS (Venkatesh et al).txt")
liu_pcos_adj_age <- fread("output/2_replicating_liu_et_al/4_mr_fat_depots/1_expo_outcome_df/vat_asat_ratio_PCOS (adj age).txt")
liu_pcos_adj_age_bmi <- fread("output/2_replicating_liu_et_al/4_mr_fat_depots/1_expo_outcome_df/vat_asat_ratio_PCOS (adj age+BMI).txt")

##############################################################
#STEP 2: let's first rerun the analyses for the meta-analysis#
##############################################################

mr_df <- liu_pcos_meta_analysis[which(((as.numeric(liu_pcos_meta_analysis$beta.exposure)^2)/(as.numeric(liu_pcos_meta_analysis$se.exposure)^2)) > 10),] #9

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

mr_df$id.exposure <- "vat_asat_ratio"
mr_df$id.outcome <- "PCOS (meta-analysis)"
mr_df$exposure <- "vat_asat_ratio"
mr_df$outcome <- "PCOS (meta-analysis)"
mr_df$units.exposure <- "SD"
mr_df$units.outcome <- "LOR"

mr_steiger_df <- TwoSampleMR::steiger_filtering(mr_df)
mr_steiger_df <- mr_steiger_df[which(mr_steiger_df$steiger_dir == TRUE),] #19
mr_steiger_df$mr_keep <- TRUE

#Let's save the data:

dir.create("output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/")
dir.create("output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/vat_asat_ratio/")
dir.create("output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/vat_asat_ratio/meta_analysis/")
saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/vat_asat_ratio/meta_analysis/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure           id.outcome              outcome       exposure                    method nsnp           b        se      pval
# 1 vat_asat_ratio PCOS (meta-analysis) PCOS (meta-analysis) vat_asat_ratio                  MR Egger    8 -0.73622297 0.7594652 0.3697851
# 2 vat_asat_ratio PCOS (meta-analysis) PCOS (meta-analysis) vat_asat_ratio           Weighted median    8 -0.06293461 0.2382330 0.7916470
# 3 vat_asat_ratio PCOS (meta-analysis) PCOS (meta-analysis) vat_asat_ratio Inverse variance weighted    8  0.12866190 0.1862547 0.4897007
# 4 vat_asat_ratio PCOS (meta-analysis) PCOS (meta-analysis) vat_asat_ratio               Simple mode    8 -0.13361333 0.3701558 0.7287734
# 5 vat_asat_ratio PCOS (meta-analysis) PCOS (meta-analysis) vat_asat_ratio             Weighted mode    8 -0.11251383 0.3562142 0.7613181

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method   Estimate         SE      CI_low    CI_upp         P
# 1  Egger fixed effects 0.06587575 0.05301224 -0.03802634 0.1697778 0.2139962
# 2 Egger random effects 0.06587575 0.05617095 -0.04421728 0.1759688 0.1204435
# 
# [[1]]$Q
# Method        Q df         P
# 1   Q_ivw 8.280498  7 0.3085110
# 2 Q_egger 6.736315  6 0.3459225
# 3  Q_diff 1.544183  1 0.2139962
# 
# [[1]]$res
# [1] "A"

strict_het <- compute_strict_het(mr_steiger_df)

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 < 0.0001 [0.0000; 1.0678]; tau = 0.0004 [0.0000; 1.0333]
# I^2 = 11.2% [0.0%; 71.2%]; H = 1.06 [1.00; 1.86]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 7.88    7  0.3429

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/vat_asat_ratio/meta_analysis/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/vat_asat_ratio/meta_analysis/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/vat_asat_ratio/meta_analysis//sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/vat_asat_ratio/meta_analysis//strict_het_test_before_outlier_extraction.RDS")

#Slight asymmetry, but very tame 

mr_post <- remove_outlier(mr_steiger_df, rucker) #no significant outliers!

#saveRDS(mr_post, "output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/vat_asat_ratio/meta_analysis/instruments_after_outlier_removal_df.RDS")

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
# ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/meta_analysis/plots_after_outlier_extraction.svg", height=16, width=16, units="in")
# 
# saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/meta_analysis/mr_res_after_outlier_extraction.RDS")
# saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/meta_analysis/sensitivity_test_after_outlier_extraction.RDS")
# saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/meta_analysis/strict_het_test_after_outlier_extraction.RDS")

###############################################################
#Let's first analyse the data for the whole set for vat_asat_ratio - PCOS#
###############################################################

#STEP 1: let's use only those that have F>10

mr_df <- liu_pcos_doctor[which(((liu_pcos_doctor$beta.exposure^2)/(liu_pcos_doctor$se.exposure^2)) > 10),] #8/8

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

dir.create("output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/")
dir.create("output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/doctor")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/doctor/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp           b        se      pval
# 1    exposure    outcome outcome exposure                  MR Egger    9  0.47774776 0.7532884 0.5460990
# 2    exposure    outcome outcome exposure           Weighted median    9  0.08613317 0.2282018 0.7058444
# 3    exposure    outcome outcome exposure Inverse variance weighted    9  0.23281484 0.2320559 0.3157304
# 4    exposure    outcome outcome exposure               Simple mode    9  0.68918784 0.5295061 0.2292924
# 5    exposure    outcome outcome exposure             Weighted mode    9 -0.31282967 0.4157515 0.4733404

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method   Estimate         SE      CI_low     CI_upp         P
# 1  Egger fixed effects -0.0213431 0.03645099 -0.09278573 0.05009953 0.5581923
# 2 Egger random effects -0.0213431 0.06204141 -0.14294204 0.10025584 0.6345820
# 
# [[1]]$Q
# Method          Q df           P
# 1   Q_ivw 20.6216683  8 0.008223312
# 2 Q_egger 20.2788244  7 0.004997890
# 3  Q_diff  0.3428439  1 0.558192270
# 
# [[1]]$res
# [1] "B"

strict_het <- compute_strict_het(mr_steiger_df)

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0.2706 [0.0147; 1.6373]; tau = 0.5202 [0.1210; 1.2796]
# I^2 = 57.9% [11.9%; 79.9%]; H = 1.54 [1.07; 2.23]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 19.02    8  0.0148

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/doctor/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/doctor/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/doctor/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/doctor/strict_het_test_before_outlier_extraction.RDS")

#Heterogeneity, though funnel plots looks decent.
#RadialMR will most likely detect outliers, but might produce issues

mr_post <- remove_outlier(mr_steiger_df, rucker) #outliers detected!

saveRDS(mr_post, "output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/vat_asat_ratio/doctor/instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure                    method nsnp           b        se      pval
# 1    exposure    outcome outcome exposure                  MR Egger    6  1.11385875 1.1714874 0.3955433
# 2    exposure    outcome outcome exposure           Weighted median    6  0.03371131 0.2670978 0.8995630
# 3    exposure    outcome outcome exposure Inverse variance weighted    6  0.12045558 0.2178777 0.5803601
# 4    exposure    outcome outcome exposure               Simple mode    6 -0.19925361 0.4220904 0.6567650
# 5    exposure    outcome outcome exposure             Weighted mode    6 -0.19925361 0.3817289 0.6239770

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_post)
 
# [[1]]$intercept
# Method   Estimate         SE     CI_low     CI_upp         P
# 1  Egger fixed effects -0.0735398 0.07351511 -0.2176268 0.07054718 0.3171481
# 2 Egger random effects -0.0735398 0.08512832 -0.2403882 0.09330865 0.8061703
# 
# [[1]]$Q
# Method        Q df         P
# 1   Q_ivw 6.364253  5 0.2723710
# 2 Q_egger 5.363581  4 0.2519836
# 3  Q_diff 1.000672  1 0.3171481
# 
# [[1]]$res
# [1] "A"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_post) #they disagree, but still

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0.0591 [0.0000; 1.6092]; tau = 0.2431 [0.0000; 1.2685]
# I^2 = 18.8% [0.0%; 63.6%]; H = 1.11 [1.00; 1.66]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 6.16    5  0.2911

#That is quite good:

mr_post$labels <- NA
mr_post$units.exposure <- "SD"
mr_post$units.outcome <- "LOR"
 
plotio=mr_plots(mr_post)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/doctor/plots_after_outlier_extraction.svg", height=16, width=16, units="in")
 
saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/doctor/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/doctor/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/doctor/strict_het_test_after_outlier_extraction.RDS")

#Heterogeneity removed and data holds.

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

dir.create("output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/broad")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/broad/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp           b        se      pval
# 1    exposure    outcome outcome exposure                  MR Egger    9  0.20169495 0.6741557 0.7734889
# 2    exposure    outcome outcome exposure           Weighted median    9 -0.14642438 0.1773546 0.4090297
# 3    exposure    outcome outcome exposure Inverse variance weighted    9  0.07374892 0.2053543 0.7194977
# 4    exposure    outcome outcome exposure               Simple mode    9 -0.27033776 0.4445875 0.5600073
# 5    exposure    outcome outcome exposure             Weighted mode    9 -0.37418319 0.2401355 0.1577984

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method    Estimate         SE      CI_low     CI_upp         P
# 1  Egger fixed effects -0.01111483 0.02504704 -0.06020613 0.03797647 0.6572174
# 2 Egger random effects -0.01111483 0.05539136 -0.11967990 0.09745024 0.5795178
# 
# [[1]]$Q
# Method          Q df            P
# 1   Q_ivw 34.4318311  8 3.392156e-05
# 2 Q_egger 34.2349097  7 1.556569e-05
# 3  Q_diff  0.1969214  1 6.572174e-01
# 
# [[1]]$res
# [1] "B"

#Here there is heterogeneity...

strict_het <- compute_strict_het(mr_steiger_df) 

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0.2709 [0.0671; 1.3489]; tau = 0.5205 [0.2590; 1.1614]
# I^2 = 73.8% [48.9%; 86.5%]; H = 1.95 [1.40; 2.72]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 30.49    8  0.0002

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/broad/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/broad/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/broad/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/broad/strict_het_test_before_outlier_extraction.RDS")

#Asymmetry and heterogeneity introduced!

mr_post <- remove_outlier(mr_steiger_df, rucker)
 
saveRDS(mr_post, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/broad/instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure                    method nsnp           b        se      pval
# 1    exposure    outcome outcome exposure                  MR Egger    5  0.99399800 0.9066164 0.3530362
# 2    exposure    outcome outcome exposure           Weighted median    5  0.05061415 0.1963130 0.7965430
# 3    exposure    outcome outcome exposure Inverse variance weighted    5  0.12251232 0.1892798 0.5174668
# 4    exposure    outcome outcome exposure               Simple mode    5 -0.12961360 0.3811118 0.7508953
# 5    exposure    outcome outcome exposure             Weighted mode    5 -0.18420450 0.3674493 0.6424886

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method    Estimate         SE     CI_low     CI_upp         P
# 1  Egger fixed effects -0.06348227 0.05075162 -0.1629536 0.03598909 0.2109921
# 2 Egger random effects -0.06348227 0.06457357 -0.1900441 0.06307959 0.8372209
# 
# [[1]]$Q
# Method        Q df         P
# 1   Q_ivw 6.421190  4 0.1698243
# 2 Q_egger 4.856584  3 0.1826051
# 3  Q_diff 1.564606  1 0.2109921
# 
# [[1]]$res
# [1] "A"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_post) #they disagree, but still

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0.0700 [0.0000; 1.3291]; tau = 0.2646 [0.0000; 1.1529]
# I^2 = 35.3% [0.0%; 75.7%]; H = 1.24 [1.00; 2.03]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 6.18    4  0.1858

#That is quite good:

mr_post$labels <- NA
mr_post$units.exposure <- "SD"
mr_post$units.outcome <- "LOR"
 
plotio=mr_plots(mr_post)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/broad/plots_after_outlier_extraction.svg", height=16, width=16, units="in")
 
saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/broad/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/broad/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/broad/strict_het_test_after_outlier_extraction.RDS")

#Heterogeneity controled -  results hold

################################################################
#Let's first analyse the data for the whole set for vat_asat_ratio - PCOS# consortium definition
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

dir.create("output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/consortium")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/consortium/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp          b         se      pval
# 1    exposure    outcome outcome exposure                  MR Egger    9 0.14927658 0.26031637 0.5842890
# 2    exposure    outcome outcome exposure           Weighted median    9 0.06071651 0.06035685 0.3144354
# 3    exposure    outcome outcome exposure Inverse variance weighted    9 0.03053303 0.07932260 0.7002950
# 4    exposure    outcome outcome exposure               Simple mode    9 0.08444226 0.09802609 0.4140758
# 5    exposure    outcome outcome exposure             Weighted mode    9 0.09690349 0.09515633 0.3383119

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method    Estimate          SE      CI_low      CI_upp         P
# 1  Egger fixed effects -0.01024685 0.009172788 -0.02822519 0.007731482 0.2639548
# 2 Egger random effects -0.01024685 0.021278812 -0.05195256 0.031458853 0.6849378
# 
# [[1]]$Q
# Method         Q df            P
# 1   Q_ivw 38.917415  8 5.091859e-06
# 2 Q_egger 37.669520  7 3.501527e-06
# 3  Q_diff  1.247895  1 2.639548e-01
# 
# [[1]]$res
# [1] "B"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_steiger_df) #they disagree, but still

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0.0463 [0.0135; 0.2601]; tau = 0.2153 [0.1160; 0.5100]
# I^2 = 75.7% [53.3%; 87.4%]; H = 2.03 [1.46; 2.81]
# 
# Test of heterogeneity:
#   Q d.f.  p-value
# 32.91    8 < 0.0001

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/consortium/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/consortium/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/consortium/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/consortium/strict_het_test_before_outlier_extraction.RDS")

#We have heterogeneity and asymmetry

#STEP 5: removing outliers!

mr_post <- remove_outlier(mr_steiger_df, rucker) #but not significant outliers!!
  
saveRDS(mr_post, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/consortium/instruments_after_outlier_removal_df.RDS")
  
#STEP 6: let's run the results for real now:
  
mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure                    method nsnp          b         se       pval
# 1    exposure    outcome outcome exposure                  MR Egger    5 0.24789928 0.14684313 0.18995804
# 2    exposure    outcome outcome exposure           Weighted median    5 0.09378872 0.05693318 0.09948668
# 3    exposure    outcome outcome exposure Inverse variance weighted    5 0.08104121 0.04598956 0.07804140
# 4    exposure    outcome outcome exposure               Simple mode    5 0.09051770 0.08713645 0.35756289
# 5    exposure    outcome outcome exposure             Weighted mode    5 0.12229097 0.08751641 0.23483545

# #Let's check the plot:
  
rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method    Estimate          SE      CI_low      CI_upp         P
# 1  Egger fixed effects -0.01457682 0.009533094 -0.03326134 0.004107701 0.1262458
# 2 Egger random effects -0.01457682 0.009533094 -0.03326134 0.004107701 0.9368771
# 
# [[1]]$Q
# Method        Q df         P
# 1   Q_ivw 3.268505  4 0.5139378
# 2 Q_egger 1.836902  3 0.6069377
# 3  Q_diff 1.431603  1 0.2315030
# 
# [[1]]$res
# [1] "A"

# #Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_post) #they disagree, but still

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0 [0.0000; 0.0533]; tau = 0 [0.0000; 0.2308]
# I^2 = 0.0% [0.0%; 79.2%]; H = 1.00 [1.00; 2.19]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 3.18    4  0.5274

#That is quite good:
  
mr_post$labels <- NA
mr_post$units.exposure <- "SD"
mr_post$units.outcome <- "LOR"
  
plotio=mr_plots(mr_post)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/consortium/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/consortium/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/consortium/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/consortium/strict_het_test_after_outlier_extraction.RDS")

#I think we might have now more asymmetry, despite heterogeneity not popping up...

#mr_post_2 <- remove_outlier(mr_post, rucker) #cannot detect it...

###############################################################
#Let's first analyse the data for the whole set for vat_asat_ratio - PCOS# age adjusted
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

dir.create("output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/age_adjusted")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/age_adjusted/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp           b        se      pval
# 1    exposure    outcome outcome exposure                  MR Egger    8 -0.11422846 0.6212758 0.8601788
# 2    exposure    outcome outcome exposure           Weighted median    8  0.09711306 0.1936258 0.6159839
# 3    exposure    outcome outcome exposure Inverse variance weighted    8  0.13591414 0.1757727 0.4393814
# 4    exposure    outcome outcome exposure               Simple mode    8  0.15371715 0.3120974 0.6374163
# 5    exposure    outcome outcome exposure             Weighted mode    8  0.12108696 0.3084817 0.7063481

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method   Estimate         SE      CI_low     CI_upp         P
# 1  Egger fixed effects 0.02082227 0.03233176 -0.04254682 0.08419137 0.5195631
# 2 Egger random effects 0.02082227 0.04931518 -0.07583371 0.11747825 0.3364291
# 
# [[1]]$Q
# Method          Q df          P
# 1   Q_ivw 14.3737418  7 0.04491953
# 2 Q_egger 13.9589812  6 0.03009774
# 3  Q_diff  0.4147605  1 0.51956308
# 
# [[1]]$res
# [1] "B"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_steiger_df) #they disagree, but still
 
# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0.1106 [0.0000; 1.0543]; tau = 0.3325 [0.0000; 1.0268]
# I^2 = 48.0% [0.0%; 76.9%]; H = 1.39 [1.00; 2.08]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 13.47    7  0.0614

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/age_adjusted//plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/age_adjusted//mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/age_adjusted//sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/age_adjusted//strict_het_test_before_outlier_extraction.RDS")

#Very slight asymmetry.

mr_post <- remove_outlier(mr_steiger_df, rucker) #no significant outliers. 

saveRDS(mr_post, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/age_adjusted//instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure                    method nsnp         b        se       pval
# 1    exposure    outcome outcome exposure                  MR Egger    7 0.4137828 0.6026356 0.52286083
# 2    exposure    outcome outcome exposure           Weighted median    7 0.2270249 0.1857192 0.22155277
# 3    exposure    outcome outcome exposure Inverse variance weighted    7 0.2866277 0.1635048 0.07959805
# 4    exposure    outcome outcome exposure               Simple mode    7 0.1822547 0.2776133 0.53585045
# 5    exposure    outcome outcome exposure             Weighted mode    7 0.1936002 0.2393524 0.44947670

# #Let's check the plot:
  
rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method     Estimate         SE      CI_low     CI_upp         P
# 1  Egger fixed effects -0.009986114 0.03487385 -0.07833760 0.05836537 0.7746103
# 2 Egger random effects -0.009986114 0.04521037 -0.09859682 0.07862459 0.5874075
# 
# [[1]]$Q
# Method          Q df         P
# 1   Q_ivw 8.48523031  6 0.2046643
# 2 Q_egger 8.40323419  5 0.1353683
# 3  Q_diff 0.08199612  1 0.7746103
# 
# [[1]]$res
# [1] "A"

# #Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_post) 

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0.0071 [0.0000; 1.0470]; tau = 0.0840 [0.0000; 1.0232]
# I^2 = 24.5% [0.0%; 66.8%]; H = 1.15 [1.00; 1.74]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 7.95    6  0.2421

#That is quite good:

mr_post$labels <- NA
mr_post$units.exposure <- "SD"
mr_post$units.outcome <- "LOR"

plotio=mr_plots(mr_post)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/age_adjusted/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/age_adjusted/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/age_adjusted/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/age_adjusted/strict_het_test_after_outlier_extraction.RDS")

#Heterogeneity removed - results hold

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

dir.create("output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/age_and_bmi_adjusted")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/age_and_bmi_adjusted/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)
 
# id.exposure id.outcome outcome exposure                    method nsnp           b        se      pval
# 1    exposure    outcome outcome exposure                  MR Egger    8 -0.44848163 0.6250223 0.5000076
# 2    exposure    outcome outcome exposure           Weighted median    8  0.12662573 0.2001997 0.5270620
# 3    exposure    outcome outcome exposure Inverse variance weighted    8  0.09876669 0.1859434 0.5953042
# 4    exposure    outcome outcome exposure               Simple mode    8 -0.02891618 0.2854756 0.9221596
# 5    exposure    outcome outcome exposure             Weighted mode    8 -0.04459100 0.2773219 0.8767996

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method   Estimate         SE      CI_low    CI_upp         P
# 1  Egger fixed effects 0.04553534 0.03799192 -0.02892746 0.1199981 0.2307017
# 2 Egger random effects 0.04553534 0.04959663 -0.05167227 0.1427429 0.1792797
# 
# [[1]]$Q
# Method        Q df         P
# 1   Q_ivw 11.66176  7 0.1122452
# 2 Q_egger 10.22523  6 0.1154819
# 3  Q_diff  1.43653  1 0.2307017
# 
# [[1]]$res
# [1] "A"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_steiger_df) #they disagree, but still

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0.0543 [0.0000; 1.2930]; tau = 0.2331 [0.0000; 1.1371]
# I^2 = 33.6% [0.0%; 70.6%]; H = 1.23 [1.00; 1.84]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 10.55    7  0.1595

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/age_and_bmi_adjusted/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/age_and_bmi_adjusted/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/age_and_bmi_adjusted/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/age_and_bmi_adjusted/strict_het_test_before_outlier_extraction.RDS")

#We have a weird outlier

mr_post <- remove_outlier(mr_steiger_df, rucker) #no significant outliers. 

saveRDS(mr_post, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/age_and_bmi_adjusted//instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure                    method nsnp            b        se      pval
# 1    exposure    outcome outcome exposure                  MR Egger    7 -0.137009511 0.5312429 0.8067599
# 2    exposure    outcome outcome exposure           Weighted median    7  0.052872007 0.1958843 0.7872260
# 3    exposure    outcome outcome exposure Inverse variance weighted    7  0.003754940 0.1492439 0.9799275
# 4    exposure    outcome outcome exposure               Simple mode    7 -0.034408228 0.2939480 0.9106361
# 5    exposure    outcome outcome exposure             Weighted mode    7 -0.006700779 0.2936470 0.9825345

# #Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method   Estimate         SE      CI_low     CI_upp         P
# 1  Egger fixed effects 0.01210544 0.04107197 -0.06839415 0.09260503 0.7681946
# 2 Egger random effects 0.01210544 0.04360623 -0.07336121 0.09757209 0.3906566
# 
# [[1]]$Q
# Method          Q df         P
# 1   Q_ivw 5.72293531  6 0.4549315
# 2 Q_egger 5.63606526  5 0.3432564
# 3  Q_diff 0.08687005  1 0.7681946
# 
# [[1]]$res
# [1] "A"

strict_het <- compute_strict_het(mr_post) 

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 < 0.0001 [0.0000; 0.6729]; tau = 0.0021 [0.0000; 0.8203]
# I^2 = 0.0% [0.0%; 70.8%]; H = 1.00 [1.00; 1.85]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 5.53    6  0.4779

#That is quite good:

mr_post$labels <- NA
mr_post$units.exposure <- "SD"
mr_post$units.outcome <- "LOR"

plotio=mr_plots(mr_post)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/age_and_bmi_adjusted/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/age_and_bmi_adjusted/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/age_and_bmi_adjusted/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/age_and_bmi_adjusted/strict_het_test_after_outlier_extraction.RDS")

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

dir.create("output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/venkatesh")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/venkatesh/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp           b        se      pval
# 1    exposure    outcome outcome exposure                  MR Egger    9 0.006225337 0.4634974 0.9896585
# 2    exposure    outcome outcome exposure           Weighted median    9 0.218802523 0.1798238 0.2236953
# 3    exposure    outcome outcome exposure Inverse variance weighted    9 0.188472927 0.1384412 0.1733895
# 4    exposure    outcome outcome exposure               Simple mode    9 0.224181801 0.2724495 0.4344398
# 5    exposure    outcome outcome exposure             Weighted mode    9 0.251060735 0.2164081 0.2794461

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method  Estimate         SE      CI_low     CI_upp         P
# 1  Egger fixed effects 0.0166462 0.03423171 -0.05044673 0.08373912 0.6267686
# 2 Egger random effects 0.0166462 0.04017334 -0.06209209 0.09538449 0.3393055
# 
# [[1]]$Q
# Method         Q df         P
# 1   Q_ivw 9.8773461  8 0.2737404
# 2 Q_egger 9.6408778  7 0.2098542
# 3  Q_diff 0.2364683  1 0.6267686
# 
# [[1]]$res
# [1] "A"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_steiger_df) #they disagree, but still

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0.0469 [0.0000; 0.7175]; tau = 0.2166 [0.0000; 0.8471]
# I^2 = 14.8% [0.0%; 56.9%]; H = 1.08 [1.00; 1.52]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 9.39    8  0.3106

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/venkatesh/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/venkatesh/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/venkatesh/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/venkatesh/strict_het_test_before_outlier_extraction.RDS")

#This is great, except for one variant highlighted in LOO?

mr_post <- remove_outlier(mr_steiger_df, rucker) #no significant outlier

#saveRDS(mr_post, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/venkatesh/instruments_after_outlier_removal_df.RDS")

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
# ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/venkatesh/plots_after_outlier_extraction.svg", height=16, width=16, units="in")
# 
# saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/venkatesh/mr_res_after_outlier_extraction.RDS")
# saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/venkatesh/sensitivity_test_after_outlier_extraction.RDS")
# saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_asat_ratio/venkatesh/strict_het_test_after_outlier_extraction.RDS")