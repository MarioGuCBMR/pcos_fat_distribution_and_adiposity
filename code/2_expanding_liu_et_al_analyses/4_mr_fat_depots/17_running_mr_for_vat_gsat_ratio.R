##############
#INTRODUCTION#
##############

#This code performs MR analyses between vat_gsat_ratioBMI and PCOS - replicating the results from Liu et al.

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

liu_pcos_doctor <- fread("output/2_replicating_liu_et_al/4_mr_fat_depots/1_expo_outcome_df/vat_gsat_ratio_PCOS.txt")
liu_pcos_broad <- fread("output/2_replicating_liu_et_al/4_mr_fat_depots/1_expo_outcome_df/vat_gsat_ratio_PCOS (Broad).txt")
liu_pcos_consortium <- fread("output/2_replicating_liu_et_al/4_mr_fat_depots/1_expo_outcome_df/vat_gsat_ratio_PCOS (Consortium).txt")
liu_pcos_meta_analysis <- fread("output/2_replicating_liu_et_al/4_mr_fat_depots/1_expo_outcome_df/vat_gsat_ratio_PCOS (Day et al).txt")
liu_pcos_venkatesh <- fread("output/2_replicating_liu_et_al/4_mr_fat_depots/1_expo_outcome_df/vat_gsat_ratio_PCOS (Venkatesh et al).txt")
liu_pcos_adj_age <- fread("output/2_replicating_liu_et_al/4_mr_fat_depots/1_expo_outcome_df/vat_gsat_ratio_PCOS (adj age).txt")
liu_pcos_adj_age_bmi <- fread("output/2_replicating_liu_et_al/4_mr_fat_depots/1_expo_outcome_df/vat_gsat_ratio_PCOS (adj age+BMI).txt")

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

mr_df$id.exposure <- "vat_gsat_ratio"
mr_df$id.outcome <- "PCOS (meta-analysis)"
mr_df$exposure <- "vat_gsat_ratio"
mr_df$outcome <- "PCOS (meta-analysis)"
mr_df$units.exposure <- "SD"
mr_df$units.outcome <- "LOR"

mr_steiger_df <- TwoSampleMR::steiger_filtering(mr_df)
mr_steiger_df <- mr_steiger_df[which(mr_steiger_df$steiger_dir == TRUE),] #19
mr_steiger_df$mr_keep <- TRUE

#Let's save the data:

dir.create("output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/")
dir.create("output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/vat_gsat_ratio/")
dir.create("output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/vat_gsat_ratio/meta_analysis/")
saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/vat_gsat_ratio/meta_analysis/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure           id.outcome              outcome       exposure                    method nsnp           b        se      pval
# 1 vat_gsat_ratio PCOS (meta-analysis) PCOS (meta-analysis) vat_gsat_ratio                  MR Egger    7  0.79893706 0.7698362 0.3469216
# 2 vat_gsat_ratio PCOS (meta-analysis) PCOS (meta-analysis) vat_gsat_ratio           Weighted median    7  0.01382877 0.2787338 0.9604309
# 3 vat_gsat_ratio PCOS (meta-analysis) PCOS (meta-analysis) vat_gsat_ratio Inverse variance weighted    7  0.10993897 0.2125205 0.6049399
# 4 vat_gsat_ratio PCOS (meta-analysis) PCOS (meta-analysis) vat_gsat_ratio               Simple mode    7 -0.07633350 0.4476597 0.8702084
# 5 vat_gsat_ratio PCOS (meta-analysis) PCOS (meta-analysis) vat_gsat_ratio             Weighted mode    7 -0.06356034 0.4269162 0.8865243

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method    Estimate         SE     CI_low     CI_upp         P
# 1  Egger fixed effects -0.05042758 0.04000356 -0.1288331 0.02797795 0.2074612
# 2 Egger random effects -0.05042758 0.04000356 -0.1288331 0.02797795 0.8962694
# 
# [[1]]$Q
# Method         Q df         P
# 1   Q_ivw 3.5954204  6 0.7312342
# 2 Q_egger 2.7283277  5 0.7417826
# 3  Q_diff 0.8670927  1 0.3517614
# 
# [[1]]$res
# [1] "A"

strict_het <- compute_strict_het(mr_steiger_df)

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0 [0.0000; 0.4636]; tau = 0 [0.0000; 0.6809]
# I^2 = 0.0% [0.0%; 70.8%]; H = 1.00 [1.00; 1.85]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 3.49    6  0.7454

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/vat_gsat_ratio/meta_analysis/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/vat_gsat_ratio/meta_analysis/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/vat_gsat_ratio/meta_analysis//sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/vat_gsat_ratio/meta_analysis//strict_het_test_before_outlier_extraction.RDS")

#Slight asymmetry, but very tame - results are what they are!

mr_post <- remove_outlier(mr_steiger_df, rucker) #no significant outliers!

#saveRDS(mr_post, "output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/vat_gsat_ratio/meta_analysis/instruments_after_outlier_removal_df.RDS")

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
# ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/meta_analysis/plots_after_outlier_extraction.svg", height=16, width=16, units="in")
# 
# saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/meta_analysis/mr_res_after_outlier_extraction.RDS")
# saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/meta_analysis/sensitivity_test_after_outlier_extraction.RDS")
# saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/meta_analysis/strict_het_test_after_outlier_extraction.RDS")

###############################################################
#Let's first analyse the data for the whole set for vat_gsat_ratio - PCOS#
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

dir.create("output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/")
dir.create("output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/doctor")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/doctor/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp         b        se      pval
# 1    exposure    outcome outcome exposure                  MR Egger    8 0.8245822 0.7241848 0.2982673
# 2    exposure    outcome outcome exposure           Weighted median    8 0.2876073 0.2589439 0.2667004
# 3    exposure    outcome outcome exposure Inverse variance weighted    8 0.4170215 0.2632130 0.1131141
# 4    exposure    outcome outcome exposure               Simple mode    8 0.1099325 0.4414480 0.8104886
# 5    exposure    outcome outcome exposure             Weighted mode    8 0.2552732 0.3954681 0.5391798

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method    Estimate         SE     CI_low     CI_upp         P
# 1  Egger fixed effects -0.03391984 0.03564198 -0.1037768 0.03593716 0.3412582
# 2 Egger random effects -0.03391984 0.05572549 -0.1431398 0.07530011 0.7286367
# 
# [[1]]$Q
# Method          Q df          P
# 1   Q_ivw 15.5724994  7 0.02932232
# 2 Q_egger 14.6668005  6 0.02301274
# 3  Q_diff  0.9056989  1 0.34125820
# 
# [[1]]$res
# [1] "B"

strict_het <- compute_strict_het(mr_steiger_df)

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0.2527 [0.0000; 2.4232]; tau = 0.5027 [0.0000; 1.5567]
# I^2 = 48.9% [0.0%; 77.2%]; H = 1.40 [1.00; 2.10]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 13.70    7  0.0568

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/doctor/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/doctor/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/doctor/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/doctor/strict_het_test_before_outlier_extraction.RDS")

#Leave one out highlight one outlier!

mr_post <- remove_outlier(mr_steiger_df, rucker) #outliers detected!

saveRDS(mr_post, "output/2_replicating_liu_et_al/4_mr_fat_depots/2_mr/vat_gsat_ratio/doctor/instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure                    method nsnp          b        se       pval
# 1    exposure    outcome outcome exposure                  MR Egger    6 1.39235762 0.4889447 0.04650481
# 2    exposure    outcome outcome exposure           Weighted median    6 0.29984258 0.2827286 0.28890285
# 3    exposure    outcome outcome exposure Inverse variance weighted    6 0.42297903 0.1950792 0.03014027
# 4    exposure    outcome outcome exposure               Simple mode    6 0.01567084 0.4406033 0.97300405
# 5    exposure    outcome outcome exposure             Weighted mode    6 0.06580641 0.4119774 0.87934361

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_post)
 
# [[1]]$intercept
# Method    Estimate          SE     CI_low      CI_upp            P
# 1  Egger fixed effects -0.08599801 0.009686265 -0.1049827 -0.06701328 6.786214e-19
# 2 Egger random effects -0.08599801 0.009686265 -0.1049827 -0.06701328 1.000000e+00
# 
# [[1]]$Q
# Method         Q df         P
# 1   Q_ivw 4.9120693  5 0.4267055
# 2 Q_egger 0.2372263  4 0.9934977
# 3  Q_diff 4.6748429  1 0.0306075
# 
# [[1]]$res
# [1] "A"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_post) #they disagree, but still

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0.0279 [0.0000; 1.0146]; tau = 0.1670 [0.0000; 1.0073]
# I^2 = 0.0% [0.0%; 74.6%]; H = 1.00 [1.00; 1.99]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 4.64    5  0.4614

#That is quite good:

mr_post$labels <- NA
mr_post$units.exposure <- "SD"
mr_post$units.outcome <- "LOR"
 
plotio=mr_plots(mr_post)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/doctor/plots_after_outlier_extraction.svg", height=16, width=16, units="in")
 
saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/doctor/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/doctor/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/doctor/strict_het_test_after_outlier_extraction.RDS")

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

dir.create("output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/broad")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/broad/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp          b        se      pval
# 1    exposure    outcome outcome exposure                  MR Egger    8  0.7614766 0.5295004 0.2004480
# 2    exposure    outcome outcome exposure           Weighted median    8  0.1130981 0.2078227 0.5863006
# 3    exposure    outcome outcome exposure Inverse variance weighted    8  0.1638590 0.2073844 0.4294566
# 4    exposure    outcome outcome exposure               Simple mode    8 -0.3774392 0.3688398 0.3402160
# 5    exposure    outcome outcome exposure             Weighted mode    8 -0.2729660 0.3761002 0.4915204

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method    Estimate         SE      CI_low      CI_upp          P
# 1  Egger fixed effects -0.04953103 0.02451839 -0.09758618 -0.00147587 0.04336696
# 2 Egger random effects -0.04953103 0.04061671 -0.12913832  0.03007626 0.88866785
# 
# [[1]]$Q
# Method        Q df           P
# 1   Q_ivw 20.54662  7 0.004502463
# 2 Q_egger 16.46558  6 0.011461627
# 3  Q_diff  4.08104  1 0.043366961
# 
# [[1]]$res
# [1] "D"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_steiger_df) #they disagree, but still

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0.2160 [0.0198; 1.4343]; tau = 0.4647 [0.1408; 1.1976]
# I^2 = 62.2% [18.6%; 82.5%]; H = 1.63 [1.11; 2.39]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 18.54    7  0.0098

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/broad/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/broad/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/broad/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/broad/strict_het_test_before_outlier_extraction.RDS")

#Asymmetry and heterogeneity introduced!

mr_post <- remove_outlier(mr_steiger_df, rucker)
 
saveRDS(mr_post, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/broad/instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure                    method nsnp           b        se       pval
# 1    exposure    outcome outcome exposure                  MR Egger    7  1.13380937 0.3355555 0.01969771
# 2    exposure    outcome outcome exposure           Weighted median    7 -0.00443158 0.1951960 0.98188700
# 3    exposure    outcome outcome exposure Inverse variance weighted    7  0.07027412 0.2006821 0.72620638
# 4    exposure    outcome outcome exposure               Simple mode    7 -0.44251434 0.3423444 0.24369219
# 5    exposure    outcome outcome exposure             Weighted mode    7 -0.42251446 0.4081079 0.34044053

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method    Estimate         SE     CI_low      CI_upp            P
# 1  Egger fixed effects -0.09367378 0.02179309 -0.1363875 -0.05096011 1.720937e-05
# 2 Egger random effects -0.09367378 0.02179309 -0.1363875 -0.05096011 9.999914e-01
# 
# [[1]]$Q
# Method         Q df            P
# 1   Q_ivw 14.910239  6 0.0209664636
# 2 Q_egger  3.175689  5 0.6729206209
# 3  Q_diff 11.734550  1 0.0006135042
# 
# [[1]]$res
# [1] "C"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_post) #they disagree, but still

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0.1531 [0.0000; 1.1742]; tau = 0.3913 [0.0000; 1.0836]
# I^2 = 56.5% [0.0%; 81.3%]; H = 1.52 [1.00; 2.31]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 13.80    6  0.0319

#That is quite good:

mr_post$labels <- NA
mr_post$units.exposure <- "SD"
mr_post$units.outcome <- "LOR"
 
plotio=mr_plots(mr_post)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/broad/plots_after_outlier_extraction.svg", height=16, width=16, units="in")
 
saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/broad/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/broad/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/broad/strict_het_test_after_outlier_extraction.RDS")

#Less SNP have made it worse...

################################################################
#Let's first analyse the data for the whole set for vat_gsat_ratio - PCOS# consortium definition
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

dir.create("output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/consortium")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/consortium/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp           b         se       pval
# 1    exposure    outcome outcome exposure                  MR Egger    8 0.306844564 0.22394425 0.21968496
# 2    exposure    outcome outcome exposure           Weighted median    8 0.074656279 0.06485881 0.24970823
# 3    exposure    outcome outcome exposure Inverse variance weighted    8 0.148622189 0.08099591 0.06651538
# 4    exposure    outcome outcome exposure               Simple mode    8 0.007880917 0.08526516 0.92894736
# 5    exposure    outcome outcome exposure             Weighted mode    8 0.047688221 0.07349492 0.53712459

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method   Estimate          SE      CI_low      CI_upp         P
# 1  Egger fixed effects -0.0129867 0.009009052 -0.03064412 0.004670713 0.1494386
# 2 Egger random effects -0.0129867 0.017054195 -0.04641231 0.020438903 0.7768196
# 
# [[1]]$Q
# Method         Q df           P
# 1   Q_ivw 23.578814  7 0.001350668
# 2 Q_egger 21.500842  6 0.001490604
# 3  Q_diff  2.077972  1 0.149438593
# 
# [[1]]$res
# [1] "B"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_steiger_df) #they disagree, but still

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0.0353 [0.0044; 0.2585]; tau = 0.1879 [0.0661; 0.5084]
# I^2 = 63.0% [20.5%; 82.8%]; H = 1.64 [1.12; 2.41]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 18.93    7  0.0084

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/consortium/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/consortium/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/consortium/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/consortium/strict_het_test_before_outlier_extraction.RDS")

#We have heterogeneity and asymmetry

#STEP 5: removing outliers!

mr_post <- remove_outlier(mr_steiger_df, rucker) #but not significant outliers!!
  
saveRDS(mr_post, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/consortium/instruments_after_outlier_removal_df.RDS")
  
#STEP 6: let's run the results for real now:
  
mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure                    method nsnp            b         se       pval
# 1    exposure    outcome outcome exposure                  MR Egger    6  0.319602936 0.13577969 0.07818723
# 2    exposure    outcome outcome exposure           Weighted median    6  0.028624499 0.06631805 0.66601430
# 3    exposure    outcome outcome exposure Inverse variance weighted    6  0.061718995 0.05008900 0.21787941
# 4    exposure    outcome outcome exposure               Simple mode    6 -0.027909603 0.11011583 0.81000701
# 5    exposure    outcome outcome exposure             Weighted mode    6 -0.007957259 0.10125760 0.94041131

# #Let's check the plot:
  
rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method    Estimate          SE      CI_low      CI_upp            P
# 1  Egger fixed effects -0.02160057 0.005742341 -0.03285536 -0.01034579 0.0001688085
# 2 Egger random effects -0.02160057 0.005742341 -0.03285536 -0.01034579 0.9999155958
# 
# [[1]]$Q
# Method        Q df          P
# 1   Q_ivw 5.307438  5 0.37952307
# 2 Q_egger 1.169691  4 0.88306307
# 3  Q_diff 4.137746  1 0.04193696
# 
# [[1]]$res
# [1] "A"

# #Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_post) #they disagree, but still

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0.0014 [0.0000; 0.0793]; tau = 0.0381 [0.0000; 0.2817]
# I^2 = 2.1% [0.0%; 75.1%]; H = 1.01 [1.00; 2.01]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 5.10    5  0.4032

#That is quite good:
  
mr_post$labels <- NA
mr_post$units.exposure <- "SD"
mr_post$units.outcome <- "LOR"
  
plotio=mr_plots(mr_post)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/consortium/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/consortium/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/consortium/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/consortium/strict_het_test_after_outlier_extraction.RDS")

#I think we might have now more asymmetry, despite heterogeneity not popping up...

#mr_post_2 <- remove_outlier(mr_post, rucker) #cannot detect it...

###############################################################
#Let's first analyse the data for the whole set for vat_gsat_ratio - PCOS# age adjusted
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

dir.create("output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/age_adjusted")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/age_adjusted/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp         b        se       pval
# 1    exposure    outcome outcome exposure                  MR Egger    7 0.7697126 0.7274049 0.33840114
# 2    exposure    outcome outcome exposure           Weighted median    7 0.4669580 0.2124856 0.02797773
# 3    exposure    outcome outcome exposure Inverse variance weighted    7 0.1822269 0.2436068 0.45443793
# 4    exposure    outcome outcome exposure               Simple mode    7 0.4768294 0.2682237 0.12577122
# 5    exposure    outcome outcome exposure             Weighted mode    7 0.5180181 0.2451394 0.07902722

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method    Estimate         SE     CI_low     CI_upp         P
# 1  Egger fixed effects -0.04574573 0.03274215 -0.1099192 0.01842770 0.1623681
# 2 Egger random effects -0.04574573 0.05321647 -0.1500481 0.05855664 0.8049996
# 
# [[1]]$Q
# Method         Q df          P
# 1   Q_ivw 15.160361  6 0.01904548
# 2 Q_egger 13.208329  5 0.02150278
# 3  Q_diff  1.952032  1 0.16236809
# 
# [[1]]$res
# [1] "B"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_steiger_df) #they disagree, but still

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0.2299 [0.0000; 2.2034]; tau = 0.4795 [0.0000; 1.4844]
# I^2 = 55.6% [0.0%; 80.9%]; H = 1.50 [1.00; 2.29]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 13.52    6  0.0355

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/age_adjusted//plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/age_adjusted//mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/age_adjusted//sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/age_adjusted//strict_het_test_before_outlier_extraction.RDS")

#We have a weird outlier! It does not affect things much, but let's see...

mr_post <- remove_outlier(mr_steiger_df, rucker) #no significant outliers. 

saveRDS(mr_post, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/age_adjusted//instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure                    method nsnp         b        se       pval
# 1    exposure    outcome outcome exposure                  MR Egger    6 0.7533092 0.5195790 0.22070322
# 2    exposure    outcome outcome exposure           Weighted median    6 0.4828559 0.2064917 0.01936759
# 3    exposure    outcome outcome exposure Inverse variance weighted    6 0.3265704 0.1824714 0.07350101
# 4    exposure    outcome outcome exposure               Simple mode    6 0.5010930 0.2806642 0.13426090
# 5    exposure    outcome outcome exposure             Weighted mode    6 0.5327598 0.2422878 0.07920637

# #Let's check the plot:
  
rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method   Estimate         SE     CI_low     CI_upp         P
# 1  Egger fixed effects -0.0337387 0.03302255 -0.0984617 0.03098431 0.3069293
# 2 Egger random effects -0.0337387 0.03833431 -0.1088726 0.04139516 0.8106022
# 
# [[1]]$Q
# Method        Q df         P
# 1   Q_ivw 6.434156  5 0.2662347
# 2 Q_egger 5.390313  4 0.2495406
# 3  Q_diff 1.043843  1 0.3069293
# 
# [[1]]$res
# [1] "A"

# #Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_post) 

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0.0249 [0.0000; 1.1801]; tau = 0.1579 [0.0000; 1.0863]
# I^2 = 18.2% [0.0%; 63.0%]; H = 1.11 [1.00; 1.64]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 6.12    5  0.2951

#That is quite good:

mr_post$labels <- NA
mr_post$units.exposure <- "SD"
mr_post$units.outcome <- "LOR"

plotio=mr_plots(mr_post)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/age_adjusted/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/age_adjusted/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/age_adjusted/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/age_adjusted/strict_het_test_after_outlier_extraction.RDS")

#Did not work...

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

dir.create("output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/age_and_bmi_adjusted")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/age_and_bmi_adjusted/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)
 
# id.exposure id.outcome outcome exposure                    method nsnp         b        se      pval
# 1    exposure    outcome outcome exposure                  MR Egger    7 0.7535386 0.8138842 0.3970105
# 2    exposure    outcome outcome exposure           Weighted median    7 0.2742089 0.2614912 0.2943460
# 3    exposure    outcome outcome exposure Inverse variance weighted    7 0.1574341 0.2682014 0.5572039
# 4    exposure    outcome outcome exposure               Simple mode    7 0.2638981 0.3630571 0.4946715
# 5    exposure    outcome outcome exposure             Weighted mode    7 0.2742450 0.3093582 0.4094724

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method    Estimate         SE     CI_low     CI_upp         P
# 1  Egger fixed effects -0.04635368 0.03848593 -0.1217847 0.02907735 0.2284227
# 2 Egger random effects -0.04635368 0.05949776 -0.1629671 0.07025978 0.7820345
# 
# [[1]]$Q
# Method         Q df          P
# 1   Q_ivw 13.400638  6 0.03709703
# 2 Q_egger 11.949981  5 0.03547962
# 3  Q_diff  1.450656  1 0.22842267
# 
# [[1]]$res
# [1] "B"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_steiger_df) #they disagree, but still

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0.2537 [0.0000; 2.3569]; tau = 0.5037 [0.0000; 1.5352]
# I^2 = 50.9% [0.0%; 79.1%]; H = 1.43 [1.00; 2.19]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 12.21    6  0.0574

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/age_and_bmi_adjusted/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/age_and_bmi_adjusted/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/age_and_bmi_adjusted/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/age_and_bmi_adjusted/strict_het_test_before_outlier_extraction.RDS")

#We have a weird outlier

mr_post <- remove_outlier(mr_steiger_df, rucker) #no significant outliers. 

saveRDS(mr_post, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/age_and_bmi_adjusted//instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure                    method nsnp         b        se      pval
# 1    exposure    outcome outcome exposure                  MR Egger    4 0.3128066 0.5477369 0.6255570
# 2    exposure    outcome outcome exposure           Weighted median    4 0.2761138 0.2681976 0.3032371
# 3    exposure    outcome outcome exposure Inverse variance weighted    4 0.2720467 0.2342946 0.2455887
# 4    exposure    outcome outcome exposure               Simple mode    4 0.3123184 0.3606540 0.4502071
# 5    exposure    outcome outcome exposure             Weighted mode    4 0.3106358 0.3129638 0.3940895

# #Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method     Estimate          SE      CI_low      CI_upp         P
# 1  Egger fixed effects -0.003463105 0.006826185 -0.01684218 0.009915972 0.6119257
# 2 Egger random effects -0.003463105 0.006826185 -0.01684218 0.009915972 0.6940371
# 
# [[1]]$Q
# Method           Q df         P
# 1   Q_ivw 0.059444770  3 0.9962133
# 2 Q_egger 0.052667043  2 0.9740102
# 3  Q_diff 0.006777727  1 0.9343867
# 
# [[1]]$res
# [1] "A"

strict_het <- compute_strict_het(mr_post) 

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0; tau = 0; I^2 = 0.0% [0.0%; 84.7%]; H = 1.00 [1.00; 2.56]
#   
#   Test of heterogeneity:
#     Q d.f. p-value
#   0.06    3  0.9963

#That is quite good:

mr_post$labels <- NA
mr_post$units.exposure <- "SD"
mr_post$units.outcome <- "LOR"

plotio=mr_plots(mr_post)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/age_and_bmi_adjusted/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/age_and_bmi_adjusted/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/age_and_bmi_adjusted/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/age_and_bmi_adjusted/strict_het_test_after_outlier_extraction.RDS")

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

dir.create("output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/venkatesh")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/venkatesh/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp         b        se       pval
# 1    exposure    outcome outcome exposure                  MR Egger    8 0.7363991 0.5188817 0.20563986
# 2    exposure    outcome outcome exposure           Weighted median    8 0.3929985 0.2057154 0.05608183
# 3    exposure    outcome outcome exposure Inverse variance weighted    8 0.2780195 0.1898616 0.14310455
# 4    exposure    outcome outcome exposure               Simple mode    8 0.2130920 0.3486266 0.56036241
# 5    exposure    outcome outcome exposure             Weighted mode    8 0.4721995 0.2286239 0.07774305

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method    Estimate         SE     CI_low     CI_upp         P
# 1  Egger fixed effects -0.04093361 0.03550012 -0.1105126 0.02864535 0.2488876
# 2 Egger random effects -0.04093361 0.04307629 -0.1253616 0.04349436 0.8290095
# 
# [[1]]$Q
# Method         Q df         P
# 1   Q_ivw 10.163756  7 0.1794819
# 2 Q_egger  8.834219  6 0.1831184
# 3  Q_diff  1.329537  1 0.2488876
# 
# [[1]]$res
# [1] "A"

#Here there is heterogeneity and pleitropy

strict_het <- compute_strict_het(mr_steiger_df) #they disagree, but still

# Quantifying heterogeneity (with 95%-CIs):
#   tau^2 = 0.0691 [0.0000; 2.2318]; tau = 0.2628 [0.0000; 1.4939]
# I^2 = 27.3% [0.0%; 67.3%]; H = 1.17 [1.00; 1.75]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 9.62    7  0.2110

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/venkatesh/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/venkatesh/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/venkatesh/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/venkatesh/strict_het_test_before_outlier_extraction.RDS")

#This is great, except for one variant highlighted in LOO?

mr_post <- remove_outlier(mr_steiger_df, rucker) #no significant outlier

#saveRDS(mr_post, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/venkatesh/instruments_after_outlier_removal_df.RDS")

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
# ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/venkatesh/plots_after_outlier_extraction.svg", height=16, width=16, units="in")
# 
# saveRDS(mr_res, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/venkatesh/mr_res_after_outlier_extraction.RDS")
# saveRDS(rucker, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/venkatesh/sensitivity_test_after_outlier_extraction.RDS")
# saveRDS(strict_het, "output/2_replicating_liu_et_al//4_mr_fat_depots/2_mr/vat_gsat_ratio/venkatesh/strict_het_test_after_outlier_extraction.RDS")