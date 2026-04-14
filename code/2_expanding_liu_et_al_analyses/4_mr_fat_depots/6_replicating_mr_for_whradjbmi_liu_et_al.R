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

liu_pcos_doctor <- fread("output/2_replicating_liu_et_al/2_mr/1_expo_outcome_df/whradjbmi_liu_PCOS.txt")
liu_pcos_broad <- fread("output/2_replicating_liu_et_al/2_mr/1_expo_outcome_df/whradjbmi_liu_PCOS (Broad).txt")
liu_pcos_consortium <- fread("output/2_replicating_liu_et_al/2_mr/1_expo_outcome_df/whradjbmi_liu_PCOS (Consortium).txt")
liu_pcos_meta_analysis <- fread("output/2_replicating_liu_et_al/2_mr/1_expo_outcome_df/whradjbmi_liu_PCOS (Day et al).txt")
liu_pcos_adj_age <- fread("output/2_replicating_liu_et_al/2_mr/1_expo_outcome_df/whradjbmi_liu_PCOS (adj age).txt")
liu_pcos_adj_age_bmi <- fread("output/2_replicating_liu_et_al/2_mr/1_expo_outcome_df/whradjbmi_liu_PCOS (adj age+BMI).txt")
liu_pcos_venkatesh <- fread("output/2_replicating_liu_et_al/2_mr/1_expo_outcome_df/whradjbmi_liu_PCOS (Venkatesh et al).txt")

liu_amenorrhea<- fread("output/2_replicating_liu_et_al/2_mr/1_expo_outcome_df/whradjbmi_liu_Amenorrhea.txt")
liu_oligomenorrhea<- fread("output/2_replicating_liu_et_al/2_mr/1_expo_outcome_df/whradjbmi_liu_Oligomenorrhea.txt")
liu_menorrhagia<- fread("output/2_replicating_liu_et_al/2_mr/1_expo_outcome_df/whradjbmi_liu_Menorrhagia.txt")
liu_anovulation<- fread("output/2_replicating_liu_et_al/2_mr/1_expo_outcome_df/whradjbmi_liu_Anovulation.txt")
liu_hirsutism<- fread("output/2_replicating_liu_et_al/2_mr/1_expo_outcome_df/whradjbmi_liu_Hirsutism.txt")

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

mr_df$id.exposure <- "WHRadjBMI"
mr_df$id.outcome <- "PCOS (meta-analysis)"
mr_df$exposure <- "WHRadjBMI"
mr_df$outcome <- "PCOS (meta-analysis)"
mr_df$units.exposure <- "SD"
mr_df$units.outcome <- "LOR"

mr_steiger_df <- TwoSampleMR::steiger_filtering(mr_df)
mr_steiger_df <- mr_steiger_df[which(mr_steiger_df$steiger_dir == TRUE),] #248
mr_steiger_df$mr_keep <- TRUE

#Let's save the data:

dir.create("output/2_replicating_liu_et_al/2_mr/2_mr/")
dir.create("output/2_replicating_liu_et_al/2_mr/2_mr/whradjbmi/")
dir.create("output/2_replicating_liu_et_al/2_mr/2_mr/whradjbmi/meta_analysis/")
saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al/2_mr/2_mr/whradjbmi/meta_analysis/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure           id.outcome              outcome  exposure                    method nsnp            b         se      pval
# 1   WHRadjBMI PCOS (meta-analysis) PCOS (meta-analysis) WHRadjBMI                  MR Egger  248 -0.003983120 0.20064501 0.9841779
# 2   WHRadjBMI PCOS (meta-analysis) PCOS (meta-analysis) WHRadjBMI           Weighted median  248 -0.044186875 0.15585011 0.7767770
# 3   WHRadjBMI PCOS (meta-analysis) PCOS (meta-analysis) WHRadjBMI Inverse variance weighted  248 -0.002750404 0.09005898 0.9756364
# 4   WHRadjBMI PCOS (meta-analysis) PCOS (meta-analysis) WHRadjBMI               Simple mode  248 -0.069186691 0.35504917 0.8456587
# 5   WHRadjBMI PCOS (meta-analysis) PCOS (meta-analysis) WHRadjBMI             Weighted mode  248 -0.178014675 0.22657335 0.4328067

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [1]]$intercept
# Method     Estimate          SE       CI_low      CI_upp         P
# 1  Egger fixed effects 3.598153e-05 0.004861501 -0.009492386 0.009564349 0.9940947
# 2 Egger random effects 3.598153e-05 0.004861501 -0.009492386 0.009564349 0.4970473
# 
# [[1]]$Q
# Method            Q  df         P
# 1   Q_ivw 2.122716e+02 247 0.9465224
# 2 Q_egger 2.122715e+02 246 0.9413048
# 3  Q_diff 4.726885e-05   1 0.9945144
# 
# [[1]]$res
# [1] "A"

strict_het <- compute_strict_het(mr_steiger_df)

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al/2_mr/2_mr/whradjbmi/meta_analysis/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al/2_mr/2_mr/whradjbmi/meta_analysis/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al/2_mr/2_mr/whradjbmi/meta_analysis//sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al/2_mr/2_mr/whradjbmi/meta_analysis//strict_het_test_before_outlier_extraction.RDS")

#This is a barran waste land. Nothing here.

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

dir.create("output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/")
dir.create("output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/doctor")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/doctor/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp         b         se         pval
# 1    exposure    outcome outcome exposure                  MR Egger  251 0.6015849 0.19338468 0.0020830223
# 2    exposure    outcome outcome exposure           Weighted median  251 0.4411334 0.14397217 0.0021838510
# 3    exposure    outcome outcome exposure Inverse variance weighted  251 0.3161276 0.08992647 0.0004390897
# 4    exposure    outcome outcome exposure               Simple mode  251 0.4973168 0.35437657 0.1617512440
# 5    exposure    outcome outcome exposure             Weighted mode  251 0.5470981 0.19738843 0.0059960930

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method     Estimate          SE      CI_low      CI_upp          P
# 1  Egger fixed effects -0.008431924 0.004920293 -0.01807552 0.001211673 0.08658317
# 2 Egger random effects -0.008431924 0.005061997 -0.01835326 0.001489407 0.95211649
# 
# [[1]]$Q
# Method         Q  df          P
# 1   Q_ivw 266.48561 250 0.22623331
# 2 Q_egger 263.54883 249 0.25161884
# 3  Q_diff   2.93678   1 0.08658317
# 
# [[1]]$res
# [1] "A"

strict_het <- compute_strict_het(mr_steiger_df)

#That is quite good:

mr_steiger_df$labels <- NA
mr_steiger_df$units.exposure <- "SD"
mr_steiger_df$units.outcome <- "LOR"

plotio=mr_plots(mr_steiger_df)
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/doctor/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/doctor/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/doctor/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/doctor/strict_het_test_before_outlier_extraction.RDS")

#A bit of asymmetry due to one variant that is highlighted in the funnel plot. I was close to not remove it, but since both loo and funnel show it, I will do it in this case.

mr_post <- remove_outlier(mr_steiger_df, rucker)

saveRDS(mr_post, "output/2_replicating_liu_et_al/2_mr/2_mr/whradjbmi/doctor/instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure                    method nsnp         b         se         pval
# 1    exposure    outcome outcome exposure                  MR Egger  237 0.5365852 0.19229436 0.0056958937
# 2    exposure    outcome outcome exposure           Weighted median  237 0.4408212 0.14563923 0.0024715511
# 3    exposure    outcome outcome exposure Inverse variance weighted  237 0.3404754 0.08977677 0.0001491586
# 4    exposure    outcome outcome exposure               Simple mode  237 0.4917318 0.37331719 0.1890504560
# 5    exposure    outcome outcome exposure             Weighted mode  237 0.5550516 0.21824702 0.0116230479

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method     Estimate          SE     CI_low     CI_upp         P
# 1  Egger fixed effects -0.005793727 0.004425886 -0.0144683 0.00288085 0.1905158
# 2 Egger random effects -0.005793727 0.004425886 -0.0144683 0.00288085 0.9047421
# 
# [[1]]$Q
# Method         Q  df         P
# 1   Q_ivw 183.71710 236 0.9950622
# 2 Q_egger 182.38713 235 0.9954021
# 3  Q_diff   1.32997   1 0.2488105
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/doctor/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/doctor/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/doctor/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/doctor/strict_het_test_after_outlier_extraction.RDS")

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

dir.create("output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/broad")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/broad/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp          b         se         pval
# 1    exposure    outcome outcome exposure                  MR Egger  251  0.6269412 0.14622964 2.587105e-05
# 2    exposure    outcome outcome exposure           Weighted median  251  0.3269547 0.10760162 2.377061e-03
# 3    exposure    outcome outcome exposure Inverse variance weighted  251  0.2230937 0.06863987 1.153215e-03
# 4    exposure    outcome outcome exposure               Simple mode  251 -0.0621571 0.28406199 8.269725e-01
# 5    exposure    outcome outcome exposure             Weighted mode  251  0.3425071 0.15348868 2.653628e-02

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method    Estimate          SE      CI_low       CI_upp            P
# 1  Egger fixed effects -0.01190013 0.003373207 -0.01851149 -0.005288762 0.0004189689
# 2 Egger random effects -0.01190013 0.003822747 -0.01939257 -0.004407678 0.9990739487
# 
# [[1]]$Q
# Method         Q  df            P
# 1   Q_ivw 332.23547 250 0.0003803720
# 2 Q_egger 319.78983 249 0.0016277677
# 3  Q_diff  12.44564   1 0.0004189689
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/broad/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/broad/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/broad/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/broad/strict_het_test_before_outlier_extraction.RDS")

#Pleitropy detected

mr_post <- remove_outlier(mr_steiger_df, rucker)
mr_post <- remove_outlier(mr_post, rucker) #ran again until outliers are removed entirely
mr_post=mr_post[which(!(mr_post$SNP=="rs72959041")),] #variant from LOO

saveRDS(mr_post, "output/2_replicating_liu_et_al/2_mr/2_mr/whradjbmi/broad/instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure                    method nsnp            b         se       pval
# 1    exposure    outcome outcome exposure                  MR Egger  226  0.344216288 0.17897927 0.05572037
# 2    exposure    outcome outcome exposure           Weighted median  226  0.080338621 0.10420307 0.44071801
# 3    exposure    outcome outcome exposure Inverse variance weighted  226  0.137231532 0.06705659 0.04070687
# 4    exposure    outcome outcome exposure               Simple mode  226 -0.069228818 0.27111058 0.79868397
# 5    exposure    outcome outcome exposure             Weighted mode  226  0.006318404 0.21691494 0.97678794

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method     Estimate          SE      CI_low      CI_upp         P
# 1  Egger fixed effects -0.005474804 0.003840953 -0.01300293 0.002053326 0.1540484
# 2 Egger random effects -0.005474804 0.003840953 -0.01300293 0.002053326 0.9229758
# 
# [[1]]$Q
# Method          Q  df         P
# 1   Q_ivw 173.089390 225 0.9957740
# 2 Q_egger 171.533567 224 0.9962284
# 3  Q_diff   1.555823   1 0.2122778
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/broad/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/broad/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/broad/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/broad/strict_het_test_after_outlier_extraction.RDS")

###############################################################
#Let's first analyse the data for the whole set for bmi - PCOS# consortium definition
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

dir.create("output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/consortium")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/consortium/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp          b         se        pval
# 1    exposure    outcome outcome exposure                  MR Egger  251 0.17890810 0.06230899 0.004439684
# 2    exposure    outcome outcome exposure           Weighted median  251 0.12465028 0.04108882 0.002415913
# 3    exposure    outcome outcome exposure Inverse variance weighted  251 0.07745389 0.02858295 0.006732518
# 4    exposure    outcome outcome exposure               Simple mode  251 0.03990223 0.09709876 0.681465074
# 5    exposure    outcome outcome exposure             Weighted mode  251 0.12947003 0.04833759 0.007886751

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method     Estimate          SE       CI_low        CI_upp          P
# 1  Egger fixed effects -0.002969866 0.001231171 -0.005382918 -0.0005568147 0.01585536
# 2 Egger random effects -0.002969866 0.001622745 -0.006150388  0.0002106555 0.96638622
# 
# [[1]]$Q
# Method          Q  df            P
# 1   Q_ivw 438.395297 250 1.762065e-12
# 2 Q_egger 432.576453 249 4.728196e-12
# 3  Q_diff   5.818845   1 1.585536e-02
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/consortium/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/consortium/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/consortium/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/consortium/strict_het_test_before_outlier_extraction.RDS")

#Very small heterogeneity, overall. Small outlier in loo and funnel, so let's check it out.

#STEP 5: removing outliers!

mr_post <- remove_outlier(mr_steiger_df, rucker)
mr_post <- remove_outlier(mr_post, rucker) #ran again until outliers are removed entirely
#mr_post=mr_post[which(!(mr_post$SNP=="rs72959041")),] #variant from LOO
  
saveRDS(mr_post, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/consortium/instruments_after_outlier_removal_df.RDS")
  
#STEP 6: let's run the results for real now:
  
mr_res <- TwoSampleMR::mr(mr_post)
  
# id.exposure id.outcome outcome exposure                    method nsnp          b         se         pval
# 1    exposure    outcome outcome exposure                  MR Egger  215 0.18848259 0.05030728 2.307704e-04
# 2    exposure    outcome outcome exposure           Weighted median  215 0.12948240 0.04234968 2.232222e-03
# 3    exposure    outcome outcome exposure Inverse variance weighted  215 0.09404616 0.02318858 4.998223e-05
# 4    exposure    outcome outcome exposure               Simple mode  215 0.04604801 0.09762489 6.376336e-01
# 5    exposure    outcome outcome exposure             Weighted mode  215 0.13346250 0.05248799 1.170592e-02

# #Let's check the plot:
#  
rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method     Estimate          SE       CI_low        CI_upp          P
# 1  Egger fixed effects -0.002790107 0.001239573 -0.005219624 -0.0003605891 0.02439428
# 2 Egger random effects -0.002790107 0.001239573 -0.005219624 -0.0003605891 0.98780286
# 
# [[1]]$Q
# Method         Q  df          P
# 1   Q_ivw 192.59206 214 0.85058700
# 2 Q_egger 188.11753 213 0.88940212
# 3  Q_diff   4.47453   1 0.03440367
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/consortium/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/consortium/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/consortium/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/consortium/strict_het_test_after_outlier_extraction.RDS")

###############################################################
#Let's first analyse the data for the whole set for bmi - PCOS# age adjusted
###############################################################

#STEP 1: let's use only those that have F>10

mr_df <- liu_pcos_adj_age[which(((liu_pcos_adj_age$beta.exposure^2)/(liu_pcos_adj_age$se.exposure^2)) > 10),] #250
liu_pcos_doctor_match=liu_pcos_doctor[which(liu_pcos_doctor$SNP%in%mr_df$SNP),]
mr_df=mr_df[which(mr_df$SNP%in%liu_pcos_doctor_match$SNP),]
mr_df=mr_df[order(match(mr_df$SNP, liu_pcos_doctor_match$SNP)),]
print(length(mr_df$SNP==liu_pcos_doctor_match$SNP))

mr_df$eaf.outcome <- liu_pcos_doctor_match$eaf.outcome #we do not have allele frequencies, so we will have to approximate them to this one

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

dir.create("output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/age_adjusted")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/age_adjusted/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp         b         se        pval
# 1    exposure    outcome outcome exposure                  MR Egger  247 0.2278567 0.16395032 0.165853861
# 2    exposure    outcome outcome exposure           Weighted median  247 0.2816386 0.12569669 0.025050650
# 3    exposure    outcome outcome exposure Inverse variance weighted  247 0.2057067 0.07427402 0.005613128
# 4    exposure    outcome outcome exposure               Simple mode  247 0.2744618 0.29000661 0.344874046
# 5    exposure    outcome outcome exposure             Weighted mode  247 0.2965026 0.16342979 0.070857332

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method      Estimate          SE       CI_low      CI_upp         P
# 1  Egger fixed effects -0.0006500846 0.004155037 -0.008793807 0.007493638 0.8756728
# 2 Egger random effects -0.0006500846 0.004287502 -0.009053434 0.007753265 0.5602579
# 
# [[1]]$Q
# Method            Q  df         P
# 1   Q_ivw 260.89498495 246 0.2456204
# 2 Q_egger 260.87050616 245 0.2320846
# 3  Q_diff   0.02447879   1 0.8756728
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/age_adjusted/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/age_adjusted/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/age_adjusted/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/age_adjusted/strict_het_test_before_outlier_extraction.RDS")

#No trace of heterogeneity or pleiotropy

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

dir.create("output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/age_and_bmi_adjusted")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/age_and_bmi_adjusted/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp         b         se       pval
# 1    exposure    outcome outcome exposure                  MR Egger  247 0.1779329 0.20262166 0.38072079
# 2    exposure    outcome outcome exposure           Weighted median  247 0.2456212 0.16230243 0.13018946
# 3    exposure    outcome outcome exposure Inverse variance weighted  247 0.2359715 0.09168183 0.01005859
# 4    exposure    outcome outcome exposure               Simple mode  247 0.1552172 0.36611927 0.67197099
# 5    exposure    outcome outcome exposure             Weighted mode  247 0.2159289 0.23291966 0.35480756

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method    Estimate          SE       CI_low     CI_upp         P
# 1  Egger fixed effects 0.001702232 0.004885018 -0.007872227 0.01127669 0.7274949
# 2 Egger random effects 0.001702232 0.005297111 -0.008679914 0.01208438 0.3739722
# 
# [[1]]$Q
# Method           Q  df          P
# 1   Q_ivw 288.2006079 246 0.03339024
# 2 Q_egger 288.0791837 245 0.03058975
# 3  Q_diff   0.1214242   1 0.72749494
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/age_and_bmi_adjusted/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/age_and_bmi_adjusted/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/age_and_bmi_adjusted/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/age_and_bmi_adjusted/strict_het_test_before_outlier_extraction.RDS")

#Slight heterogeneity and asymmetry, unlike the adj for age. 
#Removing outliers won't change a thing, but still

mr_post <- remove_outlier(mr_steiger_df, rucker)

saveRDS(mr_post, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/age_and_bmi_adjusted/instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure                    method nsnp         b         se         pval
# 1    exposure    outcome outcome exposure                  MR Egger  227 0.2811264 0.20112241 0.1635534415
# 2    exposure    outcome outcome exposure           Weighted median  227 0.2463114 0.16678861 0.1397324753
# 3    exposure    outcome outcome exposure Inverse variance weighted  227 0.3233728 0.09041888 0.0003483771
# 4    exposure    outcome outcome exposure               Simple mode  227 0.1967259 0.34472564 0.5687878044
# 5    exposure    outcome outcome exposure             Weighted mode  227 0.2349936 0.20767517 0.2590254402

# #Let's check the plot:
#  
rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method   Estimate          SE       CI_low      CI_upp         P
# 1  Egger fixed effects 0.00120303 0.004346277 -0.007315516 0.009721575 0.7819371
# 2 Egger random effects 0.00120303 0.004346277 -0.007315516 0.009721575 0.3909686
# 
# [[1]]$Q
# Method            Q  df         P
# 1   Q_ivw 162.45408850 226 0.9995067
# 2 Q_egger 162.39878935 225 0.9994186
# 3  Q_diff   0.05529915   1 0.8140863
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/age_and_bmi_adjusted/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/age_and_bmi_adjusted/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/age_and_bmi_adjusted/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/age_and_bmi_adjusted/strict_het_test_after_outlier_extraction.RDS")

###############################################################
#Let's first analyse the data for the whole set for bmi - PCOS# Venkatesh
###############################################################

#STEP 1: let's use only those that have F>10

mr_df <- liu_pcos_venkatesh[which(((liu_pcos_venkatesh$beta.exposure^2)/(liu_pcos_venkatesh$se.exposure^2)) > 10),] #256

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

dir.create("output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/venkatesh")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/venkatesh/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp         b        se        pval
# 1    exposure    outcome outcome exposure                  MR Egger  256 0.2756343 0.1707162 0.107644063
# 2    exposure    outcome outcome exposure           Weighted median  256 0.3870177 0.1471730 0.008546531
# 3    exposure    outcome outcome exposure Inverse variance weighted  256 0.2236369 0.0824527 0.006681713
# 4    exposure    outcome outcome exposure               Simple mode  256 0.3275596 0.2863953 0.253807168
# 5    exposure    outcome outcome exposure             Weighted mode  256 0.3606200 0.1444886 0.013198559

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method     Estimate          SE      CI_low      CI_upp         P
# 1  Egger fixed effects -0.001613932 0.004332621 -0.01010571 0.006877849 0.7095153
# 2 Egger random effects -0.001613932 0.004637354 -0.01070298 0.007475114 0.6360907
# 
# [[1]]$Q
# Method           Q  df          P
# 1   Q_ivw 291.1252185 255 0.05949921
# 2 Q_egger 290.9864569 254 0.05516118
# 3  Q_diff   0.1387616   1 0.70951529
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/venkatesh/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/venkatesh/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/venkatesh/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/venkatesh/strict_het_test_before_outlier_extraction.RDS")

#Some outliers in LOO, but not strong enough. Funnel looks very good and the rest of sensitivity tests look good.

# mr_post <- remove_outlier(mr_steiger_df, rucker)
# 
# saveRDS(mr_post, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/venkatesh/instruments_after_outlier_removal_df.RDS")
# 
# #STEP 6: let's run the results for real now:
# 
# mr_res <- TwoSampleMR::mr(mr_post)
# 
# # id.exposure id.outcome outcome exposure                    method nsnp         b         se        pval
# # 1    exposure    outcome outcome exposure                  MR Egger  182 0.3850080 0.21755072 0.078463714
# # 2    exposure    outcome outcome exposure           Weighted median  182 0.4687462 0.17700812 0.008093063
# # 3    exposure    outcome outcome exposure Inverse variance weighted  182 0.3182588 0.09863094 0.001251985
# # 4    exposure    outcome outcome exposure               Simple mode  182 0.4537847 0.30790687 0.142279873
# # 5    exposure    outcome outcome exposure             Weighted mode  182 0.4537847 0.19299455 0.019784703
# 
# # #Let's check the plot:
# #  
# rucker <- TwoSampleMR::mr_rucker(mr_post)
# 
# # [[1]]$intercept
# # Method     Estimate         SE      CI_low      CI_upp         P
# # 1  Egger fixed effects -0.001876149 0.00481659 -0.01131649 0.007564194 0.6968929
# # 2 Egger random effects -0.001876149 0.00481659 -0.01131649 0.007564194 0.6515535
# # 
# # [[1]]$Q
# # Method           Q  df         P
# # 1   Q_ivw 140.6969867 181 0.9882127
# # 2 Q_egger 140.5784913 180 0.9866701
# # 3  Q_diff   0.1184955   1 0.7306720
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
# ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/venkatesh/plots_after_outlier_extraction.svg", height=16, width=16, units="in")
# 
# saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/venkatesh/mr_res_after_outlier_extraction.RDS")
# saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/venkatesh/sensitivity_test_after_outlier_extraction.RDS")
# saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/venkatesh/strict_het_test_after_outlier_extraction.RDS")

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

dir.create("output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/amenorrhea")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/amenorrhea/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp          b         se      pval
# 1    exposure    outcome outcome exposure                  MR Egger  251 0.17054832 0.16656307 0.3068638
# 2    exposure    outcome outcome exposure           Weighted median  251 0.03595596 0.14836751 0.8085137
# 3    exposure    outcome outcome exposure Inverse variance weighted  251 0.03596642 0.07606825 0.6363432
# 4    exposure    outcome outcome exposure               Simple mode  251 0.34838928 0.30747757 0.2582762
# 5    exposure    outcome outcome exposure             Weighted mode  251 0.21237709 0.17873702 0.2358783

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method     Estimate          SE      CI_low      CI_upp         P
# 1  Egger fixed effects -0.003940651 0.004230878 -0.01223302 0.004351718 0.3516454
# 2 Egger random effects -0.003940651 0.004338369 -0.01244370 0.004562397 0.8181468
# 
# [[1]]$Q
# Method           Q  df         P
# 1   Q_ivw 262.6805943 250 0.2784422
# 2 Q_egger 261.8130835 249 0.2761361
# 3  Q_diff   0.8675108   1 0.3516454
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/amenorrhea/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/amenorrhea/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/amenorrhea/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/amenorrhea/strict_het_test_before_outlier_extraction.RDS")

#Very small outlier, does not affect results at all. Skipping.

# mr_post <- remove_outlier(mr_steiger_df, rucker)
# 
# saveRDS(mr_post, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/amenorrhea/instruments_after_outlier_removal_df.RDS")
# 
# #STEP 6: let's run the results for real now:
# 
# mr_res <- TwoSampleMR::mr(mr_post)
# 
# # id.exposure id.outcome outcome exposure                    method nsnp            b         se      pval
# # 1    exposure    outcome outcome exposure                  MR Egger  180  0.347208145 0.22291614 0.1211107
# # 2    exposure    outcome outcome exposure           Weighted median  180  0.059938027 0.17200999 0.7274973
# # 3    exposure    outcome outcome exposure Inverse variance weighted  180 -0.006302846 0.09408723 0.9465902
# # 4    exposure    outcome outcome exposure               Simple mode  180  0.483134024 0.38585786 0.2121657
# # 5    exposure    outcome outcome exposure             Weighted mode  180  0.300950228 0.23553442 0.2029980
# 
# #Let's check the plot:
# 
# rucker <- TwoSampleMR::mr_rucker(mr_post)
# 
# # [[1]]$intercept
# # Method     Estimate          SE      CI_low        CI_upp          P
# # 1  Egger fixed effects -0.009111712 0.004391308 -0.01771852 -0.0005049066 0.03799185
# # 2 Egger random effects -0.009111712 0.004391308 -0.01771852 -0.0005049066 0.98100407
# # 
# # [[1]]$Q
# # Method          Q  df          P
# # 1   Q_ivw 129.573427 179 0.99791333
# # 2 Q_egger 126.513376 178 0.99870014
# # 3  Q_diff   3.060051   1 0.08023916
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
# ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/amenorrhea/plots_after_outlier_extraction.svg", height=16, width=16, units="in")
# 
# saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/amenorrhea/mr_res_after_outlier_extraction.RDS")
# saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/amenorrhea/sensitivity_test_after_outlier_extraction.RDS")
# saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/amenorrhea/strict_het_test_after_outlier_extraction.RDS")

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

dir.create("output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/oligomenorrhea")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/oligomenorrhea/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp         b        se         pval
# 1    exposure    outcome outcome exposure                  MR Egger  251 0.8134481 0.2476193 0.0011661034
# 2    exposure    outcome outcome exposure           Weighted median  251 0.5047723 0.1812313 0.0053488500
# 3    exposure    outcome outcome exposure Inverse variance weighted  251 0.3767356 0.1136850 0.0009201959
# 4    exposure    outcome outcome exposure               Simple mode  251 0.6438185 0.4145875 0.1217088181
# 5    exposure    outcome outcome exposure             Weighted mode  251 0.6072385 0.2414607 0.0125365356

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method    Estimate          SE      CI_low        CI_upp          P
# 1  Egger fixed effects -0.01277912 0.006161733 -0.02485589 -0.0007023412 0.03808411
# 2 Egger random effects -0.01277912 0.006447006 -0.02541502 -0.0001432169 0.97627036
# 
# [[1]]$Q
# Method          Q  df          P
# 1   Q_ivw 276.891130 250 0.11671524
# 2 Q_egger 272.589869 249 0.14569270
# 3  Q_diff   4.301262   1 0.03808411
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/oligomenorrhea/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/oligomenorrhea/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/oligomenorrhea/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/oligomenorrhea/strict_het_test_before_outlier_extraction.RDS")

#Very controlled..., let's run it just in case - I can see some stuff...

mr_post <- remove_outlier(mr_steiger_df, rucker) #removing outliers removes variants and makes it pleiotropic!!! We are not doing it.
#mr_post <- remove_outlier(mr_post, rucker) #Run twice
#mr_post <- remove_outlier(mr_post, rucker) #does not run. the best we can do honestly

saveRDS(mr_post, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/oligomenorrhea/instruments_after_outlier_removal_df.RDS")
 
#STEP 6: let's run the results for real now:
 
mr_res <- TwoSampleMR::mr(mr_post)
 
# id.exposure id.outcome outcome exposure                    method nsnp         b        se         pval
# 1    exposure    outcome outcome exposure                  MR Egger  229 0.7958091 0.2436946 0.0012612522
# 2    exposure    outcome outcome exposure           Weighted median  229 0.5184616 0.1864235 0.0054175790
# 3    exposure    outcome outcome exposure Inverse variance weighted  229 0.3796966 0.1125620 0.0007429423
# 4    exposure    outcome outcome exposure               Simple mode  229 0.6547538 0.4085350 0.1103882134
# 5    exposure    outcome outcome exposure             Weighted mode  229 0.6291493 0.2448054 0.0108068949

#Let's check the plot:
  
rucker <- TwoSampleMR::mr_rucker(mr_post)
 
# [[1]]$intercept
# Method   Estimate          SE     CI_low       CI_upp          P
# 1  Egger fixed effects -0.0123082 0.005262951 -0.0226234 -0.001993007 0.01935354
# 2 Egger random effects -0.0123082 0.005262951 -0.0226234 -0.001993007 0.99032323
# 
# [[1]]$Q
# Method          Q  df          P
# 1   Q_ivw 157.537183 228 0.99988614
# 2 Q_egger 153.830819 227 0.99994349
# 3  Q_diff   3.706364   1 0.05420536
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/oligomenorrhea/plots_after_outlier_extraction.svg", height=16, width=16, units="in")
 
saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/oligomenorrhea/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/oligomenorrhea/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/oligomenorrhea/strict_het_test_after_outlier_extraction.RDS")

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

dir.create("output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/menorrhagia")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/menorrhagia/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp           b         se         pval
# 1    exposure    outcome outcome exposure                  MR Egger  251  0.24622005 0.08433914 0.0038282433
# 2    exposure    outcome outcome exposure           Weighted median  251  0.17574749 0.04999327 0.0004390542
# 3    exposure    outcome outcome exposure Inverse variance weighted  251  0.12809347 0.03873979 0.0009446631
# 4    exposure    outcome outcome exposure               Simple mode  251 -0.01681149 0.16073408 0.9167834879
# 5    exposure    outcome outcome exposure             Weighted mode  251  0.27530911 0.07088473 0.0001316532

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method     Estimate          SE       CI_low        CI_upp          P
# 1  Egger fixed effects -0.003464184 0.001555697 -0.006513294 -0.0004150727 0.02596251
# 2 Egger random effects -0.003464184 0.002198703 -0.007773561  0.0008451942 0.94243616
# 
# [[1]]$Q
# Method          Q  df            P
# 1   Q_ivw 502.331322 250 4.263006e-19
# 2 Q_egger 497.372808 249 1.056731e-18
# 3  Q_diff   4.958514   1 2.596251e-02
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/menorrhagia/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/menorrhagia/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/menorrhagia/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/menorrhagia/strict_het_test_before_outlier_extraction.RDS")

#Pleiotropy! Let's remove the outliers

mr_post <- remove_outlier(mr_steiger_df, rucker) 
#mr_post <- remove_outlier(mr_post, rucker) #run twice with updated rucker ofc, but no outliers are found...

saveRDS(mr_post, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/menorrhagia/instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure                    method nsnp           b         se         pval
# 1    exposure    outcome outcome exposure                  MR Egger  215  0.22338829 0.06375316 5.586237e-04
# 2    exposure    outcome outcome exposure           Weighted median  215  0.18058822 0.05112216 4.116753e-04
# 3    exposure    outcome outcome exposure Inverse variance weighted  215  0.14835549 0.02970074 5.883304e-07
# 4    exposure    outcome outcome exposure               Simple mode  215 -0.04739505 0.15090556 7.537731e-01
# 5    exposure    outcome outcome exposure             Weighted mode  215  0.27739766 0.07383375 2.216852e-04

#Let's check the plot:
#  
rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method     Estimate          SE       CI_low      CI_upp         P
# 1  Egger fixed effects -0.002194183 0.001631319 -0.005391511 0.001003144 0.1786136
# 2 Egger random effects -0.002194183 0.001631319 -0.005391511 0.001003144 0.9106932
# 
# [[1]]$Q
# Method          Q  df         P
# 1   Q_ivw 210.059305 214 0.5633157
# 2 Q_egger 208.290186 213 0.5782928
# 3  Q_diff   1.769119   1 0.1834912
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/menorrhagia/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/menorrhagia/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/menorrhagia/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/menorrhagia/strict_het_test_after_outlier_extraction.RDS")

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

dir.create("output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/anovulation")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/anovulation/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp          b         se         pval
# 1    exposure    outcome outcome exposure                  MR Egger  251 0.65631084 0.17763102 0.0002705649
# 2    exposure    outcome outcome exposure           Weighted median  251 0.33592586 0.12514447 0.0072680885
# 3    exposure    outcome outcome exposure Inverse variance weighted  251 0.22366419 0.08222098 0.0065226414
# 4    exposure    outcome outcome exposure               Simple mode  251 0.05781672 0.33492461 0.8630850231
# 5    exposure    outcome outcome exposure             Weighted mode  251 0.47844528 0.18324384 0.0095739668

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method    Estimate          SE      CI_low       CI_upp           P
# 1  Egger fixed effects -0.01266756 0.003952034 -0.02041340 -0.004921713 0.001349096
# 2 Egger random effects -0.01266756 0.004626095 -0.02173454 -0.003600579 0.996911958
# 
# [[1]]$Q
# Method         Q  df            P
# 1   Q_ivw 351.45680 250 2.388406e-05
# 2 Q_egger 341.18268 249 9.274261e-05
# 3  Q_diff  10.27411   1 1.349096e-03
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/anovulation/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/anovulation/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/anovulation/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/anovulation/strict_het_test_before_outlier_extraction.RDS")

#WHRadjBMI is linked to anovulation

mr_post <- remove_outlier(mr_steiger_df, rucker) 
mr_post <- remove_outlier(mr_post, rucker) #run twice with updated rucker ofc, but no outliers are found...

saveRDS(mr_post, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/anovulation/instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure                    method nsnp          b         se         pval
# 1    exposure    outcome outcome exposure                  MR Egger  227 0.58099519 0.15887791 0.0003179368
# 2    exposure    outcome outcome exposure           Weighted median  227 0.27872517 0.11888559 0.0190535161
# 3    exposure    outcome outcome exposure Inverse variance weighted  227 0.15472193 0.07333744 0.0348819241
# 4    exposure    outcome outcome exposure               Simple mode  227 0.08102794 0.34647993 0.8153050489
# 5    exposure    outcome outcome exposure             Weighted mode  227 0.44389162 0.17660908 0.0126541500

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method    Estimate          SE      CI_low       CI_upp           P
# 1  Egger fixed effects -0.01250512 0.003677146 -0.01971219 -0.005298045 0.000671969
# 2 Egger random effects -0.01250512 0.003677146 -0.01971219 -0.005298045 0.999664015
# 
# [[1]]$Q
# Method          Q  df          P
# 1   Q_ivw 187.115844 226 0.97218919
# 2 Q_egger 177.968109 225 0.99090341
# 3  Q_diff   9.147736   1 0.00249026
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/anovulation/plots_after_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/anovulation/mr_res_after_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/anovulation/sensitivity_test_after_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/anovulation/strict_het_test_after_outlier_extraction.RDS")

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

dir.create("output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/hirsutism")

saveRDS(mr_steiger_df, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/hirsutism/instruments_after_steiger_df.RDS")

#STEP 4: let's run the analyses:

mr_res <- TwoSampleMR::mr(mr_steiger_df)

# id.exposure id.outcome outcome exposure                    method nsnp         b        se        pval
# 1    exposure    outcome outcome exposure                  MR Egger  251 0.2886729 0.3533258 0.414700229
# 2    exposure    outcome outcome exposure           Weighted median  251 0.6295227 0.2935961 0.032018526
# 3    exposure    outcome outcome exposure Inverse variance weighted  251 0.4510175 0.1610227 0.005095153
# 4    exposure    outcome outcome exposure               Simple mode  251 0.8143606 0.5947466 0.172147953
# 5    exposure    outcome outcome exposure             Weighted mode  251 0.4953822 0.3630578 0.173645964

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)

# [[1]]$intercept
# Method    Estimate          SE      CI_low     CI_upp         P
# 1  Egger fixed effects 0.004749545 0.009064493 -0.01301653 0.02251563 0.6002977
# 2 Egger random effects 0.004749545 0.009197476 -0.01327718 0.02277627 0.3027887
# 
# [[1]]$Q
# Method           Q  df         P
# 1   Q_ivw 256.6341631 250 0.3730224
# 2 Q_egger 256.3596159 249 0.3607687
# 3  Q_diff   0.2745472   1 0.6002977
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
ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/hirsutism/plots_before_outlier_extraction.svg", height=16, width=16, units="in")

saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/hirsutism/mr_res_before_outlier_extraction.RDS")
saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/hirsutism/sensitivity_test_before_outlier_extraction.RDS")
saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/hirsutism/strict_het_test_before_outlier_extraction.RDS")

#Quite good actually, we do not need anything else

# mr_post <- remove_outlier(mr_steiger_df, rucker) 
# #mr_post <- remove_outlier(mr_post, rucker) #run twice with updated rucker ofc, but no outliers are found...
# 
# saveRDS(mr_post, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/hirsutism/instruments_after_outlier_removal_df.RDS")
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
# ggsave(plot=plotio,filename = "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/hirsutism/plots_after_outlier_extraction.svg", height=16, width=16, units="in")
# 
# saveRDS(mr_res, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/hirsutism/mr_res_after_outlier_extraction.RDS")
# saveRDS(rucker, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/hirsutism/sensitivity_test_after_outlier_extraction.RDS")
# saveRDS(strict_het, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/hirsutism/strict_het_test_after_outlier_extraction.RDS")
# 
