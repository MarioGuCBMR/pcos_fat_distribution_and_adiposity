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

liu_pcos_doctor <- fread("output/2_replicating_liu_et_al/2_mr/1_expo_outcome_df/whradjbmi_liu_PCOS.txt")
liu_pcos_broad <- fread("output/2_replicating_liu_et_al/2_mr/1_expo_outcome_df/whradjbmi_liu_PCOS (Broad).txt")
liu_pcos_consortium <- fread("output/2_replicating_liu_et_al/2_mr/1_expo_outcome_df/whradjbmi_liu_PCOS (Consortium).txt")
liu_pcos_meta_analysis <- fread("output/2_replicating_liu_et_al/2_mr/1_expo_outcome_df/whradjbmi_liu_PCOS (Day et al).txt")
liu_pcos_adj_age <- fread("output/2_replicating_liu_et_al/2_mr/1_expo_outcome_df/whradjbmi_liu_PCOS (adj age).txt")
liu_pcos_adj_age_bmi <- fread("output/2_replicating_liu_et_al/2_mr/1_expo_outcome_df/whradjbmi_liu_PCOS (adj age+BMI).txt")

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

saveRDS(mr_post, "output/2_replicating_liu_et_al/2_mr/2_mr/whradjbmi/broad/instruments_after_outlier_removal_df.RDS")

#STEP 6: let's run the results for real now:

mr_res <- TwoSampleMR::mr(mr_post)

# id.exposure id.outcome outcome exposure                    method nsnp           b         se         pval
# 1    exposure    outcome outcome exposure                  MR Egger  234  0.59445382 0.13175675 1.021817e-05
# 2    exposure    outcome outcome exposure           Weighted median  234  0.30486664 0.10551070 3.859309e-03
# 3    exposure    outcome outcome exposure Inverse variance weighted  234  0.19748897 0.06149488 1.320578e-03
# 4    exposure    outcome outcome exposure               Simple mode  234 -0.06648457 0.29472151 8.217224e-01
# 5    exposure    outcome outcome exposure             Weighted mode  234  0.37686916 0.15870129 1.837380e-02

#Let's check the plot:

rucker <- TwoSampleMR::mr_rucker(mr_post)

# [1]]$intercept
# Method    Estimate          SE      CI_low       CI_upp            P
# 1  Egger fixed effects -0.01175811 0.003188565 -0.01800758 -0.005508638 0.0002263907
# 2 Egger random effects -0.01175811 0.003188565 -0.01800758 -0.005508638 0.9998868047
# 
# [[1]]$Q
# Method         Q  df            P
# 1   Q_ivw 209.60518 233 0.8623955744
# 2 Q_egger 197.99975 232 0.9486103283
# 3  Q_diff  11.60543   1 0.0006575957
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
  
saveRDS(mr_post, "output/2_replicating_liu_et_al//2_mr/2_mr/whradjbmi/consortium/instruments_after_outlier_removal_df.RDS")
  
#STEP 6: let's run the results for real now:
  
mr_res <- TwoSampleMR::mr(mr_post)
  
# id.exposure id.outcome outcome exposure                    method nsnp          b         se         pval
# 1    exposure    outcome outcome exposure                  MR Egger  219 0.20084788 0.05017473 8.587339e-05
# 2    exposure    outcome outcome exposure           Weighted median  219 0.12512096 0.04089166 2.214727e-03
# 3    exposure    outcome outcome exposure Inverse variance weighted  219 0.08513616 0.02308083 2.254831e-04
# 4    exposure    outcome outcome exposure               Simple mode  219 0.04708880 0.08763057 5.915699e-01
# 5    exposure    outcome outcome exposure             Weighted mode  219 0.13526485 0.05043504 7.879623e-03

# #Let's check the plot:
#  
rucker <- TwoSampleMR::mr_rucker(mr_post)

# [[1]]$intercept
# Method     Estimate          SE       CI_low        CI_upp           P
# 1  Egger fixed effects -0.003397927 0.001262344 -0.005872076 -0.0009237772 0.007107631
# 2 Egger random effects -0.003397927 0.001262344 -0.005872076 -0.0009237772 0.996446184
# 
# [[1]]$Q
# Method          Q  df          P
# 1   Q_ivw 208.782505 218 0.66121057
# 2 Q_egger 202.036564 217 0.75905917
# 3  Q_diff   6.745941   1 0.00939612
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


