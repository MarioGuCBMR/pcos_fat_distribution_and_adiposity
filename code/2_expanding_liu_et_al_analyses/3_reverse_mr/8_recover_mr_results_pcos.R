##############
#INTRODUCTION#
##############

#This code recovers MR results and produces a table with them

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)
library(TwoSampleMR)

###################
#Loading functions#
###################

#source("N:/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/Mediation_analysis/code/0_functions/functions_4_mediation.R")

lower_ci_simple <- function(beta_, se_){
  #A function that takes the beta and se from summary statistics
  #and generates the lower CI.
  
  lower_ci <- beta_ - qnorm(0.975)*se_
  
  return(lower_ci)
  
}

upper_ci_simple <- function(beta_, se_){
  #A function that takes the beta and se from summary statistics
  #and generates the upper CI.
  
  upper_ci <- beta_ + qnorm(0.975)*se_
  
  return(upper_ci)
  
}

lower_ci <- function(beta_, se_, nsnp, method){
  lower_ci_end <- c()
  
  for(index in seq_along(beta_)){
    dist <- method[index]
    
    # Choose critical value based on method
    crit_val <- switch(dist,
                       "MR Egger"           = qt(0.975, df = nsnp[index] - 2),
                       "Weighted mode"      = qt(0.975, df = nsnp[index] - 1),
                       "Inverse variance weighted" = qnorm(0.975),
                       "Wald ratio" = qnorm(0.975),
                       "Weighted median"    = qnorm(0.975),
                       stop(paste("Unknown method:", dist)))
    
    lower_ci_ <- beta_[index] - crit_val * se_[index]
    lower_ci_end <- c(lower_ci_end, lower_ci_)
  }
  
  return(lower_ci_end)
}

upper_ci <- function(beta_, se_, nsnp, method){
  upper_ci_end <- c()
  
  for(index in seq_along(beta_)){
    dist <- method[index]
    
    crit_val <- switch(dist,
                       "MR Egger"           = qt(0.975, df = nsnp[index] - 2),
                       "Weighted mode"      = qt(0.975, df = nsnp[index] - 1),
                       "Inverse variance weighted" = qnorm(0.975),
                       "Wald ratio" = qnorm(0.975),
                       "Weighted median"    = qnorm(0.975),
                       stop(paste("Unknown method:", dist)))
    
    upper_ci_ <- beta_[index] + crit_val * se_[index]
    upper_ci_end <- c(upper_ci_end, upper_ci_)
  }
  
  return(upper_ci_end)
}


stats_compiler <- function(dir_, out_check){
  
  #STEP 1: do a conditional and assess:
  
  if(out_check == "outliers_removed"){
    
    mr_stats <- readRDS(paste(dir_, "/mr_res_after_outlier_extraction.RDS", sep = ""))
    
  } else {
    
    mr_stats <- readRDS(paste(dir_, "/mr_res_before_outlier_extraction.RDS", sep = ""))
    
  }
  
  return(mr_stats)
  
}

sensitivity_compiler <- function(dir_, out_check){
  
  #STEP 1: do a conditional and assess:
  
  if(out_check == "outliers_removed"){
    
    rucker <- readRDS(paste(dir_, "/sensitivity_test_after_outlier_extraction.RDS", sep = ""))
    het <- readRDS(paste(dir_, "/strict_het_test_after_outlier_extraction.RDS", sep = ""))
    
  } else {
    
    rucker <- readRDS(paste(dir_, "/sensitivity_test_before_outlier_extraction.RDS", sep = ""))
    het <- readRDS(paste(dir_, "/strict_het_test_before_outlier_extraction.RDS", sep = ""))
    
  }
  
  return(list(rucker, het))
  
}

compiling_res <- function(list_of_dirs){
  
  #STEP 0:initiate loop: 
  
  for(index in seq(1, length(list_of_dirs))){
    
    #STEP 0: get the outcome:
    
    outcome <- as.character(unlist(str_split(list_of_dirs[index], "/")))
    outcome <- outcome[length(outcome)]
    print(outcome)
    
    #STEP 1: check data
    
    list_of_files <- list.files(list_of_dirs[index])
    
    analysis_check <- ifelse("mr_res_before_outlier_extraction.RDS"%in%list_of_files, "Performed", "other")
    outlier_check <- ifelse("mr_res_after_outlier_extraction.RDS"%in%list_of_files, "outliers_removed", "standard")
    nsnps_check <- ifelse("sensitivity_test_before_outlier_extraction.RDS"%in%list_of_files, analysis_check, "nSNPs < 3")
    
    if(nsnps_check == "nSNPs < 3"){
      
      outlier_check <- "standard"
      
      mr_stats <- stats_compiler(list_of_dirs[index], outlier_check)
      mr_stats$outcome <- outcome
      mr_stats$analyses <- outlier_check
      
      mr_stats$Q_diff <- NA
      mr_stats$egger_intercept <- NA
      mr_stats$I2 <- NA
      
      if(!(exists("final_df"))){
        
        final_df <- mr_stats
        
      } else {
        
        final_df <- rbind(final_df, mr_stats)
        
      }
      
      next()
      
    }
    
    #STEP 2: get the stats!!!!
    
    if(analysis_check == "Performed"){
      
      #CAREFUL, we will start with an exception. 
      #Maybe the analysis ran with 3 SNPs, but outlier removal took 1 out.
      #We will allow it to pass, but we will check the sensitivity tests...
      
      nsnps_check <- ifelse("sensitivity_test_after_outlier_extraction.RDS"%in%list_of_files, analysis_check, "nSNPs < 3")
      
      if(nsnps_check == "nSNPs < 3"){
        
        outlier_check <- "standard"
        
        mr_stats <- stats_compiler(list_of_dirs[index], outlier_check)
        mr_stats$outcome <- outcome
        mr_stats$analyses <- outlier_check
        
        
      }
      
      #If everything is normal, we are cool:
      
      mr_stats <- stats_compiler(list_of_dirs[index], outlier_check)
      mr_stats$outcome <- outcome
      mr_stats$analyses <- outlier_check
      
      #Let's add the data if the sensitivity analyses were performed:
      
      sensitivity_list <- sensitivity_compiler(list_of_dirs[index], outlier_check)
      
      q_ <- format(sensitivity_list[[1]][[1]]$Q$Q[3], scientific = FALSE, digits = 2)
      q_p <- format(sensitivity_list[[1]][[1]]$Q$P[3], scientific = TRUE, digits = 3)
      cochran_q_diff <- paste(q_, " P=", q_p, sep = "") 
      
      #cochran_q_diff <- paste(sensitivity_list[[1]][[1]]$Q$Q[3], " P=", sensitivity_list[[1]][[1]]$Q$P[3], sep = "") #gets in the rucker list, extracts Q info for the Q_diff [3] between IVW and Egger!
      
      int_ <- format(sensitivity_list[[1]][[1]]$intercept$Estimate[1], scientific = TRUE, digits = 3)
      int_p <- format(sensitivity_list[[1]][[1]]$intercept$P[1], scientific = TRUE, digits = 3)
      
      #intercept <- paste(sensitivity_list[[1]][[1]]$intercept$Estimate[1], " P=", sensitivity_list[[1]][[1]]$intercept$P[1], sep = "")
      intercept <- paste(int_, " P=", int_p, sep = "") 
      
      
      isq <- unlist(sensitivity_list[[2]])
      
      if(length(isq) != 1){ #this only happens when it is empty
        
        isq <- paste(format(isq$I2, scientifiC=FALSE, digits=2), 
                     " [", format(isq$lower.I2, scientifiC=FALSE, digits=3), 
                     ",", format(isq$upper.I2, scientifiC=FALSE, digits=2), "]", sep = "")
        
      }
      
      mr_stats$Q_diff <- cochran_q_diff
      mr_stats$egger_intercept <- intercept
      mr_stats$I2 <- isq
      
      
    } else {
      
      mr_stats <- c("id.exposure", "id.outcome", "outcome", "exposure", "method", "nsnp", "b", "se", "pval")
      mr_stats <- as.data.frame(rbind(mr_stats, rep(NA, 9)))
      colnames(mr_stats) <- mr_stats[1,]
      mr_stats <- mr_stats[-1,]
      
      mr_stats$outcome <- outcome
      mr_stats$analyses <- analysis_check
      
      mr_stats$Q_diff <- NA
      mr_stats$egger_intercept <- NA
      mr_stats$I2 <- NA
      
    }
    
    #STEP 3: let's check now the sensitivity tests:
    
    if(!(exists("final_df"))){
      
      final_df <- mr_stats
      
    } else {
      
      final_df <- rbind(final_df, mr_stats)
      
    }
    
  }
  
  return(final_df)
  
}

mr_plotter_first <- function(mr_df, title){
  #This plot takes the dataframe of result from MR and tries to plot it so that it looks amazing and beautiful
  
  #STEP 1: make dummy data
  
  betas_rounded <- format(mr_df$b, scientific = FALSE, digits = 4)
  lower_ci_rounded <- format(mr_df$lower_ci,scientific = FALSE, digits = 4)
  upper_ci_rounded <- format(mr_df$upper_ci, scientific = FALSE, digits = 4)
  
  mr_df$betas <- betas_rounded
  mr_df$lower_ci <- lower_ci_rounded
  mr_df$upper_ci <- upper_ci_rounded
  
  #STEP 2: change easy_for_code handles to beauty_informative_for_figure:
  
  mr_df$outcome_=NA
  
  mr_df$outcome_ <- ifelse(mr_df$outcome == "doctor", "PCOS (FinnGen)", mr_df$outcome_)
  mr_df$outcome_ <- ifelse(mr_df$outcome == "broad", "PCOS (FinnGen - Broad)", mr_df$outcome_)
  mr_df$outcome_ <- ifelse(mr_df$outcome == "consortium", "PCOS (FinnGen - Consortium)", mr_df$outcome_)
  mr_df$outcome_ <- ifelse(mr_df$outcome == "meta_analysis", "PCOS (Day et al)", mr_df$outcome_)
  mr_df$outcome_ <- ifelse(mr_df$outcome == "venkatesh", "PCOS (Venkatesh et al)", mr_df$outcome_)
  mr_df$outcome_ <- ifelse(mr_df$outcome == "age_adjusted", "PCOS (Tyrmi et al)", mr_df$outcome_)
  mr_df$outcome_ <- ifelse(mr_df$outcome == "age_and_bmi_adjusted", "PCOSadjBMI (Tyrmi et al)", mr_df$outcome_)
  mr_df$outcome_ <- ifelse(mr_df$outcome == "amenorrhea", "Amenorrhea", mr_df$outcome_)
  mr_df$outcome_ <- ifelse(mr_df$outcome == "oligomenorrhea", "Oligomenorrhea", mr_df$outcome_)
  mr_df$outcome_ <- ifelse(mr_df$outcome == "menorrhagia", "Menorrhagia", mr_df$outcome_)
  mr_df$outcome_ <- ifelse(mr_df$outcome == "anovulation", "Anovulation", mr_df$outcome_)
  mr_df$outcome_ <- ifelse(mr_df$outcome == "hirsutism", "Hirsutism", mr_df$outcome_)
  
  mr_df$outcome=mr_df$outcome
  
  #Let's add the type to make it more beautiful:
  
  clinically_diagnosed=c("PCOS (FinnGen)", "PCOS (Tyrmi et al)", "PCOSadjBMI (Tyrmi et al)")
  extended_criteria=c("PCOS (FinnGen - Broad)", "PCOS (FinnGen - Consortium)")
  meta_analysis=c("PCOS (Day et al)", "PCOS (Venkatesh et al)")
  female_infertiltiy=c("Amenorrhea", "Oligomenorrhea", "Menorrhagia", "Anovulation", "Hirsutism")
  
  mr_df$type=NA
  mr_df$type=ifelse(mr_df$outcome_%in%clinically_diagnosed, "Clinically diagnosed", mr_df$type)
  mr_df$type=ifelse(mr_df$outcome_%in%extended_criteria, "Extended criteria", mr_df$type)
  mr_df$type=ifelse(mr_df$outcome_%in%meta_analysis, "Meta-analysis (with self-reported)", mr_df$type)
  mr_df$type=ifelse(mr_df$outcome_%in%female_infertiltiy, "Female infertility", mr_df$type)
  
  #Let's convert the data into factors so that we can organize the results:
  
  mr_df$outcome_ <- factor(mr_df$outcome_, levels = c("PCOS (FinnGen)", "PCOS (FinnGen - Broad)", "PCOS (FinnGen - Consortium)",  "PCOS (Day et al)", "PCOS (Venkatesh et al)", "PCOS (Tyrmi et al)", "PCOSadjBMI (Tyrmi et al)", "Amenorrhea", "Oligomenorrhea", "Menorrhagia", "Anovulation", "Hirsutism"))
  mr_df$method  <- factor(mr_df$method, levels = c("Weighted mode", "Weighted median", "Wald ratio", "MR Egger", "Inverse variance weighted", "Genetic correlation"))
  mr_df$type <- factor(mr_df$type, levels = c("Clinically diagnosed", "Extended criteria", "Meta-analysis (with self-reported)",  "Female infertility"))
  
  # Ensure pvals numeric and format
  mr_df$pvals_num <- as.numeric(mr_df$pval)
  mr_df$pval_formatted <- paste0("p=", format(mr_df$pvals_num, scientific = TRUE, digits = 2))
  
  # Create significance stars column based on p-value thresholds
  mr_df$signif_stars <- ""
  mr_df$signif_stars[mr_df$pvals_num < 0.05] <- "*"
  mr_df$signif_stars[mr_df$pvals_num < 5e-8] <- "**"
  
  # Position to place stars slightly above upper_ci
  star_y_pos <- as.numeric(mr_df$upper_ci) + 0.2  # tweak this offset if needed
  
  #Now we can make our plot!
  
  plotio <- ggplot(mr_df, aes(x = fct_rev(outcome_), y = as.numeric(betas))) + #make a plot that takes into account outcome and betas
    
    geom_point(aes(color=method, shape=method), 
               size = 7, position = position_dodge(width = 0.75)) + #color by method, make them a bit separated to distinguish
    
    #A little change so that we can have a more beatiful plot:
    
    scale_shape_manual(values = c(
      "Genetic correlation" = 23,
      "Inverse variance weighted" = 16,
      "Wald ratio" = 16,
      "Weighted median" = 16,
      "Weighted mode" = 16,
      "MR Egger" = 16
    )) +
    
    geom_errorbar(
      aes(
        ymin = as.numeric(lower_ci),
        ymax = as.numeric(upper_ci),
        group = method  # this is what was missing
      ),
      width = 0.25,
      position = position_dodge(width = 0.75),
      color = "black",
      linewidth = 1
    ) +
    
    geom_text(aes(y = star_y_pos, label = signif_stars, group = method),
              position = position_dodge(width = 0.75),
              size = 6, color = "black", fontface = "bold")+
    
    facet_wrap(~type, strip.position = "left", nrow = 7, scales = "free_y") + #separate by outcome
    
    geom_hline(yintercept = 0, color = "red", size = 0.5, linetype = 2) +
    
    coord_flip() +
    xlab("Trait") +
    ylab("MR causal effect sizes") +
    ggtitle(title) +
    theme_bw() +
    theme(
      legend.position = "none",
      legend.text = element_text(size = 11),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(),
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 11.5),
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5)
    ) +
    guides(
      color = guide_legend(reverse = TRUE),
      shape = guide_legend(reverse = TRUE)
    )
  
  plotio <- plotio + scale_y_continuous(limits = c(-2, 2))
  
  return(plotio)
  
}

mr_plotter_mid <- function(mr_df, title){
  #This plot takes the dataframe of result from MR and tries to plot it so that it looks amazing and beautiful
  
  #STEP 1: make dummy data
  
  betas_rounded <- format(mr_df$b, scientific = FALSE, digits = 4)
  lower_ci_rounded <- format(mr_df$lower_ci,scientific = FALSE, digits = 4)
  upper_ci_rounded <- format(mr_df$upper_ci, scientific = FALSE, digits = 4)
  
  mr_df$betas <- betas_rounded
  mr_df$lower_ci <- lower_ci_rounded
  mr_df$upper_ci <- upper_ci_rounded
  
  #STEP 2: change easy_for_code handles to beauty_informative_for_figure:
  
  mr_df$outcome_=NA
  
  mr_df$outcome_ <- ifelse(mr_df$outcome == "doctor", "PCOS (FinnGen)", mr_df$outcome_)
  mr_df$outcome_ <- ifelse(mr_df$outcome == "broad", "PCOS (FinnGen - Broad)", mr_df$outcome_)
  mr_df$outcome_ <- ifelse(mr_df$outcome == "consortium", "PCOS (FinnGen - Consortium)", mr_df$outcome_)
  mr_df$outcome_ <- ifelse(mr_df$outcome == "meta_analysis", "PCOS (Day et al)", mr_df$outcome_)
  mr_df$outcome_ <- ifelse(mr_df$outcome == "venkatesh", "PCOS (Venkatesh et al)", mr_df$outcome_)
  mr_df$outcome_ <- ifelse(mr_df$outcome == "age_adjusted", "PCOS (Tyrmi et al)", mr_df$outcome_)
  mr_df$outcome_ <- ifelse(mr_df$outcome == "age_and_bmi_adjusted", "PCOSadjBMI (Tyrmi et al)", mr_df$outcome_)
  mr_df$outcome_ <- ifelse(mr_df$outcome == "amenorrhea", "Amenorrhea", mr_df$outcome_)
  mr_df$outcome_ <- ifelse(mr_df$outcome == "oligomenorrhea", "Oligomenorrhea", mr_df$outcome_)
  mr_df$outcome_ <- ifelse(mr_df$outcome == "menorrhagia", "Menorrhagia", mr_df$outcome_)
  mr_df$outcome_ <- ifelse(mr_df$outcome == "anovulation", "Anovulation", mr_df$outcome_)
  mr_df$outcome_ <- ifelse(mr_df$outcome == "hirsutism", "Hirsutism", mr_df$outcome_)
  
  mr_df$outcome=mr_df$outcome
  
  #Let's add the type to make it more beautiful:
  
  clinically_diagnosed=c("PCOS (FinnGen)", "PCOS (Tyrmi et al)", "PCOSadjBMI (Tyrmi et al)")
  extended_criteria=c("PCOS (FinnGen - Broad)", "PCOS (FinnGen - Consortium)")
  meta_analysis=c("PCOS (Day et al)", "PCOS (Venkatesh et al)")
  female_infertiltiy=c("Amenorrhea", "Oligomenorrhea", "Menorrhagia", "Anovulation", "Hirsutism")
  
  mr_df$type=NA
  mr_df$type=ifelse(mr_df$outcome_%in%clinically_diagnosed, "Clinically diagnosed", mr_df$type)
  mr_df$type=ifelse(mr_df$outcome_%in%extended_criteria, "Extended criteria", mr_df$type)
  mr_df$type=ifelse(mr_df$outcome_%in%meta_analysis, "Meta-analysis (with self-reported)", mr_df$type)
  mr_df$type=ifelse(mr_df$outcome_%in%female_infertiltiy, "Female infertility", mr_df$type)
  
  #Let's convert the data into factors so that we can organize the results:
  
  mr_df$outcome_ <- factor(mr_df$outcome_, levels = c("PCOS (FinnGen)", "PCOS (FinnGen - Broad)", "PCOS (FinnGen - Consortium)",  "PCOS (Day et al)", "PCOS (Venkatesh et al)", "PCOS (Tyrmi et al)", "PCOSadjBMI (Tyrmi et al)", "Amenorrhea", "Oligomenorrhea", "Menorrhagia", "Anovulation", "Hirsutism"))
  mr_df$method  <- factor(mr_df$method, levels = c("Weighted mode", "Weighted median", "Wald ratio", "MR Egger", "Inverse variance weighted", "Genetic correlation"))
  mr_df$type <- factor(mr_df$type, levels = c("Clinically diagnosed", "Extended criteria", "Meta-analysis (with self-reported)",  "Female infertility"))
  
  # Ensure pvals numeric and format
  mr_df$pvals_num <- as.numeric(mr_df$pval)
  mr_df$pval_formatted <- paste0("p=", format(mr_df$pvals_num, scientific = TRUE, digits = 2))
  
  # Create significance stars column based on p-value thresholds
  mr_df$signif_stars <- ""
  mr_df$signif_stars[mr_df$pvals_num < 0.05] <- "*"
  mr_df$signif_stars[mr_df$pvals_num < 5e-8] <- "**"
  
  # Position to place stars slightly above upper_ci
  star_y_pos <- as.numeric(mr_df$upper_ci) + 0.2  # tweak this offset if needed
  
  #Now we can make our plot!
  
  plotio <- ggplot(mr_df, aes(x = fct_rev(outcome_), y = as.numeric(betas))) + #make a plot that takes into account outcome and betas
    
    geom_point(aes(color=method, shape=method), 
               size = 7, position = position_dodge(width = 0.75)) + #color by method, make them a bit separated to distinguish
    
    #A little change so that we can have a more beatiful plot:
    
    scale_shape_manual(values = c(
      "Genetic correlation" = 23,
      "Inverse variance weighted" = 16,
      "Wald ratio" = 16,
      "Weighted median" = 16,
      "Weighted mode" = 16,
      "MR Egger" = 16
    )) +
    
    geom_errorbar(
      aes(
        ymin = as.numeric(lower_ci),
        ymax = as.numeric(upper_ci),
        group = method  # this is what was missing
      ),
      width = 0.25,
      position = position_dodge(width = 0.75),
      color = "black",
      linewidth = 1
    ) +
    
    geom_text(aes(y = star_y_pos, label = signif_stars, group = method),
              position = position_dodge(width = 0.75),
              size = 6, color = "black", fontface = "bold")+
    
    facet_wrap(~type, strip.position = "left", nrow = 7, scales = "free_y") +
    
    geom_hline(yintercept = 0, color = "red", size = 0.5, linetype = 2) +
    
    coord_flip() +
    xlab("Trait") +
    ylab("MR causal effect sizes") +
    ggtitle(title) +
    theme_bw() +
    theme(
      legend.position = "none",
      legend.text = element_text(size = 11),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      strip.text.y.left = element_blank(), 
      axis.line = element_line(),
      axis.text.y = element_blank(),
      axis.ticks.y=element_blank(),
      axis.title.x = element_text(size = 11.5),
      axis.title.y = element_blank(),
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5)
    ) +
    guides(
      color = guide_legend(reverse = TRUE),
      shape = guide_legend(reverse = TRUE)
    )
  
  plotio <- plotio + scale_y_continuous(limits = c(-2, 2))
  
  return(plotio)
  
}

mr_plotter_last <- function(mr_df, title){
  #This plot takes the dataframe of result from MR and tries to plot it so that it looks amazing and beautiful
  
  #STEP 1: make dummy data
  
  betas_rounded <- format(mr_df$b, scientific = FALSE, digits = 4)
  lower_ci_rounded <- format(mr_df$lower_ci,scientific = FALSE, digits = 4)
  upper_ci_rounded <- format(mr_df$upper_ci, scientific = FALSE, digits = 4)
  
  mr_df$betas <- betas_rounded
  mr_df$lower_ci <- lower_ci_rounded
  mr_df$upper_ci <- upper_ci_rounded
  
  #STEP 2: change easy_for_code handles to beauty_informative_for_figure:
  
  mr_df$outcome_=NA
  
  mr_df$outcome_ <- ifelse(mr_df$outcome == "doctor", "PCOS (FinnGen)", mr_df$outcome_)
  mr_df$outcome_ <- ifelse(mr_df$outcome == "broad", "PCOS (FinnGen - Broad)", mr_df$outcome_)
  mr_df$outcome_ <- ifelse(mr_df$outcome == "consortium", "PCOS (FinnGen - Consortium)", mr_df$outcome_)
  mr_df$outcome_ <- ifelse(mr_df$outcome == "meta_analysis", "PCOS (Day et al)", mr_df$outcome_)
  mr_df$outcome_ <- ifelse(mr_df$outcome == "venkatesh", "PCOS (Venkatesh et al)", mr_df$outcome_)
  mr_df$outcome_ <- ifelse(mr_df$outcome == "age_adjusted", "PCOS (Tyrmi et al)", mr_df$outcome_)
  mr_df$outcome_ <- ifelse(mr_df$outcome == "age_and_bmi_adjusted", "PCOSadjBMI (Tyrmi et al)", mr_df$outcome_)
  mr_df$outcome_ <- ifelse(mr_df$outcome == "amenorrhea", "Amenorrhea", mr_df$outcome_)
  mr_df$outcome_ <- ifelse(mr_df$outcome == "oligomenorrhea", "Oligomenorrhea", mr_df$outcome_)
  mr_df$outcome_ <- ifelse(mr_df$outcome == "menorrhagia", "Menorrhagia", mr_df$outcome_)
  mr_df$outcome_ <- ifelse(mr_df$outcome == "anovulation", "Anovulation", mr_df$outcome_)
  mr_df$outcome_ <- ifelse(mr_df$outcome == "hirsutism", "Hirsutism", mr_df$outcome_)
  
  mr_df$outcome=mr_df$outcome
  
  #Let's add the type to make it more beautiful:
  
  clinically_diagnosed=c("PCOS (FinnGen)", "PCOS (Tyrmi et al)", "PCOSadjBMI (Tyrmi et al)")
  extended_criteria=c("PCOS (FinnGen - Broad)", "PCOS (FinnGen - Consortium)")
  meta_analysis=c("PCOS (Day et al)", "PCOS (Venkatesh et al)")
  female_infertiltiy=c("Amenorrhea", "Oligomenorrhea", "Menorrhagia", "Anovulation", "Hirsutism")
  
  mr_df$type=NA
  mr_df$type=ifelse(mr_df$outcome_%in%clinically_diagnosed, "Clinically diagnosed", mr_df$type)
  mr_df$type=ifelse(mr_df$outcome_%in%extended_criteria, "Extended criteria", mr_df$type)
  mr_df$type=ifelse(mr_df$outcome_%in%meta_analysis, "Meta-analysis (with self-reported)", mr_df$type)
  mr_df$type=ifelse(mr_df$outcome_%in%female_infertiltiy, "Female infertility", mr_df$type)
  
  #Let's convert the data into factors so that we can organize the results:
  
  mr_df$outcome_ <- factor(mr_df$outcome_, levels = c("PCOS (FinnGen)", "PCOS (FinnGen - Broad)", "PCOS (FinnGen - Consortium)",  "PCOS (Day et al)", "PCOS (Venkatesh et al)", "PCOS (Tyrmi et al)", "PCOSadjBMI (Tyrmi et al)", "Amenorrhea", "Oligomenorrhea", "Menorrhagia", "Anovulation", "Hirsutism"))
  mr_df$method  <- factor(mr_df$method, levels = c("Weighted mode", "Weighted median", "Wald ratio", "MR Egger", "Inverse variance weighted", "Genetic correlation"))
  mr_df$type <- factor(mr_df$type, levels = c("Clinically diagnosed", "Extended criteria", "Meta-analysis (with self-reported)",  "Female infertility"))
  
  # Ensure pvals numeric and format
  mr_df$pvals_num <- as.numeric(mr_df$pval)
  mr_df$pval_formatted <- paste0("p=", format(mr_df$pvals_num, scientific = TRUE, digits = 2))
  
  # Create significance stars column based on p-value thresholds
  mr_df$signif_stars <- ""
  mr_df$signif_stars[mr_df$pvals_num < 0.05] <- "*"
  mr_df$signif_stars[mr_df$pvals_num < 5e-8] <- "**"
  
  # Position to place stars slightly above upper_ci
  star_y_pos <- as.numeric(mr_df$upper_ci) + 0.2  # tweak this offset if needed
  
  #Now we can make our plot!
  
  plotio <- ggplot(mr_df, aes(x = fct_rev(outcome_), y = as.numeric(betas))) + #make a plot that takes into account outcome and betas
    
    geom_point(aes(color=method, shape=method), 
               size = 7, position = position_dodge(width = 0.75)) + #color by method, make them a bit separated to distinguish
    
    #A little change so that we can have a more beatiful plot:
    
    scale_shape_manual(values = c(
      "Genetic correlation" = 23,
      "Inverse variance weighted" = 16,
      "Wald ratio" = 16,
      "Weighted median" = 16,
      "Weighted mode" = 16,
      "MR Egger" = 16
    )) +
    
    geom_errorbar(
      aes(
        ymin = as.numeric(lower_ci),
        ymax = as.numeric(upper_ci),
        group = method  # this is what was missing
      ),
      width = 0.25,
      position = position_dodge(width = 0.75),
      color = "black",
      linewidth = 1
    ) +
    
    geom_text(aes(y = star_y_pos, label = signif_stars, group = method),
              position = position_dodge(width = 0.75),
              size = 6, color = "black", fontface = "bold")+
    
    facet_wrap(~type, strip.position = "left", nrow = 7, scales = "free_y") +
    
    geom_hline(yintercept = 0, color = "red", size = 0.5, linetype = 2) +
    
    coord_flip() +
    xlab("Trait") +
    ylab("MR causal effect sizes") +
    ggtitle(title) +
    theme_bw() +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 11),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      strip.text.y.left = element_blank(), 
      axis.line = element_line(),
      axis.text.y = element_blank(),
      axis.ticks.y=element_blank(),
      axis.title.x = element_text(size = 11.5),
      axis.title.y = element_blank(),
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5)
    ) +
    guides(
      color = guide_legend(reverse = TRUE),
      shape = guide_legend(reverse = TRUE)
    )
  
  plotio <- plotio + scale_y_continuous(limits = c(-2, 2))
  
  return(plotio)
  
}


compute_p_value <- function(z) {
  return(2 * (1 - pnorm(abs(z))))
}

result_wrapper <- function(path_2_res){
  #Cleans a systematically retrieves results for a certain path
  
  outcome_list <- list.dirs(path_2_res)
  outcome_list <- outcome_list[which(!(outcome_list == path_2_res))]
  
  res <- compiling_res(outcome_list)
  
  #Let's remove simple mode:
  
  res <- res[which(res$method != "Simple mode"),]
  #res <- res[which(res$method != "MR Egger"),]
  
  
  #Let's return the data:
  
  return(res)
  
}

results_updated <- function(res){
  
  #Let's add the confidence intervals:
  
  res$lower_ci_ <- lower_ci(res$b, res$se, res$nsnp, res$method)
  res$upper_ci_ <- upper_ci(res$b, res$se, res$nsnp, res$method)
  
  betas_rounded <- round(as.numeric(res$b), 2)
  lower_ci_rounded <- round(as.numeric(res$lower_ci_), digits = 2)
  upper_ci_rounded <- round(as.numeric(res$upper_ci_), digits = 2)
  
  res$full_effect <- paste(betas_rounded,
                                     " (",
                                     lower_ci_rounded,
                                     ",",
                                     upper_ci_rounded,
                                     ")",
                                     sep=""
  )
  
  # Define desired order for outcome and method
  desired_outcome_order <- c("doctor", "broad", "consortium", "meta_analysis", "venkatesh", "age_adjusted", "age_and_bmi_adjusted", "amenorrhea", "oligomenorrhea", "menorrhagia", "anovulation","hirsutism")
  desired_method_order <- c("Inverse variance weighted", "MR Egger", "Weighted median", "Weighted mode", "Wald ratio")
  
  # Apply ordering only if both columns exist
  if (all(c("outcome", "method") %in% names(res))) {
    res$outcome <- factor(res$outcome, levels = desired_outcome_order)
    res$method <- factor(res$method, levels = desired_method_order)
    res <- res[order(res$outcome, res$method), ]
  }
  
  # Reorder columns: bring selected ones to the front
  key_cols <- c("method", "nsnp", "full_effect", "pval", "Q_diff", "egger_intercept", "I2")
  other_cols <- setdiff(names(res), key_cols)
  res <- res[, c(key_cols, other_cols)]
  
  return(res)
  
}

gc_2_mr_df <- function(doctor = c(NA, NA), 
                       broad = c(NA, NA),
                       consortium = c(NA, NA),
                       meta = c(NA, NA), 
                       venkatesh = c(NA, NA), 
                       adj_age = c(NA, NA), 
                       adj_age_bmi = c(NA, NA),
                       amenorrhea = c(NA, NA),
                       oligomenorrhea = c(NA, NA),
                       menorrhagia = c(NA, NA),
                       anovulation = c(NA, NA),
                       hirsutism = c(NA, NA)){
  
  #This function takes GC and SEs from log file and returns a well-structured dataframe!
  
  gc <- rbind(as.data.frame(t(doctor)), as.data.frame(t(broad)), as.data.frame(t(consortium)), as.data.frame(t(meta)), as.data.frame(t(venkatesh)), as.data.frame(t(adj_age)),
                   as.data.frame(t(adj_age_bmi)), as.data.frame(t(amenorrhea)), as.data.frame(t(oligomenorrhea)), as.data.frame(t(menorrhagia)), as.data.frame(t(anovulation)), as.data.frame(t(hirsutism)))
  
  colnames(gc) <- c("b", "se")
  
  gc$outcome <- c("doctor", "broad", "consortium", "meta_analysis", "venkatesh", "age_adjusted", "age_and_bmi_adjusted", "amenorrhea", "oligomenorrhea", "menorrhagia", "anovulation", "hirsutism")
  gc$method <- rep("Genetic correlation", 12)
  gc$pval <- compute_p_value(gc$b/gc$se)
  
  #Let's add dummy columns so that we can fit in the results:
  
  gc$exposure <- "gc"
  gc$id.exposure <- gc$exposure 
  gc$id.outcome <- gc$outcome
  
  gc$nsnp <- NA
  gc$Q_diff <- NA
  gc$egger_intercept <- NA
  gc$I2 <- NA
  gc$analyses <- NA
  
  #Let's also add the CIs:
  
  gc$lower_ci_ <- lower_ci_simple(gc$b, gc$se)
  gc$upper_ci_ <- upper_ci_simple(gc$b, gc$se)
  
  betas_rounded <- format(as.numeric(gc$b), scientific = FALSE, digits = 2)
  lower_ci_rounded <- format(gc$lower_ci,scientific = FALSE, digits = 2)
  upper_ci_rounded <- format(gc$upper_ci, scientific = FALSE, digits = 2)
  
  gc$full_effect <- paste(betas_rounded,
                                     " (",
                                     lower_ci_rounded,
                                     ",",
                                     upper_ci_rounded,
                                     ")",
                                     sep=""
  )
  
  gc$pval <- format(gc$pval, scientific = TRUE, digits = 3)
  
  
  #And return!
  
  return(gc)
  
}

##############
#Loading data#
##############

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/PCOS_2026/")

bmi_res <- "output/2_replicating_liu_et_al/3_reverse_mr/2_mr/bmi/"
whr_res <- "output/2_replicating_liu_et_al/3_reverse_mr/2_mr/whr/"
whradjbmi_res <- "output/2_replicating_liu_et_al/3_reverse_mr/2_mr/whradjbmi/"

##############################################################################################################################
#The approach will be making a function that loops over the outcomes in each folder that we give and searches for the results#
##############################################################################################################################

bmi_df <- result_wrapper(bmi_res)
whr_df <- result_wrapper(whr_res)
whradjbmi_df <- result_wrapper(whradjbmi_res)

########################################################################
#We need to do a couple of adjustments, but aside from that we are good#
########################################################################

bmi_df <- results_updated(bmi_df)
whr_df <- results_updated(whr_df)
whradjbmi_df <- results_updated(whradjbmi_df)

##########################################################################################
#Let's plot PCOS results and additional traits later, because the plot is a bit too messy# BMI first and see what is up#
##########################################################################################

bmi_pcos_plot=mr_plotter_first(bmi_df, "")
whr_pcos_plot=mr_plotter_mid(whr_df, "")
whradjbmi_pcos_plot=mr_plotter_last(whradjbmi_df, "")

library(patchwork)

pcos_all=bmi_pcos_plot+whr_pcos_plot+whradjbmi_pcos_plot+plot_layout(ncol=3)

dir.create("manuscript/figures/drafts")

ggsave(
  filename = "manuscript/figures/drafts/reverse_mr_bmi_whr_whradjbmi_pcos.svg",
  plot = pcos_all,
  width = 15,
  height = 11,
  units = "in",  # or "in" depending on your preference, but cm is common
  device = "svg"
)

#####################################################################
#Let's format the data so that it easier to integrate in our results#
#####################################################################

bmi_df <- bmi_df %>%
  dplyr::select(outcome, method, nsnp, "full_effect",    "se",  "pval" ,           "Q_diff",          "egger_intercept", "I2")

bmi_df$se <- format(as.numeric(bmi_df$se), scientific = TRUE, digits = 3)
bmi_df$pval <- format(as.numeric(bmi_df$pval), scientific = TRUE, digits = 3)

whr_df <- whr_df %>%
  dplyr::select(outcome, method, nsnp, "full_effect",    "se",  "pval" ,           "Q_diff",          "egger_intercept", "I2")

whr_df$se <- format(as.numeric(whr_df$se), scientific = TRUE, digits = 3)
whr_df$pval <- format(as.numeric(whr_df$pval), scientific = TRUE, digits = 3)

whradjbmi_df <- whradjbmi_df %>%
  dplyr::select(outcome, method, nsnp, "full_effect",    "se",  "pval" ,           "Q_diff",          "egger_intercept", "I2")

whradjbmi_df$se <- format(as.numeric(whradjbmi_df$se), scientific = TRUE, digits = 3)
whradjbmi_df$pval <- format(as.numeric(whradjbmi_df$pval), scientific = TRUE, digits = 3)

#Let's save the data:

dir.create("manuscript/supplementary_tables/")
dir.create("manuscript/supplementary_tables/drafts")

fwrite(bmi_df, "manuscript/supplementary_tables/drafts/supplementary_table_4_pcos_bmi.csv")
fwrite(whr_df, "manuscript/supplementary_tables/drafts/supplementary_table_4_pcos_whr.csv")
fwrite(whradjbmi_df, "manuscript/supplementary_tables/drafts/supplementary_table_4_pcos_whradjbmi.csv")

##########################################
#OK - let's do it the other way around!!!#
##########################################

# pcos_res = "output/2_replicating_liu_et_al/5_reverse_mr/2_mr/"
# 
# ##############################################################################################################################
# #The approach will be making a function that loops over the outcomes in each folder that we give and searches for the results#
# ##############################################################################################################################
# 
# pcos_df <- result_wrapper(pcos_res)
# 
# ########################################################################
# #We need to do a couple of adjustments, but aside from that we are good#
# ########################################################################
# 
# pcos_df <- results_updated(pcos_df)
# 
# #Let's format the data so that it easier to integrate in our results:
# 
# pcos_df <- pcos_df %>%
#   dplyr::select(outcome, method, nsnp, "full_effect",    "se",  "pval" ,           "Q_diff",          "egger_intercept", "I2")
# 
# pcos_df$se <- format(as.numeric(pcos_df$se), scientific = TRUE, digits = 3)
# pcos_df$pval <- format(as.numeric(pcos_df$pval), scientific = TRUE, digits = 3)
# 
# #Let's save the data:
# 
# fwrite(pcos_df, "manuscript_thesis/supplementary_tables/drafts/supplementary_table_2_pcos_whradjbmi.csv")
# 
# 
