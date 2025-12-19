##############
#INTRODUCTION#
##############

#This code performs genetic correlations

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)
library(GenomicSEM)

###################
#Loading functions#
###################

curated_2_munging_mv <- function(curated_df){
  
  #STEP 0: let's remove those that have not rsID, cuz LDSC gets confused:
  
  curated_df <- curated_df[which(str_detect(as.character(curated_df$variant), "rs") == TRUE)]
  
  #STEP 1: get hte MAF:
  
  curated_df$MAF <- curated_df$minimum_allele_frequency
  
  #STEP 2: get the right columns:
  
  curated_df_4_ldsc <- curated_df %>%
    select(variant, effect_allele, other_allele, beta, standard_error, p_value, MAF, sample_size)
  
  colnames(curated_df_4_ldsc) <- c("SNP", "A1", "A2", "BETA", "SE", "P", "MAF", "N")
  
  #STEP 3: Save the data:
  
  return(curated_df_4_ldsc)
  
}

curated_2_munging <- function(curated_df){
  
  #STEP 0: let's remove those that have not rsID, cuz LDSC gets confused:
  
  curated_df <- curated_df[which(str_detect(as.character(curated_df$variant), "rs") == TRUE)]
  
  #STEP 1: get hte MAF:
  
  curated_df$MAF <- ifelse(as.numeric(curated_df$effect_allele_frequency) > 0.50, 1-as.numeric(curated_df$effect_allele_frequency), as.numeric(curated_df$effect_allele_frequency))
  
  #STEP 2: get the right columns:
  
  curated_df_4_ldsc <- curated_df %>%
    select(variant, effect_allele, other_allele, beta, standard_error, p_value, MAF, sample_size)
  
  colnames(curated_df_4_ldsc) <- c("SNP", "A1", "A2", "BETA", "SE", "P", "MAF", "N")
  
  #STEP 3: Save the data:
  
  return(curated_df_4_ldsc)
  
}

curated_2_binary <- function(curated_df){
  
  #STEP 0: let's remove those that have not rsID, cuz LDSC gets confused:
  
  curated_df <- curated_df[which(str_detect(as.character(curated_df$variant), "rs") == TRUE)]
  
  #STEP 1: get hte MAF:
  
  curated_df$MAF <- ifelse(as.numeric(curated_df$effect_allele_frequency) > 0.50, 1-as.numeric(curated_df$effect_allele_frequency), as.numeric(curated_df$effect_allele_frequency))
  
  #STEP 2: get the effect sample size:
  
  curated_df$Neff<-4/((2*curated_df$MAF*(1-curated_df$MAF))*curated_df$standard_error^2)
  
  #STEP 2: get the right columns:
  
  curated_df_4_ldsc <- curated_df %>%
    select(variant, effect_allele, other_allele, beta, standard_error, p_value, MAF, Neff)
  
  colnames(curated_df_4_ldsc) <- c("SNP", "A1", "A2", "BETA", "SE", "P", "MAF", "N")
  
  #STEP 3: Save the data:
  
  return(curated_df_4_ldsc)
  
}


#####################################
#Loading data for dichotomous traits#
#####################################

#Running LDSC:

path_2_input <- "/projects/kilpelainen-AUDIT/people/zlc436/PCOS_2026/"

setwd(path_2_input)

LDSCoutput_SUD <- readRDS("output/2_replicating_liu_et_al/1_genetic_correlation/3_gc/genetic_correlation_res.RDS")
correlationSUD=round(LDSCoutput_SUD$S_Stand, 3) 
        
########################################
# Check which GWAS had overinflation!! #
########################################

rownames(correlationSUD)=colnames(correlationSUD)

str(LDSCoutput_SUD)
intercept=LDSCoutput_SUD$I
colnames(intercept)=colnames(LDSCoutput_SUD$S)
rownames(intercept)=colnames(LDSCoutput_SUD$S)
checkGC=data.frame(pheno = colnames(intercept),
                   intercept = diag(intercept)) # get the LDSC intercept per sumstats

checkGC=subset(checkGC, intercept >=1 ) #some had overinflation!!

#######################
#Let's clean the data:#
#######################

dir.create("output/2_replicating_liu_et_al/1_genetic_correlation/1_data_for_munging_clean/")
dir.create("output/2_replicating_liu_et_al/1_genetic_correlation/2_data_munged_cleaned")
#dir.create("output/2_replicating_liu_et_al/1_genetic_correlation/2_data_munged_cleaned_")
dir.create("output/2_replicating_liu_et_al/1_genetic_correlation/2_data_munged_cleaned_for_correction_together/")

setwd("output/2_replicating_liu_et_al/1_genetic_correlation/2_data_munged_cleaned_for_correction_together")

if(is_empty(checkGC$pheno)){
  
  print("ALL GOOD")
  
} else {
  
  for(i in 1:length(checkGC$pheno)) {
    
    print(paste0("Read in data for ", checkGC$pheno[i]))
    phenoDirmunge=paste0("../2_data_munged_for_ldsc_together/", checkGC$pheno[i], ".sumstats.gz") 
    phenoDir=paste0("../1_data_for_munging/", checkGC$pheno[i], "_4_munging.txt") 
    data <- fread(phenoDir,header=T,data.table=F)
    
    N <- fread(phenoDirmunge,header=T,data.table=F)$N[1]
    print(paste0("Multiply SE by LDSC intercept of ", checkGC$intercept[i]))
    data$SE = data$SE * sqrt(checkGC$intercept[i])
    data$Zscore=as.numeric( data$BETA)/as.numeric( data$SE) # estimate adjusted p-values
    data$P=2*pnorm(-abs(data$Zscore))
    data$Zscore=NULL
    head(data)
    
    print("Save file on cluster")
    fwrite(data, 
           file=paste0("../1_data_for_munging_clean/", checkGC$pheno[i], "_clean"), 
           sep="\t", 
           row.names = FALSE, 
           col.names = TRUE, 
           quote=F) 
    print("Remove non-GC controlled munged file and munge again")
    
    munge(files = paste0("../1_data_for_munging_clean/", checkGC$pheno[i], "_clean"), 
          hm3 = paste0("../../../../raw_data/ldsc-master/eur_w_ld_chr.tar/eur_w_ld_chr/w_hm3.snplist"),
          trait.names=paste(checkGC$pheno[i], "_clean"),
          N = N,
          info.filter = 0, 
          maf.filter = 0.01) 
  }
  
}

##############################################################
# ===== Re-run LDSC regression including GCed sumstats ===== #
##############################################################

#The non-cleaned ones have been copied to the same data folder as the others manually!

input_vect <- c("bmi.sumstats.gz", "whr_clean.sumstats.gz","whradjbmi_clean.sumstats.gz", "pcos_doctor_clean.sumstats.gz",  "pcos_broad_clean.sumstats.gz", "pcos_consortium_clean.sumstats.gz", "pcos_meta_analysis.sumstats.gz", "pcosadjage.sumstats.gz", "pcosadjagebmi.sumstats.gz")

LDSCoutput_SUD <- ldsc(traits = input_vect, sample.prev =  c(NA, NA, NA, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5), population.prev = c(NA, NA, NA, 0.0038, 0.0084, 0.0953, 0.0700, 0.0186, 0.0186), trait.names = c("bmi", "whr", "whradjbmi", "pcos_doctor", "pcos_broad", "pcos_consortium", "pcos_meta_analysis", "pcosadjage", "pcosadjagebmi"),
                       ld=("../../../../raw_data/ldsc-master/eur_w_ld_chr.tar/eur_w_ld_chr"), ## set the directory where the LD scores used the used in the analysis are located
                       wld=("../../../../raw_data/ldsc-master/eur_w_ld_chr.tar/eur_w_ld_chr"), stand = TRUE, ldsc.log = "test")
# Format decimals
correlationSUD=round(LDSCoutput_SUD$S_Stand, 3) 

#whradjbmi pcos_doctor pcos_broad pcos_consortium pcos_meta_analysis pcosadjage pcosadjagebmi
#[1,]     1.000       0.132      0.125           0.045              0.109      0.217         0.272

dir.create("../3_gc/")

saveRDS(LDSCoutput_SUD, "../3_gc/genetic_correlation_res_corrected.RDS")

