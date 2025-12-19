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

path_2_input <- "your_path/"

setwd(path_2_input)

#############################################################################
#STEP 1: let's read the original curated data and prepare it for the munging#
#############################################################################

bmi <- fread("output/1_curated_data/bmi_curated_female.txt")
whr <- fread("output/1_curated_data/whr_curated_female.txt")
whradjbmi <- fread("output/1_curated_data/whradjbmi_curated_female.txt")
pcos_ss <- fread("output/1_curated_data/pcos_finngen_curated.txt")
pcos_broad_ss <- fread("output/1_curated_data/pcos_finngen_broad_curated.txt")
pcos_consortium_ss <- fread("output/1_curated_data/pcos_finngen_consortium_curated.txt")
pcos_day <- fread("output/1_curated_data/pcos_day_curated.txt")
pcos_adj_age <- fread("output/1_curated_data/pcos_adj_age_curated.txt")
pcos_adj_age_bmi <- fread("output/1_curated_data/pcosadjbmi_curated.txt")

#Let's change prevalences, which are all in %

pcos_ss$prevalence<- as.numeric(pcos_ss$prevalence)/100
pcos_broad_ss$prevalence<- as.numeric(pcos_broad_ss$prevalence)/100
pcos_consortium_ss$prevalence<- as.numeric(pcos_consortium_ss$prevalence)/100
pcos_day$prevalence<- as.numeric(pcos_day$prevalence)/100
pcos_adj_age$prevalence<- as.numeric(pcos_adj_age$prevalence)/100
pcos_adj_age_bmi$prevalence<- as.numeric(pcos_adj_age_bmi$prevalence)/100

#################################################
#Let's solve some issues with allele frequencies#
#################################################

finn_match <- pcos_ss[which(pcos_ss$variant%in%pcos_adj_age$variant),]
finn_match <- finn_match[which(duplicated(finn_match$variant) == FALSE),] #we do not care about duplicates here, only to have the same variant

#Let's be careful and remove duplicates in pcos_adj_age too:

pcos_adj_age = pcos_adj_age[order(as.numeric(pcos_adj_age$p_value)),]
pcos_adj_age = pcos_adj_age[which(duplicated(pcos_adj_age$variant) == FALSE),]

#Let's compute the maf

finn_match$maf <- ifelse(as.numeric(finn_match$effect_allele_frequency) > 0.5, 1-as.numeric(finn_match$effect_allele_frequency), as.numeric(finn_match$effect_allele_frequency))

#Let's order the data:

finn_match <- finn_match[order(match(finn_match$variant, pcos_adj_age$variant)),]

length(which(finn_match$variant == pcos_adj_age$variant)) #8014219! perfect

pcos_adj_age$effect_allele_frequency <- finn_match$maf

###################################################
#Let's clean the pcos data from Estonia a bit more# first for BMI and age-adjusted
###################################################

finn_match <- pcos_ss[which(pcos_ss$variant%in%pcos_adj_age_bmi$variant),]
finn_match <- finn_match[which(duplicated(finn_match$variant) == FALSE),] #we do not care about duplicates here, only to have the same variant

pcos_adj_age_bmi = pcos_adj_age_bmi[order(as.numeric(pcos_adj_age_bmi$p_value)),]
pcos_adj_age_bmi = pcos_adj_age_bmi[which(duplicated(pcos_adj_age_bmi$variant) == FALSE),]

#Let's compute the maf

finn_match$maf <- ifelse(as.numeric(finn_match$effect_allele_frequency) > 0.5, 1-as.numeric(finn_match$effect_allele_frequency), as.numeric(finn_match$effect_allele_frequency))

#Let's order the data:

finn_match <- finn_match[order(match(finn_match$variant, pcos_adj_age_bmi$variant)),]

length(which(finn_match$variant == pcos_adj_age_bmi$variant)) #8009890! perfect

pcos_adj_age_bmi$effect_allele_frequency <- finn_match$maf

######################################
#Now let's get this ready for munging#
######################################

bmi_4_munging <- curated_2_munging(bmi)
whr_4_munging <- curated_2_munging(whr)
whradjbmi_4_munging <- curated_2_munging(whradjbmi)

#Disease traits:

pcos_4_munging <- curated_2_binary(pcos_ss)
pcos_broad_4_munging <- curated_2_binary(pcos_broad_ss)
pcos_consoritum_4_munging <- curated_2_binary(pcos_consortium_ss)
pcos_day_4_munging <- curated_2_binary(pcos_day)
pcos_adj_age_4_munging <- curated_2_binary(pcos_adj_age)
pcos_adj_bmi_4_munging <- curated_2_binary(pcos_adj_age_bmi)

dir.create("output/2_replicating_liu_et_al/")
dir.create("output/2_replicating_liu_et_al/1_genetic_correlation/")
dir.create("output/2_replicating_liu_et_al/1_genetic_correlation/1_data_for_munging/")

#Let's save the data:

fwrite(bmi_4_munging, "output/2_replicating_liu_et_al/1_genetic_correlation/1_data_for_munging/bmi_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)
fwrite(whr_4_munging, "output/2_replicating_liu_et_al/1_genetic_correlation/1_data_for_munging/whr_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)
fwrite(whradjbmi_4_munging, "output/2_replicating_liu_et_al/1_genetic_correlation/1_data_for_munging/whradjbmi_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)

fwrite(pcos_4_munging, "output/2_replicating_liu_et_al/1_genetic_correlation/1_data_for_munging/pcos_doctor_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)
fwrite(pcos_broad_4_munging, "output/2_replicating_liu_et_al/1_genetic_correlation/1_data_for_munging/pcos_broad_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)
fwrite(pcos_consoritum_4_munging, "output/2_replicating_liu_et_al/1_genetic_correlation/1_data_for_munging/pcos_consortium_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)
fwrite(pcos_day_4_munging, "output/2_replicating_liu_et_al/1_genetic_correlation/1_data_for_munging/pcos_meta_analysis_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)
fwrite(pcos_adj_age_4_munging, "output/2_replicating_liu_et_al/1_genetic_correlation/1_data_for_munging/pcosadjage_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)
fwrite(pcos_adj_bmi_4_munging, "output/2_replicating_liu_et_al/1_genetic_correlation/1_data_for_munging/pcosadjbmi_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)

###############################
#STEP 2: let's run the munging#
###############################

dir.create("output/2_replicating_liu_et_al/1_genetic_correlation/2_data_munged_for_ldsc_together") #we had to change wd several times because GSEM works weirdly with tmp files.

setwd("output/2_replicating_liu_et_al/1_genetic_correlation/2_data_munged_for_ldsc_together") #this works

files<-c("../1_data_for_munging/bmi_4_munging.txt",
         "../1_data_for_munging/whr_4_munging.txt",
         "../1_data_for_munging/whradjbmi_4_munging.txt",
         "../1_data_for_munging/pcos_doctor_4_munging.txt",
         "../1_data_for_munging/pcos_broad_4_munging.txt",
         "../1_data_for_munging/pcos_consortium_4_munging.txt",
         "../1_data_for_munging/pcos_meta_analysis_4_munging.txt",
         "../1_data_for_munging/pcosadjage_4_munging.txt",
         "../1_data_for_munging/pcosadjbmi_4_munging.txt")

#Let's add the hapmap files:

hm3<-"../../../../raw_data/ldsc-master/eur_w_ld_chr.tar/eur_w_ld_chr/w_hm3.snplist"

#name the traits 
trait.names<-c("bmi", "whr", "whradjbmi", "pcos_doctor", "pcos_broad", 
               "pcos_consortium", "pcos_meta_analysis", "pcosadjage", "pcosadjagebmi")

#definte the imputation quality filter
info.filter=0

#define the MAF filter
maf.filter=0.01

#run munge

munge(files=files,hm3=hm3,trait.names=trait.names,info.filter=info.filter,maf.filter=maf.filter, overwrite = TRUE, log.name="munging_test_1")

###########################################
#STEP 3: let's run the genetic correlation#
###########################################

input_vect <- c("bmi.sumstats.gz", "whr.sumstats.gz", "whradjbmi.sumstats.gz", "pcos_doctor.sumstats.gz",  "pcos_broad.sumstats.gz", "pcos_consortium.sumstats.gz", "pcos_meta_analysis.sumstats.gz", "pcosadjage.sumstats.gz", "pcosadjagebmi.sumstats.gz")

LDSCoutput_SUD <- ldsc(traits = input_vect, sample.prev =  c(NA, NA, NA, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5), population.prev = c(NA, NA, NA, 0.0038, 0.0084, 0.0953, 0.0700, 0.0186, 0.0186), trait.names = c("bmi", "whr", "whradjbmi", "pcos_doctor", "pcos_broad", "pcos_consortium", "pcos_meta_analysis", "pcosadjage", "pcosadjagebmi"),
                       ld=("../../../../raw_data/ldsc-master/eur_w_ld_chr.tar/eur_w_ld_chr"), ## set the directory where the LD scores used the used in the analysis are located
                       wld=("../../../../raw_data/ldsc-master/eur_w_ld_chr.tar/eur_w_ld_chr"), stand = TRUE)
# Format decimals
correlationSUD=round(LDSCoutput_SUD$S_Stand, 3) 

#whradjbmi pcos_doctor pcos_broad pcos_consortium pcos_meta_analysis pcosadjage pcosadjagebmi
#[1,]     1.000       0.133      0.125           0.048              0.110      0.218         0.272

dir.create("../3_gc/")

saveRDS(LDSCoutput_SUD, "../3_gc/genetic_correlation_res.RDS")
