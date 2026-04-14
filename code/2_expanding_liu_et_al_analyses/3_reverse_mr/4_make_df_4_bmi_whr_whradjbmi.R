##############
#INTRODUCTION#
##############

#This code overlaps the PCOS leads with BMI data to perform MR.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)
library(TwoSampleMR)

###################
#Loading functions#
###################

source("/projects/kilpelainen-AUDIT/people/zlc436/PCOS_2026/code/0_functions/functions_4_mediation.R")

chr_cleaner <- function(chr_pos){
  
  chr_ <- str_split(chr_pos, ":")[[1]][1]
  chr_ <- str_split(chr_, "chr")[[1]][2]
  
  return(as.numeric(chr_))
}

pos_cleaner <- function(chr_pos){
  
  pos_ <- str_split(chr_pos, ":")[[1]][2]
  
  return(as.numeric(pos_))

}

data_aligner <- function(query_ss, other_ss){
  
  ################################################################################
  #This code uses: exposure and proxy dataframes that should be loaded beforehand#
  ################################################################################
  
  #########################################################################################
  #STEP 0: let's run the matching with TwoSampleMR so we need the data in a certain format#
  #########################################################################################
  
  #Let's first check if we have the effect_allele_frequency column:
  
  check <- which(colnames(query_ss) == "effect_allele_frequency")
  
  if(is_empty(check)){
    
    query_ss$effect_allele_frequency <- NA
    
  }
  
  #Now with the outcome:
  
  check <- which(colnames(other_ss) == "effect_allele_frequency")
  
  if(is_empty(check)){
    
    other_ss$effect_allele_frequency <- NA
    
  }
  
  #Now we can proceed
  
  outcome <- other_ss %>%
    dplyr::select(variant, chromosome, base_pair_location, effect_allele, other_allele, effect_allele_frequency, beta, standard_error, p_value, sample_size, chr_pos)
  
  colnames(outcome) <- c("SNP", "chr.outcome", "pos.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome", "chr_pos.outcome")
  
  outcome$outcome <- "outcome"
  outcome$id.outcome <- "outcome"
  
  #Now we can proceed
  
  if("sample_cases"%in%colnames(query_ss)){
    
    exposure <- query_ss %>%
      dplyr::select(variant, chromosome, base_pair_location, effect_allele, other_allele, effect_allele_frequency, beta, standard_error, p_value, sample_size, sample_cases, sample_controls, prevalence, chr_pos)
    
    colnames(exposure) <- c("SNP", "chr.exposure", "pos.exposure", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure", "samplesize.exposure", "ncase.exposure", "ncontrol.exposure", "prevalence.exposure", "chr_pos.exposure")
    
  } else {
    
    exposure <- query_ss %>%
      dplyr::select(variant, chromosome, base_pair_location, effect_allele, other_allele, effect_allele_frequency, beta, standard_error, p_value, sample_size, chr_pos)
    
    colnames(exposure) <- c("SNP", "chr.exposure", "pos.exposure", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure", "samplesize.exposure", "chr_pos.exposure")
    
  }
  
  exposure$exposure <- "exposure"
  exposure$id.exposure <- "exposure"
  
  ############################################################################################################
  #STEP 1: match the data. This will probably fail with tri-allelic SNPs. Here I think we are gonna be OK... #
  ############################################################################################################
  
  merged_df <- harmonise_data(exposure, outcome, action=3)
  print(dim(merged_df))
  
  merged_df <- merged_df[which(merged_df$remove == FALSE),] #removing incompatible alleles
  
  ######################################
  #STEP 2: do we need to query proxies?#
  ######################################
  
  missing_snps <- query_ss[which(!(query_ss$chr_pos%in%merged_df$chr_pos.exposure)),] 
  
  if(length(missing_snps$chr_pos) == 0){
    
    return(merged_df)
    
  } else {
    
    #We have to deal with proxies!!
    
    #1: check proxies for missing variants:
    
    proxies_4_missing <- proxies[which(proxies$query_snp_rsid%in%missing_snps$variant),]
    
    #2: check which are available in outcome:
    
    proxies_in_outcome <- proxies_4_missing[which(proxies_4_missing$rsID%in%other_ss$variant),]
    
    #3. Let's order the proxies by exposure data:
    
    exposure_match <- exposure_df[which(exposure_df$chr_pos%in%proxies_in_outcome$chr_pos),]
    exposure_match <- exposure_match[order(exposure_match$p_value),] #ordering by p-value to avoid triallelic issues
    exposure_match <- exposure_match[which(duplicated(exposure_match$variant) == FALSE),]
    
    #4. Let's add lead SNP so that we can properly add data:
    
    proxies_match <- proxies_in_outcome[which(proxies_in_outcome$chr_pos%in%exposure_match$chr_pos),]
    proxies_ordered <- proxies_match[order(match(proxies_match$chr_pos, exposure_match$chr_pos)),]
    
    length(which(proxies_ordered$chr_pos == exposure_match$chr_pos))
    
    exposure_match$lead_snp <- proxies_ordered$query_snp_rsid
    exposure_match$variant <- proxies_ordered$rsID #adding RSID to have ALL INFO. Maybe the summary statistics only has CHR:POS.
    exposure_match$rsq <- proxies_ordered$r2 #adding RSID to have ALL INFO. Maybe the summary statistics only has CHR:POS.
    
    #5. Let's get only one hit per signal, retaining only the best:
    
    exposure_best <- exposure_match[order(as.numeric(exposure_match$rsq), decreasing = TRUE),]
    exposure_best <- exposure_best[which(duplicated(exposure_best$lead_snp) == FALSE),] #7!
    
    #6. Let's align to the positive allele:
    
    tmp_df <- exposure_best

    #And align them to the positive allele:
    
    new_a1 <- ifelse(as.numeric(tmp_df$beta < 0), tmp_df$other_allele, tmp_df$effect_allele)
    new_a2 <- ifelse(as.numeric(tmp_df$beta < 0), tmp_df$effect_allele, tmp_df$other_allele)
    new_beta <- ifelse(as.numeric(tmp_df$beta < 0), as.numeric(tmp_df$beta)*(-1), as.numeric(tmp_df$beta))
    new_eaf <- ifelse(as.numeric(tmp_df$beta < 0), 1-as.numeric(tmp_df$effect_allele_frequency), as.numeric(tmp_df$effect_allele_frequency))
    
    tmp_df$effect_allele <- new_a1
    tmp_df$other_allele <- new_a2
    tmp_df$beta <- new_beta
    tmp_df$effect_allele_frequency <- new_eaf
    
    #7. Let's get the data ready for merging...
    
    exposure_proxies <- tmp_df %>%
      dplyr::select(variant, chromosome, base_pair_location, effect_allele, other_allele, effect_allele_frequency, beta, standard_error, p_value, sample_size, chr_pos)
    
    colnames(exposure_proxies) <- c("SNP", "chr.exposure", "pos.exposure", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure", "samplesize.exposure", "chr_pos.exposure")
    
    exposure_proxies$exposure <- "exposure"
    exposure_proxies$id.exposure <- "exposure"
    
    #And let's merge with the outcome, that should have ALL the data:
    
    merged_proxies_df <- harmonise_data(exposure_proxies, outcome, action=2)
    print(dim(merged_proxies_df))
    
    merged_proxies_df <- merged_proxies_df[which(merged_proxies_df$remove == FALSE),] #removing incompatible alleles
    
    #8. Let's add this info to the previous version:
    
    merged_df <- rbind(merged_df, merged_proxies_df)
    
    #9. Let's double-check for independence:
    
    #merged_df <- filt_ld(merged_df)
    
    return(merged_df)
    
  }
  
}

recursively_matching_data <- function(exp_df, list_of_phenos_, trait_names_, output_path, exposure){
  
  #STEP 0: let's make a list of traits that we need to be careful with:
  
  #STEP 1: let's loop and perform the matching:
  
  for(index_trait in seq(1, length(list_of_phenos_))){
    
    print(index_trait)
    
    trait_df <- list_of_phenos_[index_trait][[1]]
    
    #Let's see if we use data_aligner_37 or other:
    
    trait_name <- trait_names_[index_trait]
    
    print(trait_name)
    
    matched_data <- data_aligner(exp_df, trait_df)
    
    fwrite(matched_data, paste(output_path, exposure, trait_name, ".txt", sep =""))    
    
  }
  
}


#######################
#Loading exposure data#
#######################

#Let's get a path where the clusters are so that we can loop throught them:

path_2_input <- "/projects/kilpelainen-AUDIT/people/zlc436/PCOS_2026/"

setwd(path_2_input)

#Let's load all the exposure data:

pcos_ss <- fread("output/1_curated_data/pcos_finngen_curated.txt")
pcos_broad_ss <- fread("output/1_curated_data/pcos_finngen_broad_curated.txt")
pcos_consortium_ss <- fread("output/1_curated_data/pcos_finngen_consortium_curated.txt")
pcos_day <- fread("output/1_curated_data/pcos_day_curated.txt")
pcos_venkatesh <- fread("output/1_curated_data/pcos_venkatesh_curated.txt")
pcos_adj_age <- fread("output/1_curated_data/pcos_adj_age_curated.txt")
pcos_adj_age_bmi <- fread("output/1_curated_data/pcosadjbmi_curated.txt")

######################################################################
#STEP 1: Let's clean the data a little bit to get it through the SNPs#
######################################################################

#Let's get the prevalences:

pcos_ss$prevalence<- as.numeric(pcos_ss$prevalence)/100
pcos_broad_ss$prevalence<- as.numeric(pcos_broad_ss$prevalence)/100
pcos_consortium_ss$prevalence<- as.numeric(pcos_consortium_ss$prevalence)/100

pcos_day$prevalence<- as.numeric(pcos_day$prevalence)/100
pcos_venkatesh$prevalence<- as.numeric(pcos_venkatesh$prevalence)/100
pcos_adj_age$prevalence<- as.numeric(pcos_adj_age$prevalence)/100
pcos_adj_age_bmi$prevalence<- as.numeric(pcos_adj_age_bmi$prevalence)/100

#And also the outcome data:

bmi  <- fread("output/1_curated_data/bmi_curated_female.txt")
whr  <- fread("output/1_curated_data/whr_curated_female.txt")
whradjbmi  <- fread("output/1_curated_data/whradjbmi_curated_female.txt")

#Finally, let's get proxies:

proxies_full <- fread("output/2_replicating_liu_et_al/3_reverse_mr/0_proxies/all_proxies_dictionary_build_37_38.txt")

#####################################################################################################
#Let's prepare WHRadjBMI as if it were a list of outcomes to allow the data be used by our functions#
#####################################################################################################

list_of_phenos <- list(bmi, whr, whradjbmi)

trait_names <- c("BMI", "WHR", "WHRadjBMI")

type_of_trait <- rep("Anthropometric", 3)

#Let's also add the output path:

output_path="output/2_replicating_liu_et_al/3_reverse_mr/1_expo_outcome_df/"
dir.create(output_path)

################################
#Let's do this for PCOS Finngen#
################################

#Get clumped data and align to PCOS-increasing effects

pcos_clumped <- fread("output/2_replicating_liu_et_al/3_reverse_mr/0_ld_clumping/clumped_data/pcos_finngen_2026_clumped.txt")

data_df <- pcos_ss[which(pcos_ss$variant%in%pcos_clumped$rsid),] #7/7

#And align them to the positive allele:

new_a1 <- ifelse(as.numeric(data_df$beta) < 0, data_df$other_allele, data_df$effect_allele)
new_a2 <- ifelse(as.numeric(data_df$beta) < 0, data_df$effect_allele, data_df$other_allele)
new_beta <- ifelse(as.numeric(data_df$beta) < 0, as.numeric(data_df$beta)*(-1), as.numeric(data_df$beta))
new_eaf <- ifelse(as.numeric(data_df$beta) < 0, 1-as.numeric(data_df$effect_allele_frequency), as.numeric(data_df$effect_allele_frequency))

data_df$effect_allele <- new_a1
data_df$other_allele <- new_a2
data_df$beta <- new_beta
data_df$effect_allele_frequency <- new_eaf

#Now let's assess the "EXPOSURE" and "PROXIES" for PCOS (FinnGen):

exposure_df <- pcos_ss
proxies <- proxies_full[which(proxies_full$analysis == "PCOS (FinnGen)"),]
length(which(duplicated(proxies$query_snp_rsid) == FALSE)) #7/7

#I think we are ready!!

recursively_matching_data(exp_df = data_df, list_of_phenos_ =list_of_phenos, trait_names_ = trait_names, output_path <- output_path, exposure = "pcos_finngen_")

##############################
#Let's do this for PCOS Broad#
##############################

#Get clumped data and align to PCOS-increasing effects

pcos_clumped <- fread("output/2_replicating_liu_et_al/3_reverse_mr/0_ld_clumping/clumped_data/pcos_broad_2026_clumped.txt")

data_df <- pcos_broad_ss[which(pcos_broad_ss$variant%in%pcos_clumped$rsid),] #12/12

#And align them to the positive allele:

new_a1 <- ifelse(as.numeric(data_df$beta) < 0, data_df$other_allele, data_df$effect_allele)
new_a2 <- ifelse(as.numeric(data_df$beta) < 0, data_df$effect_allele, data_df$other_allele)
new_beta <- ifelse(as.numeric(data_df$beta) < 0, as.numeric(data_df$beta)*(-1), as.numeric(data_df$beta))
new_eaf <- ifelse(as.numeric(data_df$beta) < 0, 1-as.numeric(data_df$effect_allele_frequency), as.numeric(data_df$effect_allele_frequency))

data_df$effect_allele <- new_a1
data_df$other_allele <- new_a2
data_df$beta <- new_beta
data_df$effect_allele_frequency <- new_eaf

#Now let's assess the "EXPOSURE" and "PROXIES" for PCOS (Broad):

exposure_df <- pcos_broad_ss
proxies <- proxies_full[which(proxies_full$analysis == "PCOS (Broad)"),]
length(which(duplicated(proxies$query_snp_rsid) == FALSE)) #12/12

#I think we are ready!!

recursively_matching_data(exp_df = data_df, list_of_phenos_ =list_of_phenos, trait_names_ = trait_names, output_path <- output_path, exposure = "pcos_broad_")

###################################
#Let's do this for PCOS Consortium#
###################################

#Get clumped data and align to PCOS-increasing effects

pcos_clumped <- fread("output/2_replicating_liu_et_al/3_reverse_mr/0_ld_clumping/clumped_data/pcos_consortium_2026_clumped.txt")

data_df <- pcos_consortium_ss[which(pcos_consortium_ss$variant%in%pcos_clumped$rsid),] #5/5

#And align them to the positive allele:

new_a1 <- ifelse(as.numeric(data_df$beta) < 0, data_df$other_allele, data_df$effect_allele)
new_a2 <- ifelse(as.numeric(data_df$beta) < 0, data_df$effect_allele, data_df$other_allele)
new_beta <- ifelse(as.numeric(data_df$beta) < 0, as.numeric(data_df$beta)*(-1), as.numeric(data_df$beta))
new_eaf <- ifelse(as.numeric(data_df$beta) < 0, 1-as.numeric(data_df$effect_allele_frequency), as.numeric(data_df$effect_allele_frequency))

data_df$effect_allele <- new_a1
data_df$other_allele <- new_a2
data_df$beta <- new_beta
data_df$effect_allele_frequency <- new_eaf

#Now let's assess the "EXPOSURE" and "PROXIES" for PCOS (Consortium):

exposure_df <- pcos_consortium_ss
proxies <- proxies_full[which(proxies_full$analysis == "PCOS (Consortium)"),]
length(which(duplicated(proxies$query_snp_rsid) == FALSE)) #5/5

#I think we are ready!!

recursively_matching_data(exp_df = data_df, list_of_phenos_ =list_of_phenos, trait_names_ = trait_names, output_path <- output_path, exposure = "pcos_consortium_")

##################################
#Let's do this for PCOS Day et al#
##################################

#This is a special case, we need the IVs as they are

data_df <- fread("output/1_curated_data/pcos_day_ivs.txt")
data_df$chr_pos=paste("chr", data_df$chromosome, ":", data_df$base_pair_location, sep = "")

data_df$beta=as.numeric(data_df$beta)
data_df$standard_error=as.numeric(data_df$standard_error)
data_df$p_value=as.numeric(data_df$p_value)

#And align them to the positive allele:

new_a1 <- ifelse(as.numeric(data_df$beta) < 0, data_df$other_allele, data_df$effect_allele)
new_a2 <- ifelse(as.numeric(data_df$beta) < 0, data_df$effect_allele, data_df$other_allele)
new_beta <- ifelse(as.numeric(data_df$beta) < 0, as.numeric(data_df$beta)*(-1), as.numeric(data_df$beta))
new_eaf <- ifelse(as.numeric(data_df$beta) < 0, 1-as.numeric(data_df$effect_allele_frequency), as.numeric(data_df$effect_allele_frequency))

data_df$effect_allele <- new_a1
data_df$other_allele <- new_a2
data_df$beta <- new_beta
data_df$effect_allele_frequency <- new_eaf

#Now let's assess the "EXPOSURE" and "PROXIES" for PCOS (Consortium):

exposure_df <- pcos_day
proxies <- proxies_full[which(proxies_full$analysis == "PCOS (Day)"),]
length(which(duplicated(proxies$query_snp_rsid) == FALSE)) #14/14

#I think we are ready!!

recursively_matching_data(exp_df = data_df, list_of_phenos_ =list_of_phenos, trait_names_ = trait_names, output_path <- output_path, exposure = "pcos_day_")

########################################
#Let's do this for PCOS Venkatesh et al#
########################################

#This is a special case, we need the IVs as they are

pcos_clumped <- fread("output/2_replicating_liu_et_al/3_reverse_mr/0_ld_clumping/clumped_data/pcos_venkatesh_2026_clumped.txt")

data_df <- pcos_venkatesh[which(pcos_venkatesh$variant%in%pcos_clumped$rsid),] #17/17

#And align them to the positive allele:

new_a1 <- ifelse(as.numeric(data_df$beta) < 0, data_df$other_allele, data_df$effect_allele)
new_a2 <- ifelse(as.numeric(data_df$beta) < 0, data_df$effect_allele, data_df$other_allele)
new_beta <- ifelse(as.numeric(data_df$beta) < 0, as.numeric(data_df$beta)*(-1), as.numeric(data_df$beta))
new_eaf <- ifelse(as.numeric(data_df$beta) < 0, 1-as.numeric(data_df$effect_allele_frequency), as.numeric(data_df$effect_allele_frequency))

data_df$effect_allele <- new_a1
data_df$other_allele <- new_a2
data_df$beta <- new_beta
data_df$effect_allele_frequency <- new_eaf

#Now let's assess the "EXPOSURE" and "PROXIES" for PCOS (Venkatesh):

exposure_df <- pcos_venkatesh
proxies <- proxies_full[which(proxies_full$analysis == "PCOS (Venkatesh)"),]
length(which(duplicated(proxies$query_snp_rsid) == FALSE)) #17/17

#I think we are ready!!

recursively_matching_data(exp_df = data_df, list_of_phenos_ =list_of_phenos, trait_names_ = trait_names, output_path <- output_path, exposure = "pcos_venkatesh_")

##############################
#Let's do this for PCOS Tyrmi#
##############################

#Get clumped data and align to PCOS-increasing effects

pcos_clumped <- fread("output/2_replicating_liu_et_al/3_reverse_mr/0_ld_clumping/clumped_data/pcos_tyrmi_2026_clumped.txt")

data_df <- pcos_adj_age[which(pcos_adj_age$variant%in%pcos_clumped$rsid),] #5/5

#And align them to the positive allele:

new_a1 <- ifelse(as.numeric(data_df$beta) < 0, data_df$other_allele, data_df$effect_allele)
new_a2 <- ifelse(as.numeric(data_df$beta) < 0, data_df$effect_allele, data_df$other_allele)
new_beta <- ifelse(as.numeric(data_df$beta) < 0, as.numeric(data_df$beta)*(-1), as.numeric(data_df$beta))
new_eaf <- ifelse(as.numeric(data_df$beta) < 0, 1-as.numeric(data_df$effect_allele_frequency), as.numeric(data_df$effect_allele_frequency))

data_df$effect_allele <- new_a1
data_df$other_allele <- new_a2
data_df$beta <- new_beta
data_df$effect_allele_frequency <- new_eaf

#Now let's assess the "EXPOSURE" and "PROXIES" for PCOS (FinnGen):

exposure_df <- pcos_adj_age
proxies <- proxies_full[which(proxies_full$analysis == "PCOS (Tyrmi)"),]
length(which(duplicated(proxies$query_snp_rsid) == FALSE)) #6/6

#I think we are ready!!

recursively_matching_data(exp_df = data_df, list_of_phenos_ =list_of_phenos, trait_names_ = trait_names, output_path <- output_path, exposure = "pcos_tyrmi_")

####################################
#Let's do this for PCOSadjBMI Tyrmi#
####################################

#Get clumped data and align to PCOS-increasing effects

pcos_clumped <- fread("output/2_replicating_liu_et_al/3_reverse_mr/0_ld_clumping/clumped_data/pcosadjbmi_tyrmi_2026_clumped.txt")

data_df <- pcos_adj_age_bmi[which(pcos_adj_age_bmi$variant%in%pcos_clumped$rsid),] #3/3

#And align them to the positive allele:

new_a1 <- ifelse(as.numeric(data_df$beta) < 0, data_df$other_allele, data_df$effect_allele)
new_a2 <- ifelse(as.numeric(data_df$beta) < 0, data_df$effect_allele, data_df$other_allele)
new_beta <- ifelse(as.numeric(data_df$beta) < 0, as.numeric(data_df$beta)*(-1), as.numeric(data_df$beta))
new_eaf <- ifelse(as.numeric(data_df$beta) < 0, 1-as.numeric(data_df$effect_allele_frequency), as.numeric(data_df$effect_allele_frequency))

data_df$effect_allele <- new_a1
data_df$other_allele <- new_a2
data_df$beta <- new_beta
data_df$effect_allele_frequency <- new_eaf

#Now let's assess the "EXPOSURE" and "PROXIES" for PCOS (FinnGen):

exposure_df <- pcos_adj_age_bmi
proxies <- proxies_full[which(proxies_full$analysis == "PCOSadjBMI (Tyrmi)"),]
length(which(duplicated(proxies$query_snp_rsid) == FALSE)) #3/3

#I think we are ready!!

recursively_matching_data(exp_df = data_df, list_of_phenos_ =list_of_phenos, trait_names_ = trait_names, output_path <- output_path, exposure = "pcosadjbmi_tyrmi_")

