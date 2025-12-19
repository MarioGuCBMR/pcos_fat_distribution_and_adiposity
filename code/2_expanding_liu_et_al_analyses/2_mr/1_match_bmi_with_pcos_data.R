##############
#INTRODUCTION#
##############

#This code matches BMI female IVs utilized by Liu et al with PCOS data from 6 distinct GWAS.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)
library(TwoSampleMR)

###################
#Loading functions#
###################

source("your_path/code/0_functions/functions_4_mediation.R")

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
  
  #Now we can proceed
  
  exposure <- query_ss %>%
    select(variant, chromosome, base_pair_location, effect_allele, other_allele, effect_allele_frequency, beta, standard_error, p_value, sample_size, chr_pos)
  
  colnames(exposure) <- c("SNP", "chr.exposure", "pos.exposure", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure", "samplesize.exposure", "chr_pos.exposure")
  
  exposure$exposure <- "exposure"
  exposure$id.exposure <- "exposure"
  
  #Now with the outcome:
  
  check <- which(colnames(other_ss) == "effect_allele_frequency")
  
  if(is_empty(check)){
    
    other_ss$effect_allele_frequency <- NA
    
  }
  
  #Now we can proceed
  
  if("sample_cases"%in%colnames(other_ss)){
    
    outcome <- other_ss %>%
      select(variant, chromosome, base_pair_location, effect_allele, other_allele, effect_allele_frequency, beta, standard_error, p_value, sample_size, sample_cases, sample_controls, prevalence, chr_pos)
    
    colnames(outcome) <- c("SNP", "chr.outcome", "pos.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome", "ncase.outcome", "ncontrol.outcome", "prevalence.outcome", "chr_pos.outcome")
    
  } else {
    
    outcome <- other_ss %>%
      select(variant, chromosome, base_pair_location, effect_allele, other_allele, effect_allele_frequency, beta, standard_error, p_value, sample_size, chr_pos)
    
    colnames(outcome) <- c("SNP", "chr.outcome", "pos.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome", "chr_pos.outcome")
    
  }
  
  outcome$outcome <- "outcome"
  outcome$id.outcome <- "outcome"
  
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
      select(variant, chromosome, base_pair_location, effect_allele, other_allele, effect_allele_frequency, beta, standard_error, p_value, sample_size, chr_pos)
    
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

path_2_input <- "your_path"

setwd(path_2_input)

bmi <- fread("output/1_curated_data/bmi_curated_female.txt")
all_variants <- fread("output/1_curated_data/bmi_female_ivs.txt")

#Also the whole dataframe to get proxies if necessary:

exposure_df <- bmi

proxies <- fread("output/2_replicating_liu_et_al/2_mr/0_proxies/bmi_female_ivs_proxies.txt", fill = TRUE)

proxies <- proxies[which(proxies$query_snp_rsid%in%all_variants$variant),]

length(which(duplicated(proxies$query_snp_rsid) == FALSE)) #278

#####################################################################################
#Let's get the data in build37 for the proxies cuz it makes our lives easier, really#
#####################################################################################

library(GenomicRanges)
library(rtracklayer)

#Before converting stuff, we may need to clean some proxies. Haploreg structure is not the best:

proxies = proxies[which(is.na(proxies$chr) == FALSE),]
proxies$chr=unlist(str_replace(proxies$chr, "Array", ""))
proxies = proxies[which(proxies$chr != ""),]

summary(as.numeric(proxies$chr)) #perfectly cleaned
summary(as.numeric(proxies$pos_hg38)) #perfectly cleaned

proxies$chr_pos <- NA #we will fill this with chr_pos in build 37

chain <- import.chain("raw_data/hg38ToHg19.over.chain")

gr <- GRanges(
  seqnames = paste("chr", proxies$chr, sep = ""),
  ranges = IRanges(
    start = as.numeric(proxies$pos_hg38),
    end   = as.numeric(proxies$pos_hg38)
  ),
  strand = "*"
)

# 3) Perform liftOver
lifted <- liftOver(gr, chain)
lifted_unlisted = unlist(lifted)

len <- lengths(lifted)
idx_mapped      <- which(len == 1L)

proxies <- proxies[idx_mapped,] #all were converted

converted_df <- as.data.frame(lifted_unlisted)

proxies$chr_pos=paste(converted_df$seqnames, ":", converted_df$start, sep = "")

length(unique(proxies$query_snp_rsid)) #still 278

#fwrite(proxies, "output/2_replicating_liu_et_al/2_mr/0_proxies/bmi_proxies_df.txt")

proxies <- fread("output/2_replicating_liu_et_al/2_mr/0_proxies/bmi_proxies_df.txt")

###############################
#Let's prepare the output data#
###############################

dir.create("output/2_replicating_liu_et_al/2_mr/1_expo_outcome_df/")

output_path <- "output/2_replicating_liu_et_al/2_mr/1_expo_outcome_df/"

######################
#Loading outcome data#
######################

pcos_ss <- fread("output/1_curated_data/pcos_finngen_curated.txt")
pcos_broad_ss <- fread("output/1_curated_data/pcos_finngen_broad_curated.txt")
pcos_consortium_ss <- fread("output/1_curated_data/pcos_finngen_consortium_curated.txt")
pcos_day <- fread("output/1_curated_data/pcos_day_curated.txt")
pcos_adj_age <- fread("output/1_curated_data/pcos_adj_age_curated.txt")
pcos_adj_age_bmi <- fread("output/1_curated_data/pcosadjbmi_curated.txt")

#####################################################################################
#STEP 1: Let's clean the data a little bit so that we can run MR a bit more smoothly#
#####################################################################################

pcos_ss$prevalence<- as.numeric(pcos_ss$prevalence)/100
pcos_broad_ss$prevalence<- as.numeric(pcos_broad_ss$prevalence)/100
pcos_consortium_ss$prevalence<- as.numeric(pcos_consortium_ss$prevalence)/100

pcos_day$prevalence<- as.numeric(pcos_day$prevalence)/100
pcos_adj_age$prevalence<- as.numeric(pcos_adj_age$prevalence)/100
pcos_adj_age_bmi$prevalence<- as.numeric(pcos_adj_age_bmi$prevalence)/100

######################################
#STEP 2: Let's load the data properly#
######################################

list_of_phenos <- list(pcos_ss,
                       pcos_broad_ss,
                       pcos_consortium_ss, pcos_day, pcos_adj_age, pcos_adj_age_bmi)

trait_names <- c("PCOS",
                 "PCOS (Broad)",
                 "PCOS (Consortium)",
                 "PCOS (Day et al)", "PCOS (adj age)", "PCOS (adj age+BMI)")

type_of_trait <- rep("PCOS", 6)

########################################################
#STEP 1: let's align the cluster to the positive allele#
########################################################

data_df <- bmi[which(bmi$variant%in%all_variants$variant),]

#And align them to the positive allele:
  
new_a1 <- ifelse(as.numeric(data_df$beta) < 0, data_df$other_allele, data_df$effect_allele)
new_a2 <- ifelse(as.numeric(data_df$beta) < 0, data_df$effect_allele, data_df$other_allele)
new_beta <- ifelse(as.numeric(data_df$beta) < 0, as.numeric(data_df$beta)*(-1), as.numeric(data_df$beta))
new_eaf <- ifelse(as.numeric(data_df$beta) < 0, 1-as.numeric(data_df$effect_allele_frequency), as.numeric(data_df$effect_allele_frequency))

data_df$effect_allele <- new_a1
data_df$other_allele <- new_a2
data_df$beta <- new_beta
data_df$effect_allele_frequency <- new_eaf

#######################################
#STEP 2: recursively match ir variants#
#######################################

#Let's do it for all:

recursively_matching_data(exp_df = data_df, list_of_phenos_ =list_of_phenos, trait_names_ = trait_names, output_path <- output_path, exposure = "bmi_liu_")

