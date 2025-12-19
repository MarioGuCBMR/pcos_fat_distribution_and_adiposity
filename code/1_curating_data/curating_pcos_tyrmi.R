##############
#INTRODUCTION#
##############

#This is a code to PCOS GWAS from Tyrmi et al meta-analysis

#Let's go:

###########
#libraries#
###########

library(data.table)
library(tidyverse)

###########################
#Loading pcos_adj_age data#
###########################

#We are gonna load the pcos_adj_age from 2024. 

setwd("your_path")

pcos_adj_age <- fread("raw_data/GCST90044902_buildGRCh37.tsv") #https://www.ebi.ac.uk/gwas/efotraits/EFO_0004703

######################
#STEP 1: GENERATE IDs#
######################

#We are gonna set the chr_pos:

pcos_adj_age$chr_pos <- paste("chr", pcos_adj_age$chromosome, ":", pcos_adj_age$base_pair_location, sep = "")
pcos_adj_age$effect_allele <- toupper(pcos_adj_age$effect_allele)
pcos_adj_age$other_allele <- toupper(pcos_adj_age$other_allele)

#Remove alleles:

yes_vect <- c("A", "G", "C", "T")

pcos_adj_age <- pcos_adj_age[which(pcos_adj_age$effect_allele%in%yes_vect & pcos_adj_age$other_allele%in%yes_vect),]

###############################
#STEP 2: CHECK FOR WEIRD STUFF#
###############################

summary(as.numeric(pcos_adj_age$chromosome)) #we have some exceptions!
pcos_adj_age <- pcos_adj_age[which(pcos_adj_age$chromosome%in%seq(1,22)),]
summary(as.numeric(pcos_adj_age$chromosome)) #only autosomal

###########################
#STEP 3: REMOVE MAF < 0.01#
###########################

#We do not have allele frequencies:

pcos_adj_age$effect_allele_frequency <- NA

pcos_adj_age_corrected_eaf <- pcos_adj_age
  
##############################
#CURATION OF THE DATA IS DONE#
##############################

#Let's check now if there are any in the MHC region.

pcos_adj_age_mhc <- pcos_adj_age_corrected_eaf[which(as.numeric(pcos_adj_age_corrected_eaf$chromosome) == 6 & as.numeric(pcos_adj_age_corrected_eaf$base_pair_location) >= 26000000 & as.numeric(pcos_adj_age_corrected_eaf$base_pair_location) <= 34000000),]

summary(as.numeric(pcos_adj_age_mhc$chromosome)) #perfect!!
summary(as.numeric(pcos_adj_age_mhc$base_pair_location)) #perfect!!

pcos_adj_age_end <- pcos_adj_age_corrected_eaf[which(!(pcos_adj_age_corrected_eaf$chr_pos%in%pcos_adj_age_mhc$chr_pos)),]

#Let's adjust the data:

pcos_adj_age_end <- pcos_adj_age_end %>%
  select("chromosome", "base_pair_location", "effect_allele", "other_allele","effect_allele_frequency", "beta", "standard_error", "p_value",  "chr_pos")

colnames(pcos_adj_age_end) <- c("chromosome", "base_pair_location", "effect_allele", "other_allele","effect_allele_frequency", "beta", "standard_error", "p_value", "chr_pos")

###############################################
#Let's add rsIDS, sample sizes and other stuff#
###############################################

pcos_adj_age_end$sample_cases <- 3609
pcos_adj_age_end$sample_controls <- 229788
pcos_adj_age_end$sample_size <- 3609+229788 #Let's use this for all instead
pcos_adj_age_end$prevalence <- (0.57+3.15)/2 #mean between both populations

#Let's get as many rsIDs from the rest:

pcos <- fread("output/1_curated_data/pcos_finngen_curated.txt")

pcos_match <- pcos[which(pcos$chr_pos%in%pcos_adj_age_end$chr_pos),]
pcos_match <- pcos_match[which(duplicated(pcos_match$chr_pos) == FALSE),] #we do not care about triallelics

pcos_tyrmi_match <- pcos_adj_age_end[which(pcos_adj_age_end$chr_pos%in%pcos_match$chr_pos),] #the matching is perfect
pcos_tyrmi_match <- pcos_tyrmi_match[which(duplicated(pcos_tyrmi_match$chr_pos) == FALSE),]

pcos_match <- pcos_match[order(match(pcos_match$chr_pos, pcos_tyrmi_match$chr_pos)),]

length(which(pcos_match$chr_pos == pcos_tyrmi_match$chr_pos))

pcos_tyrmi_match$variant <- pcos_match$variant

#########################
#We can save this data!!#
#########################

fwrite(pcos_tyrmi_match, "output/1_curated_data/pcos_adj_age_curated.txt")
