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

setwd("/projects/kilpelainen-AUDIT/people/zlc436/PCOS_2026")

pcosadjbmi <- fread("raw_data/GCST90044903_buildGRCh37.tsv.gz") #https://www.ebi.ac.uk/gwas/efotraits/EFO_0004703

######################
#STEP 1: GENERATE IDs#
######################

#We are gonna set the chr_pos:

pcosadjbmi$chr_pos <- paste("chr", pcosadjbmi$chromosome, ":", pcosadjbmi$base_pair_location, sep = "")
pcosadjbmi$effect_allele <- toupper(pcosadjbmi$effect_allele)
pcosadjbmi$other_allele <- toupper(pcosadjbmi$other_allele)

#Remove alleles:

yes_vect <- c("A", "G", "C", "T")

pcosadjbmi <- pcosadjbmi[which(pcosadjbmi$effect_allele%in%yes_vect & pcosadjbmi$other_allele%in%yes_vect),]

###############################
#STEP 2: CHECK FOR WEIRD STUFF#
###############################

summary(as.numeric(pcosadjbmi$chromosome)) #we have some exceptions!
pcosadjbmi <- pcosadjbmi[which(pcosadjbmi$chromosome%in%seq(1,22)),]
summary(as.numeric(pcosadjbmi$chromosome)) #only autosomal

###########################
#STEP 3: REMOVE MAF < 0.01#
###########################

#We do not have allele frequencies:

pcosadjbmi$effect_allele_frequency <- NA

pcosadjbmi_corrected_eaf <- pcosadjbmi
  
##############################
#CURATION OF THE DATA IS DONE#
##############################

#Let's check now if there are any in the MHC region.

pcosadjbmi_mhc <- pcosadjbmi_corrected_eaf[which(as.numeric(pcosadjbmi_corrected_eaf$chromosome) == 6 & as.numeric(pcosadjbmi_corrected_eaf$base_pair_location) >= 26000000 & as.numeric(pcosadjbmi_corrected_eaf$base_pair_location) <= 34000000),]

summary(as.numeric(pcosadjbmi_mhc$chromosome)) #perfect!!
summary(as.numeric(pcosadjbmi_mhc$base_pair_location)) #perfect!!

pcosadjbmi_end <- pcosadjbmi_corrected_eaf[which(!(pcosadjbmi_corrected_eaf$chr_pos%in%pcosadjbmi_mhc$chr_pos)),]

#Let's adjust the data:

pcosadjbmi_end <- pcosadjbmi_end %>%
  select("chromosome", "base_pair_location", "effect_allele", "other_allele","effect_allele_frequency", "beta", "standard_error", "p_value",  "chr_pos")

colnames(pcosadjbmi_end) <- c("chromosome", "base_pair_location", "effect_allele", "other_allele","effect_allele_frequency", "beta", "standard_error", "p_value", "chr_pos")

###############################################
#Let's add rsIDS, sample sizes and other stuff#
###############################################

pcosadjbmi_end$sample_cases <- 2169
pcosadjbmi_end$sample_controls <- 160321 
pcosadjbmi_end$sample_size <- 2169+160321 #Let's use this for all instead
pcosadjbmi_end$prevalence <- (0.57+3.15)/2 #mean between both populations

#Let's get as many rsIDs from the rest:

pcos <- fread("output/1_curated_data/pcos_finngen_curated.txt")

pcos_match <- pcos[which(pcos$chr_pos%in%pcosadjbmi_end$chr_pos),]
pcos_match <- pcos_match[which(duplicated(pcos_match$chr_pos) == FALSE),] #we do not care about triallelics

pcos_tyrmi_match <- pcosadjbmi_end[which(pcosadjbmi_end$chr_pos%in%pcos_match$chr_pos),] #the matching is perfect
pcos_tyrmi_match <- pcos_tyrmi_match[which(duplicated(pcos_tyrmi_match$chr_pos) == FALSE),]

pcos_match <- pcos_match[order(match(pcos_match$chr_pos, pcos_tyrmi_match$chr_pos)),]

length(which(pcos_match$chr_pos == pcos_tyrmi_match$chr_pos))

pcos_tyrmi_match$variant <- pcos_match$variant

#########################
#We can save this data!!#
#########################

fwrite(pcos_tyrmi_match, "output/1_curated_data/pcosadjbmi_curated.txt")
