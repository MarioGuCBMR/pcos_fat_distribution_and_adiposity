##############
#INTRODUCTION#
##############

#This is a code to curate PCOS from FinnGen data.

#Let's go:

###########
#libraries#
###########

library(data.table)
library(tidyverse)
library(GenomicRanges)
library(rtracklayer)

#######################
#Loading pcos_day data#
#######################

#We are gonna load the pcos_day from 2024. 

setwd("/projects/kilpelainen-AUDIT/people/zlc436/PCOS_2026")

pcos=fread("raw_data/summary_stats_release_finngen_R12_E4_PCOS_CONCORTIUM.gz")

###################################################################################################################
#The variants are in build 38, let's clean the data as they are and then use the chain and Bioconductor to convert#
###################################################################################################################

colnames(pcos) = c("chromosome", "base_pair_location", "other_allele", "effect_allele", "variant", "nearest_gene", "p_value", "logp", "beta", "standard_error", "effect_allele_frequency", "alt_cases", "alt_controls")

#Let's get the columns that we want first:

pcos_clean = pcos %>%
  dplyr::select(variant, chromosome, base_pair_location, effect_allele, other_allele, effect_allele_frequency, beta, standard_error, p_value)

#Let's remove indels:

yes_vect <- c("A", "G", "C", "T")

pcos_clean_indel = pcos_clean[which(!(pcos_clean$effect_allele%in%yes_vect & pcos_clean$other_allele%in%yes_vect)),]
pcos_clean_noindel = pcos_clean[which(pcos_clean$effect_allele%in%yes_vect & pcos_clean$other_allele%in%yes_vect),]

#Let's get only autosomal data:

summary(pcos_clean_noindel$chromosome) #sexual chr as 23

pcos_autosomal = pcos_clean_noindel[which(as.numeric(pcos_clean_noindel$chromosome) != 23),]

summary(pcos_autosomal$chromosome) #sexual chr as 23

#Let's remove rare variants:

summary(as.numeric(pcos_autosomal$effect_allele_frequency))

pcos_common = pcos_autosomal[which(as.numeric(pcos_autosomal$effect_allele_frequency) > 0.01 & as.numeric(pcos_autosomal$effect_allele_frequency) < 0.99),]

summary(as.numeric(pcos_common$effect_allele_frequency)) #perfect

#####################################################################
#Before removing MHC region, we need to convert the data to build 37#
#####################################################################

chain <- import.chain("raw_data/hg38ToHg19.over.chain")

gr <- GRanges(
  seqnames = paste("chr", pcos_common$chromosome, sep = ""),
  ranges = IRanges(
    start = pcos_common$base_pair_location,
    end   = pcos_common$base_pair_location
  ),
  strand = "*"
)

# 3) Perform liftOver
lifted <- liftOver(gr, chain)
lifted_unlisted = unlist(lifted)

len <- lengths(lifted)
idx_mapped      <- which(len == 1L)

pcos_converted <- pcos_common[idx_mapped,]

converted_df <- as.data.frame(lifted_unlisted)

pcos_converted$pos_37 = converted_df$start #check - the conversion worked ;) - the order is the same since we are using the index 

pcos_converted = pcos_converted %>%
  dplyr::select(variant, chromosome, pos_37, effect_allele, other_allele, effect_allele_frequency, beta, standard_error, p_value)

colnames(pcos_converted) = c("variant", "chromosome", "base_pair_location", "effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value")

#############################
#Let's remove the MHC region#
#############################

pcos_converted$chr_pos=paste("chr", pcos_converted$chromosome, ":", pcos_converted$base_pair_location, sep="")

mhc_region = pcos_converted[which(as.numeric(pcos_converted$chromosome) == 6 & as.numeric(pcos_converted$base_pair_location) >= 26000000 & as.numeric(pcos_converted$base_pair_location) <= 34000000),]

summary(mhc_region$chromosome)
summary(mhc_region$base_pair_location)

pcos_no_mhc = pcos_converted[which(!(pcos_converted$chr_pos%in%mhc_region$chr_pos)),]

##########################################################################
#Amazing, let's add the cases, controls, total sample size and prevalence#
##########################################################################

#Data extracted from: https://r12.finngen.fi/

pcos_no_mhc$sample_cases=42630
pcos_no_mhc$sample_controls=239434
pcos_no_mhc$sample_size=42630+239434
pcos_no_mhc$prevalence=9.53 #Unadjusted period prevalence% in whole population from: https://risteys.finngen.fi/endpoints/E4_PCOS_CONCORTIUM

#We are done!

fwrite(pcos_no_mhc, "output/1_curated_data/pcos_finngen_consortium_curated.txt")



