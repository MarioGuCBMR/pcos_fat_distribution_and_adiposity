##############
#INTRODUCTION#
##############

#This is a code to curate PCOS from day et al.

#Let's go:

###########
#libraries#
###########

library(data.table)
library(tidyverse)

#######################
#Loading pcos_day data#
#######################

#We are gonna load the pcos_day from 2024. 

setwd("your_path")

pcos_day <- fread("raw_data/PCOS/PCOS_summary_data_19092018.txt") 

######################
#STEP 1: GENERATE IDs#
######################

#We are gonna set the chr_pos:

pcos_day$chr_pos <- paste("chr", pcos_day$MarkerName, sep = "")
chr_pos_dummy <- unlist(str_split(pcos_day$chr_pos, ":"))
chr_ <- chr_pos_dummy[which(str_detect(chr_pos_dummy, "chr"))]
chr_ <- unlist(str_split(chr_, "chr"))
chr_ <- chr_[which(chr_ != "")]
pos_ <- chr_pos_dummy[which(str_detect(chr_pos_dummy, "chr") == FALSE)]
pos_ <- pos_[which(pos_ != "ID")]

pcos_day$chromosome <- chr_
pcos_day$base_pair_location <- pos_
pcos_day$chr_pos <- paste("chr", pcos_day$chromosome, ":", pcos_day$base_pair_location, sep = "")

#Remove alleles:

pcos_day$effect_allele <- toupper(pcos_day$Effect_allele)
pcos_day$other_allele <- toupper(pcos_day$Other_allele)

yes_vect <- c("A", "G", "C", "T")

pcos_day <- pcos_day[which(pcos_day$effect_allele%in%yes_vect & pcos_day$other_allele%in%yes_vect),]

###############################
#STEP 2: CHECK FOR WEIRD STUFF#
###############################

summary(as.numeric(pcos_day$chromosome)) #We have chromosome X

pcos_day <- pcos_day[which(is.na(as.numeric(pcos_day$chromosome)) == FALSE),]

summary(as.numeric(pcos_day$chromosome)) #perfect

###########################
#STEP 3: REMOVE MAF < 0.01#
###########################

pcos_day$effect_allele_frequency <- pcos_day$EAF

summary(as.numeric(pcos_day$effect_allele_frequency)) #quite good!!

pcos_day_eaf_OK <- pcos_day[which(pcos_day$effect_allele_frequency > 0.01),]
pcos_day_eaf_OK <- pcos_day_eaf_OK[which(pcos_day_eaf_OK$effect_allele_frequency < 0.99),]

summary(pcos_day_eaf_OK$effect_allele_frequency) #worked like a charm.

pcos_day_corrected_eaf <- pcos_day_eaf_OK
  
##############################
#CURATION OF THE DATA IS DONE# I checked and data is indeed in build 37
##############################

#Let's check now if there are any in the MHC region.

pcos_day_mhc <- pcos_day_corrected_eaf[which(as.numeric(pcos_day_corrected_eaf$chromosome) == 6 & as.numeric(pcos_day_corrected_eaf$base_pair_location) >= 26000000 & as.numeric(pcos_day_corrected_eaf$base_pair_location) <= 34000000),]

summary(as.numeric(pcos_day_mhc$chromosome)) #perfect!!
summary(as.numeric(pcos_day_mhc$base_pair_location)) #perfect!!

pcos_day_end <- pcos_day_corrected_eaf[which(!(pcos_day_corrected_eaf$chr_pos%in%pcos_day_mhc$chr_pos)),]

#Let's adjust the data:

pcos_day_end <- pcos_day_end %>%
  select("chromosome", "base_pair_location", "effect_allele", "other_allele","effect_allele_frequency", "Effect", "StdErr", "Pvalue",  "TotalSampleSize", "chr_pos")

colnames(pcos_day_end) <- c("chromosome", "base_pair_location", "effect_allele", "other_allele","effect_allele_frequency", "beta", "standard_error", "p_value",  "sample_size", "chr_pos")

###################################################################################
#Let's add the cases, the controls and also the prevalence they state in the paper#
###################################################################################

ncases_wo_23 <- 1184+200+500+157+658+782+202+500+485+12189
ncontrols_wo_23 <- 5799+3500+2807+6774+2971+1000+403+2264

pcos_day_end$sample_cases <- ncases_wo_23
pcos_day_end$sample_controls <- ncontrols_wo_23
pcos_day_end$sample_size <- ncases_wo_23+ncontrols_wo_23 #Let's use this for all instead
pcos_day_end$prevalence <- 7 #the usual prevalence and first stated in Day et al 

#Let's get as many rsIDs from the rest:

pcos <- fread("output/1_curated_data/pcos_finngen_curated.txt")

pcos_match <- pcos[which(pcos$chr_pos%in%pcos_day_end$chr_pos),]
pcos_match <- pcos_match[which(duplicated(pcos_match$chr_pos) == FALSE),] #we do not care about triallelics

pcos_day_match <- pcos_day_end[which(pcos_day_end$chr_pos%in%pcos_match$chr_pos),] #the matching is perfect
#pcos_day_match <- pcos_day_match[which(duplicated(pcos_day_match$chr_pos) == FALSE),]

pcos_match <- pcos_match[order(match(pcos_match$chr_pos, pcos_day_match$chr_pos)),]

length(which(pcos_match$chr_pos == pcos_day_match$chr_pos))

pcos_day_match$variant <- pcos_match$variant

#########################
#We can save this data!!#
#########################

fwrite(pcos_day_match, "output/1_curated_data/pcos_day_curated.txt")
