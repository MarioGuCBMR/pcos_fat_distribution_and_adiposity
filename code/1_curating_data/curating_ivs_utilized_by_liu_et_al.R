##############
#INTRODUCTION#
##############

#This code takes Liu et al supplementary table: 12916_2022_2238_MOESM1_ESM.xlsx, downloaded from: https://static-content.springer.com/esm/art%3A10.1186%2Fs12916-022-02238-y/MediaObjects/12916_2022_2238_MOESM1_ESM.xlsx, available from: https://link.springer.com/article/10.1186/s12916-022-02238-y#Sec19

#We are going to retrieve the instruments Liu et al utilized with xlsx package. The data is the following:

#BMI: Suppl. Table 1
#WHR: Suppl. Table 2
#WHRadjBMI: Suppl. Table 3
#PCOS: Suppl. Table 4.

##################
#Loading packages#
##################

library(readxl)
library(data.table)
library(tidyverse)

##############
#Loading data#
##############

setwd("your_path/")

bmi_ivs = as.data.frame(read_xlsx("raw_data/12916_2022_2238_MOESM1_ESM.xlsx", sheet = 1)) #284
whr_ivs = as.data.frame(read_xlsx("raw_data/12916_2022_2238_MOESM1_ESM.xlsx", sheet = 2)) #206
whradjbmi_ivs = as.data.frame(read_xlsx("raw_data/12916_2022_2238_MOESM1_ESM.xlsx", sheet = 3)) #269
pcos_ivs = as.data.frame(read_xlsx("raw_data/12916_2022_2238_MOESM1_ESM.xlsx", sheet = 5)) #28

########################################################
#Let's clean the data and get it ready for our analyses#
########################################################

#First. column cleaning:

colnames(bmi_ivs) <- bmi_ivs[1,] #colnames are actually in the first row

bmi_ivs <- bmi_ivs %>%
  dplyr::select("SNP", "CHR", "BP", "A1", "A2", "EAF", "Beta", "SE", "P-value") #this removes the outcome matches

#Let's take remove the rows with extra info:

bmi_ivs <- bmi_ivs[3:length(bmi_ivs$SNP),]

#Some of them have mismatches, which show up as NAs. We can remove those.

test <- bmi_ivs[which(is.na(bmi_ivs$A2) == TRUE),] #196
bmi_ivs <- bmi_ivs[which(is.na(bmi_ivs$A2) == FALSE),] #278

#Great, let's see which variants match with the curated data. 
#The curated data from full summary statistics already has been cleaned for MHC region variants, rare variants and indels. 

bmi = fread("output/1_curated_data/bmi_curated_female.txt")

#Let's find the same SNPs and call it a day:

bmi_match <- bmi[which(bmi$variant%in%bmi_ivs$SNP),] #perfect match - lovely

fwrite(bmi_match, "output/1_curated_data/bmi_female_ivs.txt")

rm(bmi)
rm(bmi_ivs)
rm(bmi_match)

###########################
#Let's do the same for WHR#
###########################

#First. column cleaning:

colnames(whr_ivs) <- whr_ivs[1,] #colnames are actually in the first row

whr_ivs <- whr_ivs %>%
  dplyr::select("SNP", "CHR", "BP", "A1", "A2", "EAF", "Beta (SD)", "SE", "P-value") #this removes the outcome matches

#Let's take remove the rows with extra info:

whr_ivs <- whr_ivs[3:length(whr_ivs$SNP),]

#Some of them have mismatches, which show up as NAs. We can remove those.

test <- whr_ivs[which(is.na(whr_ivs$A2) == TRUE),] #196

whr_ivs <- whr_ivs[which(is.na(whr_ivs$A2) == FALSE),] #196

#Great, let's see which variants match with the curated data. 
#The curated data from full summary statistics already has been cleaned for MHC region variants, rare variants and indels. 

whr = fread("output/1_curated_data/whr_curated_female.txt")

#Let's find the same SNPs and call it a day:

whr_match <- whr[which(whr$variant%in%whr_ivs$SNP),] #perfect match - lovely

fwrite(whr_match, "output/1_curated_data/whr_female_ivs.txt")

rm(whr)
rm(whr_ivs)
rm(whr_match)

#################################
#Let's do the same for WHRadjBMI#
#################################

#First. column cleaning:

colnames(whradjbmi_ivs) <- whradjbmi_ivs[1,] #colnames are actually in the first row

whradjbmi_ivs <- whradjbmi_ivs %>%
  dplyr::select("SNP", "CHR", "BP", "A1", "A2", "EAF", "Beta (SD)", "SE", "P-value") #this removes the outcome matches

#Let's take remove the rows with extra info:

whradjbmi_ivs <- whradjbmi_ivs[3:length(whradjbmi_ivs$SNP),]

#Some of them have mismatches, which show up as NAs. We can remove those.

test <- whradjbmi_ivs[which(is.na(whradjbmi_ivs$A2) == TRUE),] #196

whradjbmi_ivs <- whradjbmi_ivs[which(is.na(whradjbmi_ivs$A2) == FALSE),] #256

#Great, let's see which variants match with the curated data. 
#The curated data from full summary statistics already has been cleaned for MHC region variants, rare variants and indels. 

whradjbmi = fread("output/1_curated_data/whradjbmi_curated_female.txt")

#Let's find the same SNPs and call it a day:

whradjbmi_match <- whradjbmi[which(whradjbmi$variant%in%whradjbmi_ivs$SNP),] #perfect match - lovely

fwrite(whradjbmi_match, "output/1_curated_data/whradjbmi_female_ivs.txt")

rm(whradjbmi)
rm(whradjbmi_ivs)
rm(whradjbmi_match)

######################################
#Let's finally do this for PCOS too!!#
######################################

#First. column cleaning:

colnames(pcos_ivs) <- pcos_ivs[1,] #colnames are actually in the first row

pcos_ivs <- pcos_ivs %>%
  dplyr::select("SNP", "CHR", "BP", "A1", "A2", "EAF", "Beta  (SD)", "SE", "P-value") #this removes the outcome matches

#Some of them have mismatches, which show up as NAs. We can remove those.

pcos_ivs <- pcos_ivs[3:length(pcos_ivs$SNP),]
pcos_ivs <- pcos_ivs[which(is.na(pcos_ivs$A2) == FALSE),] #256

#In this particular case we cannot work with full data since it does not contain meta-analysis!

colnames(pcos_ivs) <- c("variant", "chromosome", "base_pair_location", "effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value")

#Let's add the same info we added in the original datasets:

ncases_w_23 <- 10074
ncontrols_w_23 <- 103164

pcos_ivs$sample_cases <- ncases_w_23
pcos_ivs$sample_controls <- ncontrols_w_23
pcos_ivs$sample_size <- ncases_w_23+ncontrols_w_23 #Let's use this for all instead
pcos_ivs$prevalence <- 7 #the usual prevalence and first stated in Day et al 

#Let's add the data

fwrite(pcos_ivs, "output/1_curated_data/pcos_day_ivs.txt")

rm(pcos)
rm(pcos_ivs)
rm(pcos_match)

