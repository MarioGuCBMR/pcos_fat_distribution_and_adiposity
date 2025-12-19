##############
#INTRODUCTION#
##############

#This code is to curate whr data from Pulit et al 2019.

###################
#Loading libraries#
###################

library(tidyverse)
library(data.table)
library(TwoSampleMR)

###################
#Loading functions#
###################

parse_rsid <- function(rsid){
  
  rsid_ <- strsplit(rsid, ":")[[1]][1] #this works for those that do not present the weird formatting :)
  
  return(rsid_)
  
}

##############
#Loading data#
##############

path_input <- "your_path/raw_data/"

setwd(path_input)

whr<-fread("whr.giant-ukbb.meta-analysis.females.23May2018.txt.gz", stringsAsFactors = FALSE, header=TRUE)

#We need to be careful, some SNPs from Pulit et al do not present correct chr_pos.
#The mistake is in those without INFO data, thus, the ones in GIANT data.
#We are gonna remove those.

summary(whr$INFO) 
summary(whr$Freq_Tested_Allele) 
summary(whr$CHR) 
summary(whr$N) 

whr <- whr[which(is.na(whr$INFO) == FALSE),]

summary(whr$INFO) #all good
summary(whr$Freq_Tested_Allele) #all good
summary(whr$CHR) #all good
summary(whr$N) #all good

#It seems all is OK.
#Now we need to clean the RSIDs, because the have the effect allele and the other allele.
#Which is actually fine.

whr$RSID <- as.character(unlist(sapply(whr$SNP, parse_rsid)))

#Let's make a copy, just in case:

whr_data_copy <- whr

##########
#CURATION#
##########

whr_data_copy <- whr_data_copy[which(whr_data_copy$N > 10000),]

#We don't have INFO, so we cannot rely on that. 
#The MAF cannot be done, we do not have this data.
#Proceed carefully. 

whr_data_maf <- whr_data_copy[which(whr_data_copy$Freq_Tested_Allele > 0.01),]
whr_data_maf <- whr_data_maf[which(whr_data_maf$Freq_Tested_Allele < 0.99),]

#2. Now we remove the MHC region:

whr_data_maf_mhc <- whr_data_maf[which(as.numeric(whr_data_maf$CHR) == 6 & as.numeric(whr_data_maf$POS) >= 26000000 & as.numeric(whr_data_maf$POS) <= 34000000),]

summary(as.numeric(whr_data_maf_mhc$CHR)) #perfect.
summary(as.numeric(whr_data_maf_mhc$POS)) #perfect.

#Now let's check if we had any interesting variants there:

whr_data_maf_no_mhc <- whr_data_maf[which(!(whr_data_maf$SNP%in%whr_data_maf_mhc$SNP)),]

#3. #Let's get chr_pos:

whr_data_maf_no_mhc$chr_pos <- paste("chr", whr_data_maf_no_mhc$CHR, ":", whr_data_maf_no_mhc$POS, sep = "")

whr_end <- whr_data_maf_no_mhc %>%
  select(RSID, CHR, POS, Tested_Allele, Other_Allele, Freq_Tested_Allele,  BETA, SE, P, N, chr_pos)

colnames(whr_end) <-c("variant", "chromosome","base_pair_location", "effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value", "sample_size", "chr_pos")

######################################################################################
#Thus, let's save the data and match as much as possible with chromosome and position#
######################################################################################

#In the whole process we will be very careful if one of the SNPs RSID is not there.
#Because if that is the case we will have a problem when clumping.

fwrite(whr_end, "../output/1_curated_data/whr_curated_female.txt")

