##############
#INTRODUCTION#
##############

#This code is to curate female BMI data from Pulit et al 2019.

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

bmi<-fread("bmi.giant-ukbb.meta-analysis.females.23May2018.txt.gz", stringsAsFactors = FALSE, header=TRUE)

#We need to be careful, some SNPs from Pulit et al do not present correct chr_pos.
#The mistake is in those without INFO data, thus, the ones in GIANT data.
#We are gonna remove those.

summary(bmi$INFO) 
summary(bmi$Freq_Tested_Allele) 
summary(bmi$CHR) 
summary(bmi$N) 

bmi <- bmi[which(is.na(bmi$INFO) == FALSE),]

summary(bmi$INFO) #all good
summary(bmi$Freq_Tested_Allele) #all good
summary(bmi$CHR) #all good
summary(bmi$N) #all good

#It seems all is OK.
#Now we need to clean the RSIDs, because the have the effect allele and the other allele.
#Which is actually fine.

bmi$RSID <- as.character(unlist(sapply(bmi$SNP, parse_rsid)))

#Let's make a copy, just in case:

bmi_data_copy <- bmi

##########
#CURATION#
##########

bmi_data_copy <- bmi_data_copy[which(bmi_data_copy$N > 10000),]

#We don't have INFO, so we cannot rely on that. 
#The MAF cannot be done, we do not have this data.
#Proceed carefully. 

bmi_data_maf <- bmi_data_copy[which(bmi_data_copy$Freq_Tested_Allele > 0.01),]
bmi_data_maf <- bmi_data_maf[which(bmi_data_maf$Freq_Tested_Allele < 0.99),]

#2. Now we remove the MHC region:

bmi_data_maf_mhc <- bmi_data_maf[which(as.numeric(bmi_data_maf$CHR) == 6 & as.numeric(bmi_data_maf$POS) >= 26000000 & as.numeric(bmi_data_maf$POS) <= 34000000),]

summary(as.numeric(bmi_data_maf_mhc$CHR)) #perfect.
summary(as.numeric(bmi_data_maf_mhc$POS)) #perfect.

#Now let's check if we had any interesting variants there:

bmi_data_maf_no_mhc <- bmi_data_maf[which(!(bmi_data_maf$SNP%in%bmi_data_maf_mhc$SNP)),]

#3. #Let's get chr_pos:

bmi_data_maf_no_mhc$chr_pos <- paste("chr", bmi_data_maf_no_mhc$CHR, ":", bmi_data_maf_no_mhc$POS, sep = "")

bmi_end <- bmi_data_maf_no_mhc %>%
  select(RSID, CHR, POS, Tested_Allele, Other_Allele, Freq_Tested_Allele,  BETA, SE, P, N, chr_pos)

colnames(bmi_end) <-c("variant", "chromosome","base_pair_location", "effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value", "sample_size", "chr_pos")

######################################################################################
#Thus, let's save the data and match as much as possible with chromosome and position#
######################################################################################

#In the whole process we will be very careful if one of the SNPs RSID is not there.
#Because if that is the case we will have a problem when clumping.

fwrite(bmi_end, "../output/1_curated_data/bmi_curated_female.txt")
