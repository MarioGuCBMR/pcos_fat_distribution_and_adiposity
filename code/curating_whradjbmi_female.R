##############
#INTRODUCTION#
##############

#This code is to curate whradjbmi data from Pulit et al 2019.

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

whradjbmi<-fread("whradjbmi.giant-ukbb.meta-analysis.females.23May2018.txt.gz", stringsAsFactors = FALSE, header=TRUE)

#We need to be careful, some SNPs from Pulit et al do not present correct chr_pos.
#The mistake is in those without INFO data, thus, the ones in GIANT data.
#We are gonna remove those.

summary(whradjbmi$INFO) 
summary(whradjbmi$Freq_Tested_Allele) 
summary(whradjbmi$CHR) 
summary(whradjbmi$N) 

whradjbmi <- whradjbmi[which(is.na(whradjbmi$INFO) == FALSE),]

summary(whradjbmi$INFO) #all good
summary(whradjbmi$Freq_Tested_Allele) #all good
summary(whradjbmi$CHR) #all good
summary(whradjbmi$N) #all good

#It seems all is OK.
#Now we need to clean the RSIDs, because the have the effect allele and the other allele.
#Which is actually fine.

whradjbmi$RSID <- as.character(unlist(sapply(whradjbmi$SNP, parse_rsid)))

#Let's make a copy, just in case:

whradjbmi_data_copy <- whradjbmi

##########
#CURATION#
##########

whradjbmi_data_copy <- whradjbmi_data_copy[which(whradjbmi_data_copy$N > 10000),]

#We don't have INFO, so we cannot rely on that. 
#The MAF cannot be done, we do not have this data.
#Proceed carefully. 

whradjbmi_data_maf <- whradjbmi_data_copy[which(whradjbmi_data_copy$Freq_Tested_Allele > 0.01),]
whradjbmi_data_maf <- whradjbmi_data_maf[which(whradjbmi_data_maf$Freq_Tested_Allele < 0.99),]

#2. Now we remove the MHC region:

whradjbmi_data_maf_mhc <- whradjbmi_data_maf[which(as.numeric(whradjbmi_data_maf$CHR) == 6 & as.numeric(whradjbmi_data_maf$POS) >= 26000000 & as.numeric(whradjbmi_data_maf$POS) <= 34000000),]

summary(as.numeric(whradjbmi_data_maf_mhc$CHR)) #perfect.
summary(as.numeric(whradjbmi_data_maf_mhc$POS)) #perfect.

#Now let's check if we had any interesting variants there:

whradjbmi_data_maf_no_mhc <- whradjbmi_data_maf[which(!(whradjbmi_data_maf$SNP%in%whradjbmi_data_maf_mhc$SNP)),]

#3. #Let's get chr_pos:

whradjbmi_data_maf_no_mhc$chr_pos <- paste("chr", whradjbmi_data_maf_no_mhc$CHR, ":", whradjbmi_data_maf_no_mhc$POS, sep = "")

whradjbmi_end <- whradjbmi_data_maf_no_mhc %>%
  select(RSID, CHR, POS, Tested_Allele, Other_Allele, Freq_Tested_Allele,  BETA, SE, P, N, chr_pos)

colnames(whradjbmi_end) <-c("variant", "chromosome","base_pair_location", "effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value", "sample_size", "chr_pos")

######################################################################################
#Thus, let's save the data and match as much as possible with chromosome and position#
######################################################################################

#In the whole process we will be very careful if one of the SNPs RSID is not there.
#Because if that is the case we will have a problem when clumping.

fwrite(whradjbmi_end, "../output/1_curated_data/whradjbmi_curated_female.txt")
