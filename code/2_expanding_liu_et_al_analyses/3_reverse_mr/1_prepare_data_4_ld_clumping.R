##############
#INTRODUCTION#
##############

#This code prepares the PCOS data for LD-clumping

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)
library(TwoSampleMR)

###################
#Loading functions#
###################

formatting_data_4_clumping <- function(ss_gw){
  #This function formats the data for the clumping.
  
  #STEP 2: now make the SNP column which is SNP:A1:A2
  
  rsid <- ifelse(str_detect(ss_gw$variant, "rs") == TRUE, ss_gw$variant, ss_gw$chr_pos)
  SNP <- paste(rsid, ":", ss_gw$other_allele, ":", ss_gw$effect_allele, sep = "")
  
  #Now we select the columns that we want:
  
  ss_gw$SNP <- SNP
  ss_gw$rsid <- rsid
  
  final_df <- ss_gw %>%
    select(SNP, p_value, chr_pos, effect_allele, other_allele, rsid)
  
  colnames(final_df) <- c("SNP", "pval", "chr_pos", "effect_allele", "other_allele", "rsid")
  
  return(final_df)
  
}

###############################
#Let's prepare the output data#
###############################

setwd("/projects/kilpelainen-AUDIT/people/zlc436/PCOS_2026/")

dir.create("output/2_replicating_liu_et_al")
dir.create("output/2_replicating_liu_et_al/3_reverse_mr")
dir.create("output/2_replicating_liu_et_al/3_reverse_mr/0_ld_clumping/")

#Let's load the data:

pcos_ss <- fread("output/1_curated_data/pcos_finngen_curated.txt")
pcos_broad_ss <- fread("output/1_curated_data/pcos_finngen_broad_curated.txt")
pcos_consortium_ss <- fread("output/1_curated_data/pcos_finngen_consortium_curated.txt")
pcos_day <- fread("output/1_curated_data/pcos_day_curated.txt")
pcos_adj_age <- fread("output/1_curated_data/pcos_adj_age_curated.txt")
pcos_adj_age_bmi <- fread("output/1_curated_data/pcosadjbmi_curated.txt")

######################################################################
#STEP 1: Let's clean the data a little bit to get it through the SNPs#
######################################################################

#Let's get the prevalences right:

pcos_ss$prevalence<- as.numeric(pcos_ss$prevalence)/100
pcos_broad_ss$prevalence<- as.numeric(pcos_broad_ss$prevalence)/100
pcos_consortium_ss$prevalence<- as.numeric(pcos_consortium_ss$prevalence)/100

pcos_day$prevalence<- as.numeric(pcos_day$prevalence)/100
pcos_adj_age$prevalence<- as.numeric(pcos_adj_age$prevalence)/100
pcos_adj_age_bmi$prevalence<- as.numeric(pcos_adj_age_bmi$prevalence)/100


################################################
#Let's get the genome-wide significant variants#
################################################

pcos_ss_gw <- pcos_ss[which(as.numeric(pcos_ss$p_value) < 5e-08),]
pcos_broad_ss_gw <- pcos_broad_ss[which(as.numeric(pcos_broad_ss$p_value) < 5e-08),]
pcos_consortium_ss_gw <- pcos_consortium_ss[which(as.numeric(pcos_consortium_ss$p_value) < 5e-08),]
pcos_day_gw <- pcos_day[which(as.numeric(pcos_day$p_value) < 5e-08),] #0 - we will use the IVs utilized by Liu et al
pcos_adj_age_gw <- pcos_adj_age[which(as.numeric(pcos_adj_age$p_value) < 5e-08),]
pcos_adj_age_bmi_gw <- pcos_adj_age_bmi[which(as.numeric(pcos_adj_age_bmi$p_value) < 5e-08),]

#Let's format the data - I think it should be OK with how the data is formatted:

pcos_ss_gw_4_clump <- formatting_data_4_clumping(pcos_ss_gw)
pcos_broad_ss_gw_4_clump <- formatting_data_4_clumping(pcos_broad_ss_gw)
pcos_consortium_ss_gw_4_clump <- formatting_data_4_clumping(pcos_consortium_ss_gw)

pcos_adj_age_gw_4_clump <- formatting_data_4_clumping(pcos_adj_age_gw)
pcos_adj_age_bmi_gw_4_clump <- formatting_data_4_clumping(pcos_adj_age_bmi_gw)

#Let's order by p-value and make sure we are not including triallelic variants:

pcos_ss_gw_4_clump <- pcos_ss_gw_4_clump[order(pcos_ss_gw_4_clump$pval),]
pcos_broad_ss_gw_4_clump <- pcos_broad_ss_gw_4_clump[order(pcos_broad_ss_gw_4_clump$pval),]
pcos_consortium_ss_gw_4_clump <- pcos_consortium_ss_gw_4_clump[order(pcos_consortium_ss_gw_4_clump$pval),]
pcos_adj_age_gw_4_clump <- pcos_adj_age_gw_4_clump[order(pcos_adj_age_gw_4_clump$pval),]
pcos_adj_age_bmi_gw_4_clump <- pcos_adj_age_bmi_gw_4_clump[order(pcos_adj_age_bmi_gw_4_clump$pval),]

#Removing duplicates:

pcos_ss_gw_4_clump <- pcos_ss_gw_4_clump[which(duplicated(pcos_ss_gw_4_clump$chr_pos) == FALSE),]
pcos_broad_ss_gw_4_clump <- pcos_broad_ss_gw_4_clump[which(duplicated(pcos_broad_ss_gw_4_clump$chr_pos) == FALSE),]
pcos_consortium_ss_gw_4_clump <- pcos_consortium_ss_gw_4_clump[which(duplicated(pcos_consortium_ss_gw_4_clump$chr_pos) == FALSE),]
pcos_adj_age_gw_4_clump <- pcos_adj_age_gw_4_clump[which(duplicated(pcos_adj_age_gw_4_clump$chr_pos) == FALSE),]
pcos_adj_age_bmi_gw_4_clump <- pcos_adj_age_bmi_gw_4_clump[which(duplicated(pcos_ss_gw_4_clump$chr_pos) == FALSE),]

###########################
#Let's save all the data!!#
###########################

dir.create("output/2_replicating_liu_et_al/3_reverse_mr/0_ld_clumping/data_4_clumping")
dir.create("output/2_replicating_liu_et_al/3_reverse_mr/0_ld_clumping/clumped_data")
dir.create("output/2_replicating_liu_et_al/3_reverse_mr/0_ld_clumping/qsub_files")

write.table(pcos_ss_gw_4_clump, "output/2_replicating_liu_et_al/3_reverse_mr/0_ld_clumping/data_4_clumping/pcos_finngen_4_clumping.txt", quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
write.table(pcos_broad_ss_gw_4_clump, "output/2_replicating_liu_et_al/3_reverse_mr/0_ld_clumping/data_4_clumping/pcos_broad_4_clumping.txt", quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
write.table(pcos_consortium_ss_gw_4_clump, "output/2_replicating_liu_et_al/3_reverse_mr/0_ld_clumping/data_4_clumping/pcos_consortium_4_clumping.txt", quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
write.table(pcos_adj_age_gw_4_clump, "output/2_replicating_liu_et_al/3_reverse_mr/0_ld_clumping/data_4_clumping/pcos_tyrmi_4_clumping.txt", quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
write.table(pcos_adj_age_bmi_gw_4_clump, "output/2_replicating_liu_et_al/3_reverse_mr/0_ld_clumping/data_4_clumping/pcosadjbmi_tyrmi_4_clumping.txt", quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

#CHECKED!! Associations match with og sumstats


