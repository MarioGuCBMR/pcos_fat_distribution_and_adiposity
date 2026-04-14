##############
#INTRODUCTION#
##############

#This code will prepare the data to obtain proxies for all PCOS GWAS.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)
library(TwoSampleMR)

######################
#Loading clumped data#
######################

#Let's get a path where the clusters are so that we can loop throught them:

path_2_input <- "/projects/kilpelainen-AUDIT/people/zlc436/PCOS_2026/"

setwd(path_2_input)

pcos_finngen_clumped <- fread("output/2_replicating_liu_et_al/3_reverse_mr/0_ld_clumping/clumped_data/pcos_finngen_2026_clumped.txt") #7
pcos_broad_clumped <- fread("output/2_replicating_liu_et_al/3_reverse_mr/0_ld_clumping/clumped_data/pcos_broad_2026_clumped.txt") #12
pcos_consortium_clumped <- fread("output/2_replicating_liu_et_al/3_reverse_mr/0_ld_clumping/clumped_data/pcos_consortium_2026_clumped.txt") #5
pcos_day_clumped <- fread("output/1_curated_data/pcos_day_ivs.txt") #14
pcos_venkatesh_clumped <- fread("output/2_replicating_liu_et_al/3_reverse_mr/0_ld_clumping/clumped_data/pcos_venkatesh_2026_clumped.txt") #17
pcos_tyrmi_clumped <- fread("output/2_replicating_liu_et_al/3_reverse_mr/0_ld_clumping/clumped_data/pcos_tyrmi_2026_clumped.txt") #6
pcosadjbmi_tyrmi_clumped <- fread("output/2_replicating_liu_et_al/3_reverse_mr/0_ld_clumping/clumped_data/pcosadjbmi_tyrmi_2026_clumped.txt") #3

#We are going to get the proxies separately to get lead SNPs independently, but we will need to be careful with this - OK?

dir.create("output/2_replicating_liu_et_al/3_reverse_mr/0_proxies/")

fwrite(as.data.frame(pcos_finngen_clumped$rsid), "output/2_replicating_liu_et_al/3_reverse_mr/0_proxies/proxies_4_pcos_finngen.txt", quote=FALSE, row.names = FALSE, col.names=FALSE, sep =" ")
fwrite(as.data.frame(pcos_broad_clumped$rsid), "output/2_replicating_liu_et_al/3_reverse_mr/0_proxies/proxies_4_pcos_broad.txt", quote=FALSE, row.names = FALSE, col.names=FALSE, sep =" ")
fwrite(as.data.frame(pcos_consortium_clumped$rsid), "output/2_replicating_liu_et_al/3_reverse_mr/0_proxies/proxies_4_pcos_consortium.txt", quote=FALSE, row.names = FALSE, col.names=FALSE, sep =" ")
fwrite(as.data.frame(pcos_day_clumped$variant), "output/2_replicating_liu_et_al/3_reverse_mr/0_proxies/proxies_4_pcos_day.txt", quote=FALSE, row.names = FALSE, col.names=FALSE, sep =" ")
fwrite(as.data.frame(pcos_venkatesh_clumped$rsid), "output/2_replicating_liu_et_al/3_reverse_mr/0_proxies/proxies_4_pcos_venkatesh.txt", quote=FALSE, row.names = FALSE, col.names=FALSE, sep =" ")
fwrite(as.data.frame(pcos_tyrmi_clumped$rsid), "output/2_replicating_liu_et_al/3_reverse_mr/0_proxies/proxies_4_pcos_tyrmi.txt", quote=FALSE, row.names = FALSE, col.names=FALSE, sep =" ")
fwrite(as.data.frame(pcosadjbmi_tyrmi_clumped$rsid), "output/2_replicating_liu_et_al/3_reverse_mr/0_proxies/proxies_4_pcosadjbmi_tyrmi.txt", quote=FALSE, row.names = FALSE, col.names=FALSE, sep =" ")

######################################################
#Let's load the data and prepare the proxies properly#
######################################################

pcos_finngen_proxies <- fread("output/2_replicating_liu_et_al/3_reverse_mr/0_proxies/pcos_finngen_haploreg_output.txt", fill = TRUE)
pcos_broad_proxies <- fread("output/2_replicating_liu_et_al/3_reverse_mr/0_proxies/pcos_broad_haploreg_output.txt", fill = TRUE)
pcos_consortium_proxies <- fread("output/2_replicating_liu_et_al/3_reverse_mr/0_proxies/pcos_consortium_haploreg_output.txt", fill = TRUE)
pcos_day_proxies <- fread("output/2_replicating_liu_et_al/3_reverse_mr/0_proxies/pcos_day_haploreg_output.txt", fill = TRUE)
pcos_venkatesh_proxies <- fread("output/2_replicating_liu_et_al/3_reverse_mr/0_proxies/pcos_venkatesh_haploreg_output.txt", fill = TRUE)
pcos_tyrmi_proxies <- fread("output/2_replicating_liu_et_al/3_reverse_mr/0_proxies/pcos_tyrmi_haploreg_output.txt", fill = TRUE)
pcosadjbmi_tyrmi_proxies <- fread("output/2_replicating_liu_et_al/3_reverse_mr/0_proxies/pcosadjbmi_tyrmi_haploreg_output.txt", fill = TRUE)

#Let's check the data:

print(length(unique(pcos_finngen_proxies$query_snp_rsid))) #7
print(length(unique(pcos_broad_proxies$query_snp_rsid))) #12
print(length(unique(pcos_consortium_proxies$query_snp_rsid))) #5

print(length(unique(pcos_day_proxies$query_snp_rsid))) #15 -one of them is fake
print(length(unique(pcos_venkatesh_proxies$query_snp_rsid))) #18 -one of them is fake
print(length(unique(pcos_tyrmi_proxies$query_snp_rsid))) #6
print(length(unique(pcosadjbmi_tyrmi_proxies$query_snp_rsid))) #3

pcos_day_proxies=pcos_day_proxies[which(pcos_day_proxies$query_snp_rsid!=""),] #14
pcos_venkatesh_proxies=pcos_venkatesh_proxies[which(pcos_venkatesh_proxies$query_snp_rsid!=""),] #17

#Let's add an "analysis pipeline":

pcos_finngen_proxies$analysis <- "PCOS (FinnGen)"
pcos_broad_proxies$analysis <- "PCOS (Broad)"
pcos_consortium_proxies$analysis <- "PCOS (Consortium)"
pcos_day_proxies$analysis <- "PCOS (Day)"
pcos_venkatesh_proxies$analysis <- "PCOS (Venkatesh)"
pcos_tyrmi_proxies$analysis <- "PCOS (Tyrmi)"
pcosadjbmi_tyrmi_proxies$analysis <- "PCOSadjBMI (Tyrmi)"

#Let's add the data:

proxies <- rbind(pcos_finngen_proxies, pcos_broad_proxies, pcos_consortium_proxies, pcos_day_proxies, pcos_venkatesh_proxies, pcos_tyrmi_proxies, pcosadjbmi_tyrmi_proxies)

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

fwrite(proxies, "output/2_replicating_liu_et_al/3_reverse_mr/0_proxies/all_proxies_dictionary_build_37_38.txt")


