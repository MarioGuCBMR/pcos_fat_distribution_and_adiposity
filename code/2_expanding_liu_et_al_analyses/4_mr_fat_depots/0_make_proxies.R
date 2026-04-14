##############
#INTRODUCTION#
##############

#This code prepares the lead SNPs for all fat depots MRI phenotypes (female-specific) to run haploreg and retrieve proxies.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

#######################
#Loading exposure data#
#######################

#Let's get a path where the clusters are so that we can loop throught them:

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/PCOS_2026/")

asat <- fread("output/1_curated_data/asat_female_ivs.txt") #1
gsat <- fread("output/1_curated_data/gsat_female_ivs.txt") #5
vat <- fread("output/1_curated_data/vat_female_ivs.txt") #2

asatadj <- fread("output/1_curated_data/asatadj_female_ivs.txt") #3
gsatadj <- fread("output/1_curated_data/gsatadj_female_ivs.txt") #21
vatadj <- fread("output/1_curated_data/vatadj_female_ivs.txt") #10

asat_gsat_ratio <- fread("output/1_curated_data/asat_gsat_ratio_female_ivs.txt") #7
vat_gsat_ratio <- fread("output/1_curated_data/vat_gsat_ratio_female_ivs.txt") #8
vat_asat_ratio <- fread("output/1_curated_data/vat_asat_ratio_female_ivs.txt") #9

#####################
#Let's save the data#
#####################

dir.create("output/")
dir.create("output/2_replicating_liu_et_al/4_mr_fat_depots")
dir.create("output/2_replicating_liu_et_al/4_mr_fat_depots/0_proxies")

fwrite(as.data.frame(asat$variant), "output/2_replicating_liu_et_al/4_mr_fat_depots/0_proxies/asat_female_ivs_4_haploreg.txt", quote=FALSE, row.names = FALSE, col.names=FALSE, sep =" ")
fwrite(as.data.frame(gsat$variant), "output/2_replicating_liu_et_al/4_mr_fat_depots/0_proxies/gsat_female_ivs_4_haploreg.txt", quote=FALSE, row.names = FALSE, col.names=FALSE, sep =" ")
fwrite(as.data.frame(vat$variant), "output/2_replicating_liu_et_al/4_mr_fat_depots/0_proxies/vat_female_ivs_4_haploreg.txt", quote=FALSE, row.names = FALSE, col.names=FALSE, sep =" ")

fwrite(as.data.frame(asatadj$variant), "output/2_replicating_liu_et_al/4_mr_fat_depots/0_proxies/asatadj_female_ivs_4_haploreg.txt", quote=FALSE, row.names = FALSE, col.names=FALSE, sep =" ")
fwrite(as.data.frame(gsatadj$variant), "output/2_replicating_liu_et_al/4_mr_fat_depots/0_proxies/gsatadj_female_ivs_4_haploreg.txt", quote=FALSE, row.names = FALSE, col.names=FALSE, sep =" ")
fwrite(as.data.frame(vatadj$variant), "output/2_replicating_liu_et_al/4_mr_fat_depots/0_proxies/vatadj_female_ivs_4_haploreg.txt", quote=FALSE, row.names = FALSE, col.names=FALSE, sep =" ")

fwrite(as.data.frame(asat_gsat_ratio$variant), "output/2_replicating_liu_et_al/4_mr_fat_depots/0_proxies/asat_gsat_female_ivs_4_haploreg.txt", quote=FALSE, row.names = FALSE, col.names=FALSE, sep =" ")
fwrite(as.data.frame(vat_gsat_ratio$variant), "output/2_replicating_liu_et_al/4_mr_fat_depots/0_proxies/vat_gsat_ratio_female_ivs_4_haploreg.txt", quote=FALSE, row.names = FALSE, col.names=FALSE, sep =" ")
fwrite(as.data.frame(vat_asat_ratio$variant), "output/2_replicating_liu_et_al/4_mr_fat_depots/0_proxies/vat_asat_ratio_female_ivs_4_haploreg.txt", quote=FALSE, row.names = FALSE, col.names=FALSE, sep =" ")
