##############
#INTRODUCTION#
##############

#This code prepares the lead SNPs for BMI, WHR and WHRadjBMI (female-specific) to run haploreg and retrieve proxies.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

#######################
#Loading exposure data#
#######################

#Let's get a path where the clusters are so that we can loop throught them:

setwd("your_path")

bmi <- fread("output/1_curated_data/bmi_female_ivs.txt") #278
whr <- fread("output/1_curated_data/whr_female_ivs.txt") #196
whradjbmi <- fread("output/1_curated_data/whradjbmi_female_ivs.txt") #256

#####################
#Let's save the data#
#####################

dir.create("output//")
dir.create("output/2_replicating_liu_et_al/2_mr")
dir.create("output/2_replicating_liu_et_al/2_mr/0_proxies")

fwrite(as.data.frame(bmi$variant), "output/2_replicating_liu_et_al/2_mr/0_proxies/bmi_female_ivs_4_haploreg.txt", quote=FALSE, row.names = FALSE, col.names=FALSE, sep =" ")
fwrite(as.data.frame(whr$variant), "output/2_replicating_liu_et_al/2_mr/0_proxies/whr_female_ivs_4_haploreg.txt", quote=FALSE, row.names = FALSE, col.names=FALSE, sep =" ")
fwrite(as.data.frame(whradjbmi$variant), "output/2_replicating_liu_et_al/2_mr/0_proxies/whradjbmi_female_ivs_4_haploreg.txt", quote=FALSE, row.names = FALSE, col.names=FALSE, sep =" ")
