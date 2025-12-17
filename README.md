# Genetic predisposition of fat distribution is associated with PCOS risk 

This repository contains the code for reproducing the code for the study "".
In this code, genetic analyses like genetic correlations and mendelian randomization are performed to assess whether the genetic predisposition to abdominal fat accumuulation is associated with PCOS risk.
In particular, we take into consideration the effects of BMI as sensitivity analyses in several ways. 

## Downloading data

The analyses have been all performed utilizing GWAS summary statistics. Most GWAS summary statistics are female-specific, European GWAS. See Supplementary Table 1 in manuscript/extended_tables.xlsx for information on all GWAS utilized and links to download them. 

## Performing the analyses:

All analyses can be reproduced once all GWAS summary statistics are downloaded in the raw_data/ folder. Following the code in the order found in the code/ folder to run the analyses. 
All analyses were performed in R/4.3.1 in a HPC environment with some specific exceptions that are detailed below. 
