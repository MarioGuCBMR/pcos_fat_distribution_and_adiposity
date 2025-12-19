##############
#INTRODUCTION#
##############

#This is just a collection of functions to run unweighted GRS and MR analysis.

###################
#Loading functions#
###################

rsid_cleaner <- function(rsid_1k){
  #This function takes an rsID from 1000G and retrieves only the rsID:
  
  rsid_ <- strsplit(rsid_1k, ":")[[1]][1]
  
  return(rsid_)
  
}

position_cleaner <- function(rsid_1k){
  #This function takes an rsID from 1000G and retrieves only the positions:
  
  pos_ <- strsplit(rsid_1k, ":")[[1]][2]
  
  return(pos_)
  
}

ref_cleaner <- function(rsid_1k){
  #This function takes an rsID from 1000G and retrieves only the reference allele:
  
  ref_ <- strsplit(rsid_1k, ":")[[1]][3]
  
  return(ref_)
  
}

alt_cleaner <- function(rsid_1k){
  #This function takes an rsID from 1000G and retrieves only the alternative allele:
  
  alt_ <- strsplit(rsid_1k, ":")[[1]][4]
  
  return(alt_)
  
}

proxies_cleaner <- function(prox_str){
  #This function transforms the string in which the proxies are delivered in
  #and transforms it into a vector of clean proxies that can go through the rest of the 
  #cleaner functions above.
  
  proxies_vect_tmp_1 <- strsplit(prox_str, "[,]")[[1]] #removes the comma
  proxies_vect_tmp_2 <- as.character(unlist(sapply(proxies_vect_tmp_1, strsplit, "[(]")))
  proxies_vect <- proxies_vect_tmp_2[which(str_detect(proxies_vect_tmp_2, "[)]") == FALSE)]
  
  return(proxies_vect)
  
}


proxies_obtainer <- function(proxy_df){
  #This function cleans the proxy data and returns a dataframe with rsID, CHR, BP and chr_pos.
  
  for(i in seq(1, length(proxy_df$SNP))){
    
    snp_of_ref <- proxy_df$SNP[i]
    chr_of_ref <- proxy_df$CHR[i]
    
    proxies <- proxy_df$SP2[i]
    proxies_clean <- proxies_cleaner(proxies)
    
    #Now that we have the data, we can obtain the rest:
    
    pos_proxies <- as.numeric(as.character(sapply(proxies_clean, position_cleaner)))
    rsid_proxies <- as.character(sapply(proxies_clean, rsid_cleaner))
    
    #Now let's make columns for the data that we want:
    
    snp_of_ref_clean <- rsid_cleaner(snp_of_ref)
    pos_of_ref_clean <- position_cleaner(snp_of_ref)
    
    #Let's check this out:
    
    lead_snp_clean<- rep(snp_of_ref_clean, length(proxies_clean))
    position_snp_clean<- rep(pos_of_ref_clean, length(proxies_clean))
    chr_clean<- rep(chr_of_ref, length(proxies_clean))
    chr_pos_clean <- paste("chr", chr_clean, ":", pos_proxies, sep = "")
    
    proxies_df <- cbind(proxies_clean, rsid_proxies, as.numeric(chr_clean), as.numeric(pos_proxies), chr_pos_clean, lead_snp_clean, as.numeric(position_snp_clean))
    
    colnames(proxies_df) <- c("ID", "variant", "chromosome", "base_pair_location", "chr_pos", "lead_variant", "lead_variant_position")
    
    if(i==1){
      
      final_proxies_df <- proxies_df
      
    } else {
      
      final_proxies_df <- rbind(final_proxies_df, proxies_df)
      
    }
    
  }
  
  return(final_proxies_df)
  
  
}


smoking_aligner <- function(query_ss, other_ss){
  
  #We are gonna put the example so that we can know what happens:
  
  #fiadjbmi_ss <- exp_df_found_pos
  #other_ss <- hdl_005
  
  #STEP 0: let's run the matching with TwoSampleMR so we need the data:
  
  #Let's first check if we have the effect_allele_frequency column:
  
  check <- which(colnames(query_ss) == "effect_allele_frequency")
  
  if(is_empty(check)){
    
    query_ss$effect_allele_frequency <- NA
    
  }
  
  #Now we can proceed
  
  exposure <- query_ss %>%
    select(chr_pos, chromosome, base_pair_location, effect_allele, other_allele, effect_allele_frequency, beta, standard_error, p_value, sample_size, variant)
  
  colnames(exposure) <- c("SNP", "chr.exposure", "pos.exposure", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure", "samplesize.exposure", "rsid.exposure")
  
  exposure$exposure <- "smoking"
  exposure$id.exposure <- "smoking"
  
  #Now with the outcome:
  
  check <- which(colnames(other_ss) == "effect_allele_frequency")
  
  if(is_empty(check)){
    
    other_ss$effect_allele_frequency <- NA
    
  }
  
  #Now we can proceed
  
  outcome <- other_ss %>%
    select(chr_pos, chromosome, base_pair_location, effect_allele, other_allele, effect_allele_frequency, beta, standard_error, p_value, sample_size, variant)
  
  colnames(outcome) <- c("SNP", "chr.outcome", "pos.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome", "rsid.outcome")
  
  outcome$outcome <- "idps"
  outcome$id.outcome <- "idps"
  
  #STEP 1: match the data. This will probably fail with tri-allelic SNPs. Here I think we are gonna be OK...
  
  merged_df <- harmonise_data(exposure, outcome, action=3)
  
  merged_df <- merged_df[which(merged_df$remove == FALSE),] #removing incompatible alleles
  merged_df <- merged_df[which(merged_df$palindromic == FALSE & merged_df$ambiguous == FALSE | merged_df$palindromic == TRUE & merged_df$ambiguous == FALSE),]
  
  #merged_df <- merged_df[which(merged_df$mr_keep),] #removing incompatible alleles
  
  #I checked that all is working great. FIadjBMI betas is positive. The rest is great.
  
  #STEP 3: reorganize the dataframe, cuz we need something clean:
  
  other_ss_aligned <- merged_df %>%
    select("rsid.outcome", "chr.outcome", "pos.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome", "SNP")
  
  colnames(other_ss_aligned) <- c("variant", "chromosome", "base_pair_location",  "effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value", "sample_size", "chr_pos")
  
  return(other_ss_aligned)
  
}



lower_ci <- function(beta_, se_){
  #A function that takes the beta and se from summary statistics
  #and generates the lower CI.
  
  lower_ci <- beta_ - qnorm(0.975)*se_
  
  return(lower_ci)
  
}

upper_ci <- function(beta_, se_){
  #A function that takes the beta and se from summary statistics
  #and generates the upper CI.
  
  upper_ci <- beta_ + qnorm(0.975)*se_
  
  return(upper_ci)
  
}



proxy_obtainer_and_aligner_exposure <- function(exp_df, out_df, proxy_df, exp_full_ss){
  
  #We are gonna put the example so that we can know what happens:
  
  #fiadjbmi_ss <- exp_df_found_pos
  #other_ss <- hdl_005
  
  #STEP 0: let's run the matching with TwoSampleMR so we need the data:
  
  #Let's first check if we have the effect_allele_frequency column:
  
  check <- which(colnames(exp_df) == "effect_allele_frequency")
  
  if(is_empty(check)){
    
    query_ss$effect_allele_frequency <- NA
    
  }
  
  #Now we can proceed
  
  exposure <- exp_df %>%
    select(chr_pos, chromosome, base_pair_location, effect_allele, other_allele, effect_allele_frequency, beta, standard_error, p_value, sample_size, variant)
  
  colnames(exposure) <- c("SNP", "chr.exposure", "pos.exposure", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure", "samplesize.exposure", "rsid.exposure")
  
  exposure$exposure <- "exposure"
  exposure$id.exposure <- "exposure"
  
  #Now with the outcome:
  
  check <- which(colnames(out_df) == "effect_allele_frequency")
  
  if(is_empty(check)){
    
    out_df$effect_allele_frequency <- NA
    
  }
  
  #Now we can proceed
  
  outcome <- out_df %>%
    select(chr_pos, chromosome, base_pair_location, effect_allele, other_allele, effect_allele_frequency, beta, standard_error, p_value, sample_size, variant)
  
  colnames(outcome) <- c("SNP", "chr.outcome", "pos.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome", "rsid.outcome")
  
  outcome$outcome <- "outcome"
  outcome$id.outcome <- "outcome"
  
  #STEP 1: match the data. This will probably fail with tri-allelic SNPs. Here I think we are gonna be OK...
  
  merged_df <- harmonise_data(exposure, outcome, action=2)
  
  merged_df <- merged_df[which(merged_df$remove == FALSE),] #removing incompatible alleles
  #merged_df <- merged_df[which(merged_df$mr_keep),] #removing incompatible alleles
  
  #I checked that all is working great. FIadjBMI betas is positive. The rest is great.
  
  #STEP 3: reorganize the dataframe, cuz we need something clean:
  
  other_ss_aligned <- merged_df %>%
    select("rsid.outcome", "chr.outcome", "pos.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome", "SNP")
  
  colnames(other_ss_aligned) <- c("variant", "chromosome", "base_pair_location",  "effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value", "sample_size", "chr_pos")
  
  #STEP 4: which are the variants that are not found??
  
  exp_df_missing <- exp_df[which(!(exp_df$chr_pos%in%other_ss_aligned$chr_pos)),]
  
  #STEP 5: check proxies for those SNPs:
  
  proxy_df_4_exp <- proxy_df[which(proxy_df$lead_variant_chr_pos%in%exp_df_missing$chr_pos),]
  
  #STEP 6: we are going to loop over the lead variants and go case by case to obtain the one with the lowest p-value:
  
  lead_variants_missing <- proxy_df_4_exp$lead_variant_chr_pos[which(duplicated(proxy_df_4_exp$lead_variant_chr_pos) == FALSE)]
  
  for(i in seq(1, length(lead_variants_missing))){
    
    proxy_tmp_df <- proxy_df_4_exp[which(proxy_df_4_exp$lead_variant_chr_pos == lead_variants_missing[i]),]
    
    out_in_tmp_df <- out_df[which(out_df$chr_pos%in%proxy_tmp_df$chr_pos),]
    
    if(length(out_in_tmp_df$chr_pos) == 0){
      
      print(i)
      print("No matches for this fella")
      print(lead_variants_missing[i])
      
      next()
      
    } 
    
    final_out_in_tmp_df <- out_in_tmp_df[which.min(as.numeric(out_in_tmp_df$p_value)),]
    final_out_in_tmp_df <- final_out_in_tmp_df[which.max(as.numeric(final_out_in_tmp_df$sample_size)),] #to avoid issues with P that are equal. Rare but might happen
    final_out_in_tmp_df <- final_out_in_tmp_df[which(duplicated(final_out_in_tmp_df$chr_pos) == FALSE),] #worst case scenario. Extremely rare. At this point we do not care which proxy we obtain.
    
    if(i==1){
      
      final_out_proxy <- final_out_in_tmp_df
      
    }
    
    if(i != 1 & exists("final_out_proxy") == FALSE){
      
      final_out_proxy <- final_out_in_tmp_df
      
    }
    
    if(i != 1 & exists("final_out_proxy") == TRUE){
      
      final_out_proxy <- rbind(final_out_proxy, final_out_in_tmp_df)
      
    }
    
  }
  
  #STEP 7: we are going to find the variants in the FIadjBMI, trait and align them to the increasing-risk allele. Then we will merge them to the exp_df original and run the other function that only makes the merge with the original out_df.
  
  exp_df_found <- exp_full_ss[which(exp_full_ss$chr_pos%in%final_out_proxy$chr_pos),]
  
  #Let's align this data to the risk-increasing allele:
  
  new_a1 <- ifelse(as.numeric(exp_df_found$beta) < 0, exp_df_found$other_allele, exp_df_found$effect_allele)
  new_a2 <- ifelse(as.numeric(exp_df_found$beta) < 0, exp_df_found$effect_allele, exp_df_found$other_allele)
  new_beta <- ifelse(as.numeric(exp_df_found$beta) < 0, as.numeric(exp_df_found$beta)*(-1), as.numeric(exp_df_found$beta))
  new_eaf <- ifelse(as.numeric(exp_df_found$beta) < 0, 1-as.numeric(exp_df_found$effect_allele_frequency), as.numeric(exp_df_found$effect_allele_frequency))
  
  exp_df_found$final_effect_allele <- new_a1
  exp_df_found$final_other_allele <- new_a2
  exp_df_found$final_beta <- new_beta
  exp_df_found$final_effect_allele_frequency <- new_eaf
  
  head(exp_df_found)
  
  summary(exp_df_found$final_beta) #all positive! Great.
  
  exp_df_found_pos <- exp_df_found %>%
    select("variant", "chromosome","base_pair_location", "final_effect_allele", "final_other_allele", "final_effect_allele_frequency", "final_beta", "standard_error", "p_value", "sample_size", "chr_pos")
  
  colnames(exp_df_found_pos) <- c("variant", "chromosome","base_pair_location", "effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value", "sample_size", "chr_pos")
  
  exp_df <- rbind(exp_df, exp_df_found_pos)
  
  #Finally, let's use that to obtain the final set of data:
  
  final_merge <- fiadjbmi_aligner(exp_df, out_df)#269!! with the ones we recoveredd!!
  
  exp_df_end <- exp_df[which(exp_df$chr_pos%in%final_merge$chr_pos),]
  
  return(exp_df_end)
  
}


proxy_obtainer_and_aligner_outcome <- function(exp_df, out_df, proxy_df, exp_full_ss){
  
  #We are gonna put the example so that we can know what happens:
  
  #fiadjbmi_ss <- exp_df_found_pos
  #other_ss <- hdl_005
  
  #STEP 0: let's run the matching with TwoSampleMR so we need the data:
  
  #Let's first check if we have the effect_allele_frequency column:
  
  check <- which(colnames(exp_df) == "effect_allele_frequency")
  
  if(is_empty(check)){
    
    query_ss$effect_allele_frequency <- NA
    
  }
  
  #Now we can proceed
  
  exposure <- exp_df %>%
    select(chr_pos, chromosome, base_pair_location, effect_allele, other_allele, effect_allele_frequency, beta, standard_error, p_value, sample_size, variant)
  
  colnames(exposure) <- c("SNP", "chr.exposure", "pos.exposure", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure", "samplesize.exposure", "rsid.exposure")
  
  exposure$exposure <- "exposure"
  exposure$id.exposure <- "exposure"
  
  #Now with the outcome:
  
  check <- which(colnames(out_df) == "effect_allele_frequency")
  
  if(is_empty(check)){
    
    out_df$effect_allele_frequency <- NA
    
  }
  
  #Now we can proceed
  
  outcome <- out_df %>%
    select(chr_pos, chromosome, base_pair_location, effect_allele, other_allele, effect_allele_frequency, beta, standard_error, p_value, sample_size, variant)
  
  colnames(outcome) <- c("SNP", "chr.outcome", "pos.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome", "rsid.outcome")
  
  outcome$outcome <- "outcome"
  outcome$id.outcome <- "outcome"
  
  #STEP 1: match the data. This will probably fail with tri-allelic SNPs. Here I think we are gonna be OK...
  
  merged_df <- harmonise_data(exposure, outcome, action=2)
  
  merged_df <- merged_df[which(merged_df$remove == FALSE),] #removing incompatible alleles
  #merged_df <- merged_df[which(merged_df$mr_keep),] #removing incompatible alleles
  
  #I checked that all is working great. FIadjBMI betas is positive. The rest is great.
  
  #STEP 3: reorganize the dataframe, cuz we need something clean:
  
  other_ss_aligned <- merged_df %>%
    select("rsid.outcome", "chr.outcome", "pos.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome", "SNP")
  
  colnames(other_ss_aligned) <- c("variant", "chromosome", "base_pair_location",  "effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value", "sample_size", "chr_pos")
  
  #STEP 4: which are the variants that are not found??
  
  exp_df_missing <- exp_df[which(!(exp_df$chr_pos%in%other_ss_aligned$chr_pos)),]
  
  #STEP 5: check proxies for those SNPs:
  
  proxy_df_4_exp <- proxy_df[which(proxy_df$lead_variant_chr_pos%in%exp_df_missing$chr_pos),]
  
  #STEP 6: we are going to loop over the lead variants and go case by case to obtain the one with the lowest p-value:
  
  lead_variants_missing <- proxy_df_4_exp$lead_variant_chr_pos[which(duplicated(proxy_df_4_exp$lead_variant_chr_pos) == FALSE)]
  
  for(i in seq(1, length(lead_variants_missing))){
    
    proxy_tmp_df <- proxy_df_4_exp[which(proxy_df_4_exp$lead_variant_chr_pos == lead_variants_missing[i]),]
    
    out_in_tmp_df <- out_df[which(out_df$chr_pos%in%proxy_tmp_df$chr_pos),]
    
    if(length(out_in_tmp_df$chr_pos) == 0){
      
      print(i)
      print("No matches for this fella")
      print(lead_variants_missing[i])
      
      next()
      
    } 
    
    final_out_in_tmp_df <- out_in_tmp_df[which.min(as.numeric(out_in_tmp_df$p_value)),]
    final_out_in_tmp_df <- final_out_in_tmp_df[which.max(as.numeric(final_out_in_tmp_df$sample_size)),] #to avoid issues with P that are equal. Rare but might happen
    final_out_in_tmp_df <- final_out_in_tmp_df[which(duplicated(final_out_in_tmp_df$chr_pos) == FALSE),] #worst case scenario. Extremely rare. At this point we do not care which proxy we obtain.
    
    if(i==1){
      
      final_out_proxy <- final_out_in_tmp_df
      
    }
    
    if(i != 1 & exists("final_out_proxy") == FALSE){
      
      final_out_proxy <- final_out_in_tmp_df
      
    }
    
    if(i != 1 & exists("final_out_proxy") == TRUE){
      
      final_out_proxy <- rbind(final_out_proxy, final_out_in_tmp_df)
      
    }
    
  }
  
  #STEP 7: we are going to find the variants in the FIadjBMI, trait and align them to the increasing-risk allele. Then we will merge them to the exp_df original and run the other function that only makes the merge with the original out_df.
  
  exp_df_found <- exp_full_ss[which(exp_full_ss$chr_pos%in%final_out_proxy$chr_pos),]
  
  #Let's align this data to the risk-increasing allele:
  
  new_a1 <- ifelse(as.numeric(exp_df_found$beta) < 0, exp_df_found$other_allele, exp_df_found$effect_allele)
  new_a2 <- ifelse(as.numeric(exp_df_found$beta) < 0, exp_df_found$effect_allele, exp_df_found$other_allele)
  new_beta <- ifelse(as.numeric(exp_df_found$beta) < 0, as.numeric(exp_df_found$beta)*(-1), as.numeric(exp_df_found$beta))
  new_eaf <- ifelse(as.numeric(exp_df_found$beta) < 0, 1-as.numeric(exp_df_found$effect_allele_frequency), as.numeric(exp_df_found$effect_allele_frequency))
  
  exp_df_found$final_effect_allele <- new_a1
  exp_df_found$final_other_allele <- new_a2
  exp_df_found$final_beta <- new_beta
  exp_df_found$final_effect_allele_frequency <- new_eaf
  
  head(exp_df_found)
  
  summary(exp_df_found$final_beta) #all positive! Great.
  
  exp_df_found_pos <- exp_df_found %>%
    select("variant", "chromosome","base_pair_location", "final_effect_allele", "final_other_allele", "final_effect_allele_frequency", "final_beta", "standard_error", "p_value", "sample_size", "chr_pos")
  
  colnames(exp_df_found_pos) <- c("variant", "chromosome","base_pair_location", "effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value", "sample_size", "chr_pos")
  
  #Let's merge them with the exp_df. Since those that are missing are not in the out_df, there is no need to worry about including them here. The hamonization will properly remove them:
  
  exp_df <- rbind(exp_df, exp_df_found_pos)
  
  #Finally, let's use that to obtain the final set of data:
  
  final_merge <- fiadjbmi_aligner(exp_df, out_df)#269!! with the ones we recoveredd!!
  
  return(final_merge)
  
}


grs_producer <- function(query_ss, list_of_traits_, type_of_trait_, trait_names_, weighted=FALSE, proxies, exp_full_ss){
  
  #We are gonna put the example so that we can know what happens:
  
  for(index_trait in seq(1, length(list_of_traits_))){
    
    print(trait_names_[index_trait])
    
    #Get the data that aligns for our trait:
    
    aligned_outcome <- fiadjbmi_aligner(query_ss, as.data.frame(list_of_traits_[index_trait]))
    
    #WAIT: let's check if all variant are there:
    
    check <- length(which(query_ss$chr_pos%in%aligned_outcome$chr_pos))
    
    #Match and align just in case something went wrong:
    
    query_match <- query_ss[which(query_ss$chr_pos%in%aligned_outcome$chr_pos),]
    aligned_outcome <- aligned_outcome[order(match(aligned_outcome$chr_pos, query_match$chr_pos)),]
    
    #There is the possibility there is a mismatch due to duplicated variants
    #because some GWAS ss are shit (looking at you FI from Lagou).
    
    aligned_outcome <- aligned_outcome[order(aligned_outcome$sample_size, decreasing = TRUE),] #we will remove the duplicates according to sample size. Those with lower, go out
    aligned_outcome <- aligned_outcome[which(duplicated(aligned_outcome$chr_pos) == FALSE),]
    aligned_outcome <- aligned_outcome[order(match(aligned_outcome$chr_pos, query_match$chr_pos)),] #we order again so that we have no issues.
    
    if(check < length(query_ss$chr_pos)){ #if not all variants are there
      
      query_ss_ <- proxy_obtainer_and_aligner_exposure(exp_df = query_ss, out_df = as.data.frame(list_of_traits_[index_trait]), proxy_df= proxies, exp_full_ss = exp_full_ss)
      aligned_outcome <- proxy_obtainer_and_aligner_outcome(exp_df = query_ss, out_df = as.data.frame(list_of_traits_[index_trait]), proxy_df= proxies, exp_full_ss = exp_full_ss)
      query_match <- query_ss_[which(query_ss_$chr_pos%in%aligned_outcome$chr_pos),]
      aligned_outcome <- aligned_outcome[order(match(aligned_outcome$chr_pos, query_match$chr_pos)),]
      
      #There is the possibility there is a mismatch due to duplicated variants
      #because some GWAS ss are shit (looking at you FI from Lagou).
      
      aligned_outcome <- aligned_outcome[order(aligned_outcome$sample_size, decreasing = TRUE),] #we will remove the duplicates according to sample size. Those with lower, go out
      aligned_outcome <- aligned_outcome[which(duplicated(aligned_outcome$chr_pos) == FALSE),]
      aligned_outcome <- aligned_outcome[order(match(aligned_outcome$chr_pos, query_match$chr_pos)),] #we order again so that we have no issues.
      
    }
    
    #Finally produce the genetic risk scores:
    
    if(weighted == TRUE){
      
      grs_res <- gtx::grs.summary(query_match$beta, aligned_outcome$beta, aligned_outcome$standard_error, aligned_outcome$sample_size)
      
    } else {
      
      weights <- rep(1, length(aligned_outcome$beta))
      
      grs_res <- gtx::grs.summary(weights, aligned_outcome$beta, aligned_outcome$standard_error, aligned_outcome$sample_size)
      
    }
    
    #Let's make the results for the GRS parseable:
    
    #Finally, merge the results:
    
    if(index_trait == 1){
      
      beta_ <- grs_res$ahat
      se_ <- grs_res$aSE
      lower_ci_ <- lower_ci(beta_, se_)
      upper_ci_ <- upper_ci(beta_, se_)
      pval_ <- grs_res$pval
      trait_ <- trait_names_[index_trait]
      type_ <- type_of_trait_[index_trait]
      nsnps_ <- length(aligned_outcome$chr_pos)
      
      grs_vect <- c(beta_, se_, lower_ci_, upper_ci_, pval_, trait_, type_, nsnps_)
      grs_end <- as.data.frame(t(grs_vect))
      colnames(grs_end) <- c("beta", "se", "lower_ci", "upper_ci", "pval", "trait", "type", "nsnps")
      
    } else {
      
      beta_ <- grs_res$ahat
      se_ <- grs_res$aSE
      lower_ci_ <- lower_ci(beta_, se_)
      upper_ci_ <- upper_ci(beta_, se_)
      pval_ <- grs_res$pval
      trait_ <- trait_names_[index_trait]
      type_ <- type_of_trait_[index_trait]
      nsnps_ <- length(aligned_outcome$chr_pos)
      
      grs_vect <- c(beta_, se_, lower_ci_, upper_ci_, pval_, trait_, type_, nsnps_)
      grs_end <- rbind(grs_end, grs_vect)
      
    }
    
  }
  
  #Once everything has been looped we can just return the dataframe
  
  return(grs_end)
  
}
