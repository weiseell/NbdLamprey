#gt_SNP_stats
#function to calculate heterozygosity and percent genotyped for each SNP
#also incorporates MAF data to quality filter (can turn this off)

#inputs:
#genotype file with individuals in columns, SNPs in rows
#MAF for each SNP as a vector
#total number of individuals in gt set
#heterozygosity filter - logical
gt_SNP_stats <- function(gt,n_indiv,hetfilter  = T){
  require(tidyverse)
  #getting heterozygosity and missing indiv sums
  gt$het_counts <- rowSums(gt == "0/1")
  gt$gt_missing <- rowSums(gt == "./.")
  #calculating heterozygosity and number of genotyped individuals for all loci
  #also adding MAF values at this  step
  out <- data.frame(ID = gt$ID,stringsAsFactors = F)
  out$het <- gt$het_counts/n_indiv
  out$pGT <- (1-(gt$gt_missing/n_indiv))
  #heterozygosity filter
  if(hetfilter == T){
    subset(out,out$het < 0.5)
  }
}

