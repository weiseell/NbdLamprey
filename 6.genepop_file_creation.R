#Converting .vcf genotypes to GENEPOP format
#Created by: Ellie Weise
#Originally Created on: Oct 25th, 2019
#Last Edited on: Nov 4th, 2019
#Goals:
#1. convert selected SNPs from vcf to COLONY format
#2. subset formatted file into cohorts and create COLONY files
#libraries
library(tidyverse)
#homebrew functions

#load in data
load(file = "Input/gt_filtered.rda")
load(file = "Input/tag_selected_SNPs.rda")

#selecting genotypes that match selected SNPs
gt <- merge(SNP_selected,geno)
gt <- gt %>% 
  select(-het:-MAF) %>% 
  mutate(SNP_name = paste0("SNP",1:nrow(gt))) %>% 
  select(-CHROM:-target) %>% 
  select(SNP_name,PM_OCQ_001:PM_BMR_1089) %>% 
  gather(key = "Indiv",value = gt,-SNP_name)

genepop <- vcf_genepop(vcf = vcf)









