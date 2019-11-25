#Quality filtering for NeEstimator friendly SNP set
#Created by: Ellie Weise
#Originally Created on: Nov 5th, 2019
#Last Edited on: Nov 5th, 2019

#Goals:
#Create a SNP set with the following additional filters:
#genotype coverage of over 80%
#target loci only
#one SNP per tag

#libraries
library(tidyverse)

#homebrew functions
source("Homebrew/tag_filter.R")

#load in data
load(file = "Input/gt_summary_with_targets.rda")

#filtering for target and gt coverage
geno2 <- geno1 %>% 
  filter(target != "NonTarget") %>% 
  filter(MAF > 0.05)
#selecting SNP set
SNP_selected <- tag_filter(df = geno2)

#saving SNP set
save(SNP_selected,file = "Input/tag_selected_SNPs.rda")

