#Prepping more input files
#Created by: Ellie Weise
#Originally Created on: Oct 25th, 2019
#Last Edited on: Nov 4th, 2019
#Goals:
#1. Calculate percent genotype coverage and heterozygosity for all SNPs
#2. Create a file with only SNP names and summary stats
#3. Match SNPs to Rapture panel

#load libraries
library(tidyverse)

#homebrew functions
source("Homebrew/gt_SNP_stats.R")
source("Homebrew/match_tags.R")

#load in data:
load(file = "Input/MAF_clean.rda")
load(file = "Input/gt_filtered.rda")
load(file = "Input/rapture_panel_all_SNPs.rda")

#Goal 1
#calculating %coverage, heterozygosity for all loci
#merging with MAF data
#should have a large table with all genotypes and three summary criteria

#Checking MAF file: if there are MAF > 0.5, fix so they're < 0.5
af$Min_AF2 <- af$Min_AF
for(i in 1:length(af$Maj_AF)){
  if (af$Min_AF[i] > 0.5) {
    af$Min_AF2[i] <- 1-af$Min_AF[i]
  }
}
max(af$Min_AF2)
MAF <- af$Min_AF2

geno1 <- gt_SNP_stats(gt = geno,MAF = MAF, n_indiv = 1536,hetfilter = T)

##matching acceptable loci with Rapture panel
#matching genotypes to rapture panel

#getting unique tags
unique_tags <- unique(rapture$rad_tag_name)
rapture1 <- as.data.frame(unique_tags)

#separating tag IDs to make location ranges on reference genome scaffolds for all tags
rapture1 <- rapture1 %>% 
  separate(unique_tags, into = c("CHROM","pos"),sep = ":") %>% 
  separate(pos, into = c("min","max"),sep = "-") %>% 
  mutate(ID = paste0(CHROM,":",min,"-",max)) %>% 
  select(ID, CHROM, min, max)
SNPs <- geno1 %>% select(CHROM,POS)

geno2 <- match_tags(SNPs = SNPs,tags = rapture1,target = T)


#save two main outputs from this script into the input folder
#summary stats and on target genotypes
save(geno1,file = "Input/gt_summary_all.rda")
save(geno2,file = "Input/gt_summary_targets.rda")