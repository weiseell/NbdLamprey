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
source("Homebrew/vcf_genepop.R")

#load in data
load(file = "Input/gt_filtered.rda")
load(file = "Input/tag_selected_SNPs.rda")
df <- read.table(file = "Input/exp_lengths_weights.txt",header = T,sep = "\t",stringsAsFactors = F)
#selecting genotypes that match selected SNPs
gt <- merge(SNP_selected,geno)
gt1 <- gt %>% 
  select(-het:-MAF) %>% 
  mutate(SNP_name = paste0("SNP",1:nrow(gt))) %>% 
  select(-CHROM:-target) %>% 
  select(SNP_name,PM_OCQ_001:PM_BMR_1089) %>% 
  gather(key = "Indiv",value = gt,-SNP_name)

genepop <- vcf_genepop(vcf = gt1)

#using cat to create a chromosome file and a GENEPOP file
#chromosome file for all locations  (using the same SNPs)
names <- gt %>% 
  select(-het:-MAF) %>% 
  mutate(SNP_name = paste0("SNP",1:nrow(gt))) %>% 
  select(CHROM,SNP_name)

write.table(names,file = "Output/scaf_SNP_key.txt",append = F, sep = " ",row.names = F,col.names = T,quote = F)
#subset for each location/cohort first
##BMR
df <- df %>% rename(ID = ID_indiv)
genepop <- genepop %>% rename(ID = Indiv)
genepop1 <- merge(genepop,df,on = ID)
bmr_gp1 <- genepop1 %>% 
  mutate(ID2 = ID) %>% 
  separate(ID2,into = c("spp","loc","num"),sep = "_") %>% 
  filter(loc == "BMR" & Year_collect != 2019) %>% 
  mutate(pop = "pop1") %>% 
  select(-spp:-num) %>% 
  select(pop,ID,SNP1:SNP999)

bmr_gp2 <- genepop1 %>% 
  mutate(ID2 = ID) %>% 
  separate(ID2,into = c("spp","loc","num"),sep = "_") %>% 
  filter(loc == "BMR" & Year_collect == 2019) %>% 
  mutate(pop = "pop2") %>% 
  select(-spp:-num) %>% 
  select(pop,ID,SNP1:SNP999)

bmr_gp <- rbind(bmr_gp1,bmr_gp2)
bmr_snps <- bmr_gp %>% select(-pop:-ID)

genepop_create(SNPs = bmr_snps,df = bmr_gp,output_file = "pm.bmr.genepop.txt",title = "Genotyped SNPs on Rapture tag panel for P.Marinus larvae - Black Mallard River")

#cheboygan
che_gp <- genepop1 %>% 
  mutate(ID2 = ID) %>% 
  separate(ID2,into = c("spp","loc","num"),sep = "_") %>% 
  filter(loc == "CHE") %>% 
  mutate(pop = "pop1") %>% 
  select(-spp:-num) %>% 
  select(pop,ID,SNP1:SNP999)
che_snps <- che_gp %>% select(-pop:-ID)
genepop_create(SNPs = che_snps,df = che_gp,output_file = "pm.che.genepop.txt",
               title = "Genotyped SNPs on Rapture tag panel for P.Marinus larvae - Cheboygan River")



