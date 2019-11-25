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
source("Homebrew/genepop_create.R")
#load in data
load(file = "Input/gt_filtered.rda")
load(file = "Input/tag_selected_SNPs.rda")
load(file = "Input/bmr_cohort_data.txt")
bmr <- bmr %>% mutate(ID = paste0(species,"_",loc,"_",num)) %>% select(ID,cohort)
df <- read.table(file = "Input/exp_lengths_weights.txt",header = T,sep = "\t",stringsAsFactors = F)
locs <- read.table("Input/locs_and_sublocations.txt",header = T,sep = "\t",stringsAsFactors = F)

#selecting genotypes that match selected SNPs
gt <- merge(SNP_selected,geno)
gt1 <- gt %>% 
  select(-het:-MAF) %>% 
  mutate(SNP = paste0("SNP",1:nrow(gt))) %>% 
  select(-CHROM:-target) %>% 
  select(SNP,PM_OCQ_001:PM_BMR_1089) %>% 
  gather(key = "ID",value = gt,-SNP)

genepop <- vcf_genepop(vcf = gt1)

#using cat to create a chromosome file and a GENEPOP file
#chromosome file for all locations  (using the same SNPs)
names <- gt %>% 
  select(-het:-MAF) %>% 
  mutate(SNP_name = paste0("SNP",1:nrow(gt))) %>% 
  select(CHROM,SNP_name)

write.table(names,file = "Output/scaf_SNP_key.txt",append = F, sep = " ",row.names = F,col.names = T,quote = F)
#subset for each location/cohort first
df <- df %>% rename(ID = ID_indiv)
genepop1 <- merge(genepop,df,on = ID)
colnames(bmr) <- c("ID","cohort")
genepop1 <- merge(genepop1,bmr)
##BMR - above river #####
bmr_gp1 <- genepop1 %>% 
  mutate(ID2 = ID) %>% 
  separate(ID2,into = c("spp","loc","num"),sep = "_") %>% 
  filter(loc == "BMR" & Year_collect != 2019) %>% 
  mutate(pop = "pop1") %>% 
  select(-spp:-num) %>% 
  select(-cohort) %>% 
  select(pop,ID,SNP1:SNP999)
#BMR - below river #####
bmr_gp2 <- genepop1 %>% 
  mutate(ID2 = ID) %>% 
  separate(ID2,into = c("spp","loc","num"),sep = "_") %>% 
  filter(loc == "BMR" & Year_collect == 2019) %>% 
  mutate(pop = "pop2") %>% 
  select(-spp:-num) %>% 
  select(-cohort) %>% 
  select(pop,ID,SNP1:SNP999)

#BMR - 2015 cohort ####
bmr_2015 <- genepop1 %>% 
  mutate(ID2 = ID) %>% 
  separate(ID2,into = c("spp","loc","num"),sep = "_") %>% 
  filter(loc == "BMR" & Year_collect != 2019 & cohort == "2015") %>% 
  mutate(pop = "pop3") %>% 
  select(-spp:-num) %>% 
  select(-cohort) %>% 
  select(pop,ID,SNP1:SNP999)

#BMR - 2016 cohort #####
bmr_2016 <- genepop1 %>% 
  mutate(ID2 = ID) %>% 
  separate(ID2,into = c("spp","loc","num"),sep = "_") %>% 
  filter(loc == "BMR" & Year_collect != 2019 & cohort == "2016") %>% 
  mutate(pop = "pop4") %>% 
  select(-spp:-num) %>% 
  select(-cohort) %>% 
  select(pop,ID,SNP1:SNP999)

#BMR - 2017 cohort #####
bmr_2017 <- genepop1 %>% 
  mutate(ID2 = ID) %>% 
  separate(ID2,into = c("spp","loc","num"),sep = "_") %>% 
  filter(loc == "BMR" & Year_collect != 2019 & cohort == "2017") %>% 
  mutate(pop = "pop5") %>% 
  select(-spp:-num) %>% 
  select(-cohort) %>% 
  select(pop,ID,SNP1:SNP999)
#binding black mallard cohorts togeter and writing genepop file#####
bmr_gp <- rbind(bmr_gp1,bmr_gp2,bmr_2015,bmr_2016,bmr_2017)
bmr_snps <- bmr_gp %>% select(-pop:-ID)

genepop_create(SNPs = bmr_snps,df = bmr_gp,output_file = "Output/pm.bmr.genepop.txt",
               title = "Genotyped SNPs on Rapture tag panel for P.Marinus larvae - Black Mallard River")

##Cheboygan#####
colnames(locs) <- c("ID","loc","sub.loc")
locs <- locs %>% select(ID,sub.loc)
che <- che %>% select(ID,cohort)
genepop1 <- merge(genepop,locs)
genepop2 <- merge(genepop1,che)
#Pigeon River - 2016 cohort#####
che_pr16 <- genepop2 %>% 
  mutate(ID2 = ID) %>% 
  separate(ID2,into = c("spp","loc","num"),sep = "_") %>% 
  filter(loc == "CHE" & sub.loc == "Pigeon" & cohort == 2016) %>% 
  mutate(pop = "pop1") %>% 
  select(-spp:-num) %>% 
  select(-sub.loc:-cohort) %>% 
  select(pop,ID,everything())
#Pigeon River - 2017 cohort#####
che_pr17 <- genepop2 %>% 
  mutate(ID2 = ID) %>% 
  separate(ID2,into = c("spp","loc","num"),sep = "_") %>% 
  filter(loc == "CHE" & sub.loc == "Pigeon" & cohort == 2017) %>% 
  mutate(pop = "pop2") %>% 
  select(-spp:-num) %>% 
  select(-sub.loc:-cohort) %>% 
  select(pop,ID,everything())
#CHE - cheboygan river - 2016 cohort#####
che_cr16 <- genepop2 %>% 
  mutate(ID2 = ID) %>% 
  separate(ID2,into = c("spp","loc","num"),sep = "_") %>% 
  filter(loc == "CHE" & sub.loc == "Cheboygan" & cohort == 2016) %>% 
  mutate(pop = "pop3") %>% 
  select(-spp:-num) %>% 
  select(-sub.loc:-cohort) %>% 
  select(pop,ID,everything())

#CHE - cheboygan river - 2016 cohort#####
che_cr17 <- genepop2 %>% 
  mutate(ID2 = ID) %>% 
  separate(ID2,into = c("spp","loc","num"),sep = "_") %>% 
  filter(loc == "CHE" & sub.loc == "Cheboygan" & cohort == 2017) %>% 
  mutate(pop = "pop4") %>% 
  select(-spp:-num) %>% 
  select(-sub.loc:-cohort) %>% 
  select(pop,ID,everything())

che_gp <- rbind(che_pr16,che_pr17,che_cr16,che_cr17)
che_snps <- che_gp %>% select(-pop:-ID)
genepop_create(SNPs = che_snps,df = che_gp,output_file = "Output/pm.che.genepop.txt",
               title = "Genotyped SNPs on Rapture tag panel for P.Marinus larvae - Cheboygan River")


#OCQ - 2015 cohort#####
ocq <- ocq %>% select(ID,cohort)
genepop2 <- merge(genepop1,ocq)
ocq_gp15 <- genepop2 %>% 
  mutate(ID2 = ID) %>%
  separate(ID2,into = c("spp","loc","num"),sep = "_") %>% 
  filter(loc == "OCQ" & cohort == 2015) %>% 
  mutate(pop = "pop1") %>% 
  select(-spp:-num) %>% 
  select(pop,ID,SNP1:SNP999)

#OCQ - 2016 cohort#####
ocq_gp16 <- genepop2 %>% 
  mutate(ID2 = ID) %>%
  separate(ID2,into = c("spp","loc","num"),sep = "_") %>% 
  filter(loc == "OCQ" & cohort == 2016) %>% 
  mutate(pop = "pop2") %>% 
  select(-spp:-num) %>% 
  select(pop,ID,SNP1:SNP999)
ocq_gp <- rbind(ocq_gp15,ocq_gp16)
ocq_snps <- ocq_gp %>% select(-pop:-ID)
genepop_create(SNPs = ocq_snps,df = ocq_gp,output_file = "Output/pm.ocq.genepop.txt",
               title = "Genotyped SNPs on Rapture tag panel for P.Marinus larvae - Ocqueoc River")


