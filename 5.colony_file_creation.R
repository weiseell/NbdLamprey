#Converting .vcf genotypes to COLONY
#Created by: Ellie Weise
#Originally Created on: Oct 25th, 2019
#Last Edited on: Nov 4th, 2019
#Goals:
#1. convert selected SNPs from vcf to COLONY format
#2. subset formatted file into cohorts and create COLONY files
#libraries
library(tidyverse)
#homebrew functions
source("Homebrew/marker_create.R")
source("Homebrew/vcf_colony.R")
source("Homebrew/colonydat_create.R")

#load in data
load(file = "Input/gt_filtered.rda")
load(file = "Input/selected_loci_ld_filter.rda")
locs <- read.table("Input/locs_and_sublocations.txt",header = T,sep = "\t",stringsAsFactors = F)
#convert data from vcf to COLONY format
#removing specialized names for generic loci names
#gather file so theres only three columns: locus, individual, and genotype
#save this file so you can match loci back later if you want
col1 <- loci_select %>%
  mutate(locus = paste0("L",1:nrow(loci_select))) %>%
  select(-CHROM:-target) %>%
  gather(key = "id",value = "gt",-locus)

#converting vcf to colony
col2 <- vcf_colony(vcf = col1)
col2 <- merge(col2,locs)
bmr <- bmr %>% 
  select(-Sample_number) %>% 
  mutate(id = paste0(species,"_",loc,"_",num)) %>% 
  select(-species:-num) %>% 
  select(-Year_collect:-V2) %>% 
  select(id,everything())
col3 <- merge(col2,bmr)
col3 <- col3 %>% select(-loc)
#making a marker file
markers <- marker_create(SNPs = colnames(col2), cod = 0, gte = 0.02, ote = 0.001)

#subsetting and make colony files for all cohorts
#subset colony file by cohort
#BMR below lake#####
bmr_1718 <- col3 %>% 
  mutate(id2 = id) %>% 
  separate(id2,into = c("spp","loc","indiv"),sep = "_") %>% 
  filter(loc == "BMR" & sub.loc == "BelowLake") %>% 
  select(-spp:-indiv) %>% 
  select(-sub.loc)
colonydat_create(moms = NA,dads = NA,kids = bmr_1718,markers = markers,update.alfs = 0,spp.type = 2,inbreeding = 0,
                 ploidy = 0,fem.gamy = 0,mal.gamy = 0,clone = 0,sib.scale = 0,sib.prior = 0,known.alfs = 0,
                 run.number = 1,run.length = 2,monitor = 1,windows.version = 1,full.likelihood = 1,likelihood.precision = 3,
                 prob.mom = 0,prob.dad = 0,output_file = "bmr1718_colony2.dat")
#BMR above lake - 2019#####
bmr_19 <- col3 %>% 
  mutate(id2 = id) %>% 
  separate(id2,into = c("spp","loc","indiv"),sep = "_") %>% 
  filter(loc == "BMR" & sub.loc == "AboveLake") %>% 
  select(-spp:-indiv) %>% 
  select(-sub.loc)
markers <- marker_create(SNPs = colnames(bmr_19), cod = 0, gte = 0.02, ote = 0.001)
colonydat_create(moms = NA,dads = NA,kids = bmr_19,markers = markers,update.alfs = 0,spp.type = 2,inbreeding = 0,
                 ploidy = 0,fem.gamy = 0,mal.gamy = 0,clone = 0,sib.scale = 0,sib.prior = 0,known.alfs = 0,
                 run.number = 1,run.length = 2,monitor = 1,windows.version = 1,full.likelihood = 1,likelihood.precision = 3,
                 prob.mom = 0,prob.dad = 0,output_file = "bmr19_colony2.dat")
#BMR - 2015 cohort #####
bmr_2015 <- col3 %>% 
  mutate(id2 = id) %>% 
  separate(id2,into = c("spp","loc","indiv"),sep = "_") %>% 
  filter(loc == "BMR" & sub.loc == "BelowLake" & cohort == "2015") %>% 
  select(-spp:-indiv) %>% 
  select(-sub.loc) %>% 
  select(-cohort) %>% 
  select(-Year_collect) %>% 
  select(-V1:-V2)
markers <- marker_create(SNPs = colnames(bmr_2015), cod = 0, gte = 0.02, ote = 0.001)
colonydat_create(moms = NA,dads = NA,kids = bmr_2015,markers = markers,update.alfs = 0,spp.type = 2,inbreeding = 0,
                 ploidy = 0,fem.gamy = 0,mal.gamy = 0,clone = 0,sib.scale = 0,sib.prior = 0,known.alfs = 0,
                 run.number = 1,run.length = 2,monitor = 1,windows.version = 1,full.likelihood = 1,likelihood.precision = 3,
                 prob.mom = 0,prob.dad = 0,output_file = "Output/bmr2015_colony2.dat")
#BMR - 2016 cohort #####
bmr_2016 <- col2 %>% 
  mutate(id2 = id) %>% 
  separate(id2,into = c("spp","loc","indiv"),sep = "_") %>% 
  filter(loc == "BMR" & sub.loc == "BelowLake" & cohort == "2016") %>% 
  select(-spp:-indiv) %>% 
  select(-sub.loc) %>% 
  select(-cohort) %>% 
  select(-Year_collect) %>% 
  select(-V1:-V2)
markers <- marker_create(SNPs = colnames(bmr_2016), cod = 0, gte = 0.02, ote = 0.001)
colonydat_create(moms = NA,dads = NA,kids = bmr_2016,markers = markers,update.alfs = 0,spp.type = 2,inbreeding = 0,
                 ploidy = 0,fem.gamy = 0,mal.gamy = 0,clone = 0,sib.scale = 0,sib.prior = 0,known.alfs = 0,
                 run.number = 1,run.length = 2,monitor = 1,windows.version = 1,full.likelihood = 1,likelihood.precision = 3,
                 prob.mom = 0,prob.dad = 0,output_file = "Output/bmr2016_colony2.dat")
#BMR - 2017 cohort #####
bmr_2017 <- col3 %>% 
  mutate(id2 = id) %>% 
  separate(id2,into = c("spp","loc","indiv"),sep = "_") %>% 
  filter(loc == "BMR" & sub.loc == "BelowLake" & cohort == "2017") %>% 
  select(-spp:-indiv) %>% 
  select(-sub.loc) %>% 
  select(-cohort) %>% 
  select(-Year_collect) %>% 
  select(-V1:-V2)
markers <- marker_create(SNPs = colnames(bmr_2017), cod = 0, gte = 0.02, ote = 0.001)
colonydat_create(moms = NA,dads = NA,kids = bmr_2017,markers = markers,update.alfs = 0,spp.type = 2,inbreeding = 0,
                 ploidy = 0,fem.gamy = 0,mal.gamy = 0,clone = 0,sib.scale = 0,sib.prior = 0,known.alfs = 0,
                 run.number = 1,run.length = 2,monitor = 1,windows.version = 1,full.likelihood = 1,likelihood.precision = 3,
                 prob.mom = 0,prob.dad = 0,output_file = "Output/bmr2017_colony2.dat")
#Ocqueoc - 2015 cohort#####
ocq <- ocq %>% select(ID,cohort)
col2 <- col2 %>% rename(ID = id)
col3 <- merge(col2,ocq)
ocq2015 <- col3 %>% 
  mutate(id2 = ID) %>% 
  separate(id2,into = c("spp","loc","indiv"),sep = "_") %>% 
  filter(loc == "OCQ" & cohort == 2015) %>% 
  select(-spp:-indiv) %>% 
  select(-sub.loc:-cohort)
colonydat_create(moms = NA,dads = NA,kids = ocq2015,markers = markers,update.alfs = 0,spp.type = 2,inbreeding = 0,
                 ploidy = 0,fem.gamy = 0,mal.gamy = 0,clone = 0,sib.scale = 0,sib.prior = 0,known.alfs = 0,
                 run.number = 1,run.length = 2,monitor = 1,windows.version = 1,full.likelihood = 1,likelihood.precision = 3,
                 prob.mom = 0,prob.dad = 0,output_file = "Output/ocq2015_colony2.dat")
#Ocqueoc - 2016 cohort#####
ocq <- ocq %>% select(ID,cohort)
col2 <- col2 %>% rename(ID = id)
col3 <- merge(col2,che)
ocq2016 <- col3 %>% 
  mutate(id2 = ID) %>% 
  separate(id2,into = c("spp","loc","indiv"),sep = "_") %>% 
  filter(loc == "OCQ" & cohort == 2016) %>% 
  select(-spp:-indiv) %>% 
  select(-sub.loc:-cohort)
colonydat_create(moms = NA,dads = NA,kids = ocq2016,markers = markers,update.alfs = 0,spp.type = 2,inbreeding = 0,
                 ploidy = 0,fem.gamy = 0,mal.gamy = 0,clone = 0,sib.scale = 0,sib.prior = 0,known.alfs = 0,
                 run.number = 1,run.length = 2,monitor = 1,windows.version = 1,full.likelihood = 1,likelihood.precision = 3,
                 prob.mom = 0,prob.dad = 0,output_file = "Output/ocq2016_colony2.dat")
#Cheboygan ####
che <- che %>% select(ID,cohort)
col2 <- col2 %>% rename(ID = id)
col3 <- merge(col2,che)
#Pigeon River - 2016 #####
chePR_2016 <- col3 %>% 
  mutate(id2 = ID) %>% 
  separate(id2,into = c("spp","loc","indiv"),sep = "_") %>% 
  filter(loc == "CHE" & sub.loc == "Pigeon" & cohort == 2016) %>% 
  select(-spp:-indiv) %>% 
  select(-sub.loc:-cohort)
markers <- marker_create(SNPs = colnames(chePR_2016), cod = 0, gte = 0.02, ote = 0.001)

colonydat_create(moms = NA,dads = NA,kids = chePR_2016,markers = markers,update.alfs = 0,spp.type = 2,inbreeding = 0,
                 ploidy = 0,fem.gamy = 0,mal.gamy = 0,clone = 0,sib.scale = 0,sib.prior = 0,known.alfs = 0,
                 run.number = 1,run.length = 2,monitor = 1,windows.version = 1,full.likelihood = 1,likelihood.precision = 3,
                 prob.mom = 0,prob.dad = 0,output_file = "Output/chePR_2016_colony2.dat")
#Pigeon River - 2017####
chePR_2017 <- col3 %>% 
  mutate(id2 = ID) %>% 
  separate(id2,into = c("spp","loc","indiv"),sep = "_") %>% 
  filter(loc == "CHE" & sub.loc == "Pigeon" & cohort == 2017) %>% 
  select(-spp:-indiv) %>% 
  select(-sub.loc:-cohort)
markers <- marker_create(SNPs = colnames(chePR_2017), cod = 0, gte = 0.02, ote = 0.001)

colonydat_create(moms = NA,dads = NA,kids = chePR_2017,markers = markers,update.alfs = 0,spp.type = 2,inbreeding = 0,
                 ploidy = 0,fem.gamy = 0,mal.gamy = 0,clone = 0,sib.scale = 0,sib.prior = 0,known.alfs = 0,
                 run.number = 1,run.length = 2,monitor = 1,windows.version = 1,full.likelihood = 1,likelihood.precision = 3,
                 prob.mom = 0,prob.dad = 0,output_file = "Output/chePR_2017_colony2.dat")
#Cheboygan River - 2016 cohort #####
cheCR_2016 <- col3 %>% 
  mutate(id2 = ID) %>% 
  separate(id2,into = c("spp","loc","indiv"),sep = "_") %>% 
  filter(loc == "CHE" & sub.loc == "Cheboygan" & cohort == 2016) %>% 
  select(-spp:-indiv) %>% 
  select(-sub.loc:-cohort)
markers <- marker_create(SNPs = colnames(cheCR_2016), cod = 0, gte = 0.02, ote = 0.001)

colonydat_create(moms = NA,dads = NA,kids = cheCR_2016,markers = markers,update.alfs = 0,spp.type = 2,inbreeding = 0,
                 ploidy = 0,fem.gamy = 0,mal.gamy = 0,clone = 0,sib.scale = 0,sib.prior = 0,known.alfs = 0,
                 run.number = 1,run.length = 2,monitor = 1,windows.version = 1,full.likelihood = 1,likelihood.precision = 3,
                 prob.mom = 0,prob.dad = 0,output_file = "Output/cheCR_2016_colony2.dat")

#Cheboygan River - 2017 cohort #####
cheCR_2017 <- col3 %>% 
  mutate(id2 = ID) %>% 
  separate(id2,into = c("spp","loc","indiv"),sep = "_") %>% 
  filter(loc == "CHE" & sub.loc == "Cheboygan" & cohort == 2017) %>% 
  select(-spp:-indiv) %>% 
  select(-sub.loc:-cohort)
markers <- marker_create(SNPs = colnames(cheCR_2017), cod = 0, gte = 0.02, ote = 0.001)

colonydat_create(moms = NA,dads = NA,kids = cheCR_2017,markers = markers,update.alfs = 0,spp.type = 2,inbreeding = 0,
                 ploidy = 0,fem.gamy = 0,mal.gamy = 0,clone = 0,sib.scale = 0,sib.prior = 0,known.alfs = 0,
                 run.number = 1,run.length = 2,monitor = 1,windows.version = 1,full.likelihood = 1,likelihood.precision = 3,
                 prob.mom = 0,prob.dad = 0,output_file = "Output/cheCR_2017_colony2.dat")



