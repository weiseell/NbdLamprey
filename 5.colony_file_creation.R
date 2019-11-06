#Converting .vcf genotypes to COLONY
#Created by: Ellie Weise
#Originally Created on: Oct 25th, 2019
#Last Edited on: Nov 4th, 2019
#Goals:
#1. convert selected SNPs from vcf to COLONY format
#2. subset formatted file into cohorts and create COLONY files
#libraries

#homebrew functions
source("Homebrew/marker_create.R")
source("Homebrew/vcf_colony.R")
source("Homebrew/colonydat_create.R")

#load in data
load(file = "Input/gt_filtered.rda")
load(file = "Input/gt_summary_targets.rda")

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

#subset colony file by cohort
bmr <- col2 %>% 
  mutate(id2 = id) %>% 
  separate(id2,into = c("spp","loc","indiv"),sep = "_") %>% 
  filter(loc == "BMR") %>% 
  select(-spp:-indiv)

ocq <- col2 %>% 
  mutate(id2 = id) %>% 
  separate(id2,into = c("spp","loc","indiv"),sep = "_") %>% 
  filter(loc == "OCQ") %>% 
  select(-spp:-indiv)

che <- col2 %>% 
  mutate(id2 = id) %>% 
  separate(id2,into = c("spp","loc","indiv"),sep = "_") %>% 
  filter(loc == "CHE") %>% 
  select(-spp:-indiv)
#making a marker file
markers <- marker_create(SNPs = colnames(col2), cod = 0, gte = 0.02, ote = 0.001)

#make colony files for all cohorts
colonydat_create(moms = NA,dads = NA,kids = bmr,markers = markers,update.alfs = 0,spp.type = 2,inbreeding = 0,
                    ploidy = 0,fem.gamy = 0,mal.gamy = 0,clone = 0,sib.scale = 0,sib.prior = 0,known.alfs = 0,
                    run.number = 1,run.length = 2,monitor = 1,windows.version = 1,full.likelihood = 1,likelihood.precision = 3,
                    prob.mom = 0,prob.dad = 0,output_file = "bmr_1MB_colony2.dat")





