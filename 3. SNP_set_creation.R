#SNP set creation
#Main Objectives:
#1. Make a SNP data set for LD method and pedigree methods of Nb
#2. Make GENEPOP and Colony files for all populations
#load libraries:
library(tidyverse)

#load homebrew functions:
source("Homebrew/LD_filter.R")
source("Homebrew/ped_filter.R")
source("Homebrew/marker_create.R")
source("Homebrew/vcf_colony.R")
source("Homebrew/colonydat_create.R")
source("Homebrew/vcf_genepop.R")
source("Homebrew/genepop_create.R")
#load input data:
load("Summary_Stats/SNP_summaries_targets.rda")
load("Input/test_genotypes.rda")
#Goal 1####
#Selecting loci for pedigree methods
#selecting for independent loci
#also want to maximize MAF and coverage for selected SNPs
#prefilter for Rapture tag SNPs and minimum MAF of 0.05
prefilter <- SNPsumm %>% 
  filter(target != "NonTarget") %>% 
  filter(MAF > 0.05)
Col_select <- ped_filter(prefilter,window = 1000000,pGT_min = 0.8)
#a few summary histograms to look at the SNP set
hist(Col_select$MAF)
hist(Col_select$pGT)

#selcting loci for LD method
#select one SNP per tag
LD_select <- LD_filter(SNPsumm)
#a few summary histograms to look at the SNP set
hist(LD_select$MAF)
hist(LD_select$pGT)

#Goal 2####
#merging SNP sets with the corresponding genotype calls
gt_Col <- merge(Col_select,gt)
gt_LD <- merge(LD_select,gt)

#making a COLONY input file
#changing gt_Col into a three col data frame
#only ID, Individual, and gt
col_input <- gt_Col %>% 
  #making generic loci names for COLONY here
  mutate(locus = paste0("L",1:nrow(gt_Col))) %>% 
  select(-ID:-rvar) %>% 
  gather(key = "id",value = "gt",-locus)
#converting gts from .vcf format to COLONY format
col_format <- vcf_colony(col_input)
#making the markers file for COLONY
SNPs <- colnames(col_format)
SNPs <- SNPs[-1]
markers <- marker_create(SNPs,cod = 0,gte = 0.02,ote = 0.001)
#writing COLONY file
colonydat_create(moms = NA,dads = NA,kids = col_format,markers = markers,update.alfs = 0,spp.type = 2,inbreeding = 0,
                 ploidy = 0,fem.gamy = 0,mal.gamy = 0,clone = 0,sib.scale = 0,sib.prior = 0,known.alfs = 0,
                 run.number = 1,run.length = 2,monitor = 1,windows.version = 1,full.likelihood = 1,likelihood.precision = 3,
                 prob.mom = 0,prob.dad = 0,output_file = "SNPsets/test_colony2.dat")

#making an NeEstimator file
#formatting for input
gp_input <- gt_LD %>% 
  select(-CHROM:-MAF) %>% 
  gather(key = "indiv",value = "gt",-ID) %>% 
  rename(SNP=ID)
#changing formatting from .vcf to genepop
gp_format <- vcf_genepop(gp_input)
gp_format <- gp_format %>% 
  mutate(pop = "pop1") %>% 
  gather(key = "SNP",value = "gt",-indiv,-pop)

#making the NeEstimator file
genepop_create(gp_format,output_file = "SNPsets/test_Neestimator.txt",title = "Test file for NeEstimator")





