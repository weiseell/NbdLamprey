#SNP set creation
#Main Objectives:
#2. Make GENEPOP and Colony files for each location
#load libraries:
library(tidyverse)

#load homebrew functions:
source("Homebrew/COLONY_filter.R")
source("Homebrew/marker_create.R")
source("Homebrew/vcf_colony.R")
source("Homebrew/colonydat_create.R")
#load input data:
load("Summary_Stats/SNP_summaries_targets.rda")
load("Input_fulldata/gt_summary_stats.rda")
load("Input_fulldata/GTs.all.rda")

#Goal 1####
#Selecting loci for pedigree methods
#selecting for independent loci
#also want to maximize MAF and coverage for selected SNPs
#prefilter for Rapture tag SNPs and minimum MAF of 0.05
prefilter <- SNPsumm %>% 
  filter(target != "NonTarget")
#prefilter <- prefilter[grep(pattern = "scaf",SNPsumm$CHROM),]
Col_select <- COLONY_filter(SNPsumm,window = 500000,pGT_min = 0.8,MAF_min = 0.05)
#a few summary histograms to look at the SNP set
hist(Col_select$MAF,xlim = c(0,0.5))
hist(Col_select$pGT)

#Goal 2####
#merging SNP sets with the corresponding genotype calls
gt_Col <- merge(Col_select,geno)

#making a COLONY input file
#changing gt_Col into a three col data frame
#only ID, Individual, and gt
col_input <- gt_Col %>% 
  #making generic loci names for COLONY here
  mutate(locus = paste0("L",1:nrow(gt_Col))) %>% 
  select(-CHROM:-MAF) %>% 
  gather(key = "id",value = "gt",-locus)
#converting gts from .vcf format to COLONY format
col_format <- vcf_colony(col_input)
#making the markers file for COLONY
SNPs <- colnames(col_format)
SNPs <- SNPs[-1]
markers <- marker_create(SNPs,cod = 0,gte = 0.02,ote = 0.001)
#writing COLONY files
#Note: file appends
col_BMR <- col_format[grep(pattern = "BMR",x = col_format$id),]
colonydat_create(moms = NA,dads = NA,kids = col_BMR,markers = markers,update.alfs = 0,spp.type = 2,inbreeding = 0,
                 ploidy = 0,fem.gamy = 0,mal.gamy = 0,clone = 0,sib.scale = 0,sib.prior = 0,known.alfs = 0,
                 run.number = 1,run.length = 2,monitor = 1,windows.version = 1,full.likelihood = 1,likelihood.precision = 3,
                 prob.mom = 0,prob.dad = 0,output_file = "SNPsets/BMR_092520_colony2.dat")

col_OCQ <- col_format[grep(pattern = "OCQ",x = col_format$id),]
colonydat_create(moms = NA,dads = NA,kids = col_OCQ,markers = markers,update.alfs = 0,spp.type = 2,inbreeding = 0,
                 ploidy = 0,fem.gamy = 0,mal.gamy = 0,clone = 0,sib.scale = 0,sib.prior = 0,known.alfs = 0,
                 run.number = 1,run.length = 2,monitor = 1,windows.version = 1,full.likelihood = 1,likelihood.precision = 3,
                 prob.mom = 0,prob.dad = 0,output_file = "SNPsets/OCQ_092520_colony2.dat")

col_CHE <- col_format[grep(pattern = "CHE",x = col_format$id),]
colonydat_create(moms = NA,dads = NA,kids = col_CHE,markers = markers,update.alfs = 0,spp.type = 2,inbreeding = 0,
                 ploidy = 0,fem.gamy = 0,mal.gamy = 0,clone = 0,sib.scale = 0,sib.prior = 0,known.alfs = 0,
                 run.number = 1,run.length = 2,monitor = 1,windows.version = 1,full.likelihood = 1,likelihood.precision = 3,
                 prob.mom = 0,prob.dad = 0,output_file = "SNPsets/CHE_092520_colony2.dat")

#save SNP set for use in COLONY cohort runs (script 4)
save(Col_select,file = "SNPsets/COLONY_SNPset.rda")

