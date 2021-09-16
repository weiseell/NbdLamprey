#NbdLamprey - Script 4: SNP Sets
#Main Objectives:
#1. #1. Make a SNP data set for LD method and pedigree methods of Nb
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
load("Input_fulldata/gt_summary_with_targets.rda")
load("Input_fulldata/gt_filtered.rda")
load("SNPsets/COLONY_SNPsets_062321.rda")

#Goal 2####
#making a COLONY input file
#changing gt_Col into a three col data frame
#only ID, Individual, and gt
gts <- read.table(paste0("Input_fulldata/BMR_GT_8X.GT.FORMAT"),header = T,sep = "\t",stringsAsFactors = F)
gts$ID <- paste0(gts$CHROM,"-",gts$POS)
gt_Col <- merge(Colony_out[["BMR"]],gts)
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
##writing COLONY files for each cohort
col_BMR16 <- col_format[which(col_format$id %in% bmr_cohort16$OffspringID),]
colonydat_create(moms = NA,dads = NA,kids = col_BMR16,markers = markers,update.alfs = 0,spp.type = 2,inbreeding = 0,
                 ploidy = 0,fem.gamy = 0,mal.gamy = 0,clone = 0,sib.scale = 0,sib.prior = 0,known.alfs = 0,
                 run.number = 1,run.length = 2,monitor = 1,windows.version = 1,full.likelihood = 1,likelihood.precision = 3,
                 prob.mom = 0,prob.dad = 0,output_file = "SNPsets/BMR16_colony2_062521.dat")

col_BMR15 <- col_format[which(col_format$id %in% bmr_cohort15$OffspringID),]
colonydat_create(moms = NA,dads = NA,kids = col_BMR15,markers = markers,update.alfs = 0,spp.type = 2,inbreeding = 0,
                 ploidy = 0,fem.gamy = 0,mal.gamy = 0,clone = 0,sib.scale = 0,sib.prior = 0,known.alfs = 0,
                 run.number = 1,run.length = 2,monitor = 1,windows.version = 1,full.likelihood = 1,likelihood.precision = 3,
                 prob.mom = 0,prob.dad = 0,output_file = "SNPsets/BMR15_colony2_062521.dat")

col_BMRAL <- col_format[which(col_format$id %in% bmr19$OffspringID),]
colonydat_create(moms = NA,dads = NA,kids = col_BMRAL,markers = markers,update.alfs = 0,spp.type = 2,inbreeding = 0,
                 ploidy = 0,fem.gamy = 0,mal.gamy = 0,clone = 0,sib.scale = 0,sib.prior = 0,known.alfs = 0,
                 run.number = 1,run.length = 2,monitor = 1,windows.version = 1,full.likelihood = 1,likelihood.precision = 3,
                 prob.mom = 0,prob.dad = 0,output_file = "SNPsets/BMRAL_colony2_062521.dat")

#making a COLONY input file
#changing gt_Col into a three col data frame
#only ID, Individual, and gt
gts <- read.table(paste0("Input_fulldata/CHE_GT_8X.GT.FORMAT"),header = T,sep = "\t",stringsAsFactors = F)
gts$ID <- paste0(gts$CHROM,"-",gts$POS)
gt_Col <- merge(Colony_out[["CHE"]],gts)
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

col_PR <- col_format[which(col_format$id %in% chePR$OffspringID),]
colonydat_create(moms = NA,dads = NA,kids = col_PR,markers = markers,update.alfs = 0,spp.type = 2,inbreeding = 0,
                 ploidy = 0,fem.gamy = 0,mal.gamy = 0,clone = 0,sib.scale = 0,sib.prior = 0,known.alfs = 0,
                 run.number = 1,run.length = 2,monitor = 1,windows.version = 1,full.likelihood = 1,likelihood.precision = 3,
                 prob.mom = 0,prob.dad = 0,output_file = "SNPsets/chePR_colony2_062521.dat")

#Goal 1####
#selcting loci for LD method
#select one SNP per tag
geno_targets <- HDPlotSNPs %>% filter(target != "NonTarget")
LD_select <- LD_filter(geno_targets)
#a few summary histograms to look at the SNP set
hist(LD_select$MAF)
hist(LD_select$pGT)

#save LD SNP set
save(LD_select,file = "SNPsets/NeEst_SNPset.rda")

#making an NeEstimator file
gt_LD <- merge(LD_select,geno)
#formatting for input 
snpframe <- gt_LD %>% 
  select(CHROM) %>% 
  mutate(SNP=paste0("SNP",seq(1,length(gt_LD$CHROM),1)))
write.table(snpframe,file = "SNPsets/scaf_SNP_key_060921.txt",append = F,quote = F,sep = "\t")

gp_input <- gt_LD %>% 
  mutate(SNP=paste0("SNP",seq(1,length(gt_LD$CHROM),1))) %>% 
  select(-CHROM:-target) %>% 
  gather(key = "indiv",value = "gt",-SNP)
#changing formatting from .vcf to genepop
gp_input <- vcf_genepop(gp_input)
gp_format <- gp_input[grep("OCQ",gp_input$indiv),]
gp_format <- gp_format %>% 
  mutate(pop = "OCQ") %>% 
  gather(key = "SNP",value = "gt",-indiv,-pop)

#making the NeEstimator file
genepop_create(gp_format,output_file = "SNPsets/Neestimator_060921.txt",title = "Test file for NeEstimator")





