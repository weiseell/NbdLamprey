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
load("SNPsets/COLONY_SNPsets_062321.rda")
#!# Note: Script requires cohort summary files from script 3b to run

#Goal 1####
#making a COLONY input file
#changing gt_Col into a three col data frame
#only ID, Individual, and gt
gts <- read.table(paste0("Input_fulldata/BMR_GT_8X.GT.FORMAT"),header = T,sep = "\t",stringsAsFactors = F)
gts$ID <- paste0(gts$CHROM,"-",gts$POS)
gt_Col <- merge(Colony_out[["BMR"]],gts)
col_input <- gt_Col %>% 
  #making generic loci names for COLONY here
  mutate(locus = paste0("L",1:nrow(gt_Col))) %>% 
  select(-ID:-MAF) %>% 
  gather(key = "ID",value = "gt",-locus)
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
                 prob.mom = 0,prob.dad = 0,output_file = "SNPsets/BMR16_colony2_110121.dat")

col_BMR15 <- col_format[which(col_format$id %in% bmr_cohort15$OffspringID),]
colonydat_create(moms = NA,dads = NA,kids = col_BMR15,markers = markers,update.alfs = 0,spp.type = 2,inbreeding = 0,
                 ploidy = 0,fem.gamy = 0,mal.gamy = 0,clone = 0,sib.scale = 0,sib.prior = 0,known.alfs = 0,
                 run.number = 1,run.length = 2,monitor = 1,windows.version = 1,full.likelihood = 1,likelihood.precision = 3,
                 prob.mom = 0,prob.dad = 0,output_file = "SNPsets/BMR15_colony2_110121.dat")

col_BMRAL <- col_format[which(col_format$id %in% bmr19$OffspringID),]
colonydat_create(moms = NA,dads = NA,kids = col_BMRAL,markers = markers,update.alfs = 0,spp.type = 2,inbreeding = 0,
                 ploidy = 0,fem.gamy = 0,mal.gamy = 0,clone = 0,sib.scale = 0,sib.prior = 0,known.alfs = 0,
                 run.number = 1,run.length = 2,monitor = 1,windows.version = 1,full.likelihood = 1,likelihood.precision = 3,
                 prob.mom = 0,prob.dad = 0,output_file = "SNPsets/BMRAL_colony2_110121.dat")

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

library(tidyverse)


#Goal 2 ####
#Create LD SNP sets for each population
# and write NeEstimator input files
#load homebrew functions:
source("Homebrew/match_tags.R")
source("Homebrew/ld_filter.R")
source("Homebrew/vcf_genepop.R")
source("Homebrew/genepop_create.R")

load("Input/rapture_panel_all_SNPs.rda")
#getting unique tags
unique_tags <- unique(rapture$rad_tag_name)
rapture1 <- as.data.frame(unique_tags)

#separating tag IDs to make location ranges on reference genome scaffolds for all tags
rapture1 <- rapture1 %>% 
  separate(unique_tags, into = c("CHROM","pos"),sep = ":") %>% 
  separate(pos, into = c("min","max"),sep = "-") %>% 
  mutate(ID = paste0(CHROM,":",min,"-",max)) %>% 
  select(ID, CHROM, min, max)

pops <- c("BMR","CHE","OCQ")

LD_out <- list()
for (i in 2:length(pops)) {
  print(i)
  pop <- pops[i]
  gts <- read.table(paste0("Input_fulldata/",pop,"_GT_8X.GT.FORMAT"),header = T,sep = "\t",stringsAsFactors = F)
  gts$ID <- paste0(gts$CHROM,"-",gts$POS)
  MAF <- read.table(paste0("Input_fulldata/",pop,"_allele.frq"),header = T,sep = "\t",stringsAsFactors = F)
  MAF$ID <- paste(MAF$CHROM,MAF$POS,sep = "-")
  #calculate het and pGT for SNP selection
  #calculate heterozygosity
  #n_indiv <- ncol(gts)-3
  #stats <- data.frame(ID = gts$ID,stringsAsFactors = F)
  #het_counts <- rowSums(gts == "0/1")
  #stats$het <- het_counts/n_indiv
  #calculate percent coverage
  #gt_missing <- rowSums(gts == "./.")
  #stats$pGT <- (1-(gt_missing/n_indiv))
  #merge with MAF
  #stats <- merge(stats,MAF)
  
  #SNPs <- gts %>% 
  #  separate(ID,into = c("CHROM","POS"),sep = "-") %>% 
  #  select(CHROM,POS)
  
  #target <- match_tags(SNPs = SNPs,tags = rapture1,target = T)
  
  #merging target results with summary statistics
  #rapture_class <- target[[1]]
  #rapture_class <- rapture_class %>% 
  #  mutate(ID = paste(CHROM,POS,sep = "-")) %>% 
  #  select(ID,everything())
  #SNPsumm <- merge(rapture_class,stats)
  
  #prefilter <- SNPsumm %>%
  #  rename(MAF=MinAF) %>% 
  #  filter(target != "NonTarget")
  #LD_out[[pop]] <- LD_filter(prefilter)
  
  #making an NeEstimator file
  gt_LD <- merge(LD_out[[pop]],gts)
  #formatting for input 
  snpframe <- gt_LD %>% 
    select(CHROM) %>% 
    mutate(SNP=paste0("SNP",seq(1,length(gt_LD$CHROM),1)))
  write.table(snpframe,file = paste0("SNPsets/",pop,"_scaf_SNP_key_110221.txt"),append = F,quote = F,sep = "\t")
  
  gp_input <- gt_LD %>% 
    mutate(SNP=paste0("SNP",seq(1,length(gt_LD$CHROM),1))) %>% 
    select(-ID:-MAF) %>% 
    gather(key = "indiv",value = "gt",-SNP)
  #changing formatting from .vcf to genepop
  gp_input <- vcf_genepop(gp_input)
  #CHEPR cohort subselection
  #gp_input <- gp_input[which(gp_input$indiv%in%chePR$OffspringID),]
  gp_format <- gp_input[grep(pop,gp_input$indiv),]
  
  gp_format$pop <- pop
  gp_format <- gp_format %>% 
    gather(key = "SNP",value = "gt",-indiv,-pop)
  
  #making the NeEstimator file
  genepop_create(gp_format,output_file = paste0("SNPsets/",pop,"_Neestimator_110221.txt"),title = paste0(pop,"NeEstimator"))
  
}

save(LD_out,file = "SNPsets/LD_SNPsets_062421.rda")

#get ld for bmr subpopulations
cohorts <- c("bmr_cohort15","bmr_cohort16","bmr19")
gts <- read.table(paste0("Input_fulldata/BMR_GT_8X.GT.FORMAT"),header = T,sep = "\t",stringsAsFactors = F)
gts$ID <- paste0(gts$CHROM,"-",gts$POS)
MAF <- read.table(paste0("Input_fulldata/BMR_allele.frq"),header = T,sep = "\t",stringsAsFactors = F)
MAF$ID <- paste(MAF$CHROM,MAF$POS,sep = "-")

#making an NeEstimator file
gt_LD <- merge(LD_out[["BMR"]],gts)
#formatting for input 
snpframe <- gt_LD %>% 
  select(CHROM) %>% 
  mutate(SNP=paste0("SNP",seq(1,length(gt_LD$CHROM),1)))
write.table(snpframe,file = paste0("SNPsets/BMR_scaf_SNP_key_110121.txt"),append = F,quote = F,sep = "\t")

gp_input <- gt_LD %>% 
  mutate(SNP=paste0("SNP",seq(1,length(gt_LD$CHROM),1))) %>% 
  select(-CHROM:-MAF) %>% 
  gather(key = "indiv",value = "gt",-SNP)
#changing formatting from .vcf to genepop
gp_input <- vcf_genepop(gp_input)
bmr_cohort15$pop <- "BMR15"
bmr_cohort16$pop <- "BMR16"
bmr19$pop <- "BMRal"
gp_format <- rbind(bmr_cohort15,bmr_cohort16,bmr19)
gp_format <- gp_format %>% 
  select(OffspringID,pop) %>% 
  rename(indiv=OffspringID) %>% 
  full_join(gp_input,by = "indiv")
gp_format <- gp_format %>% 
  gather(key = "SNP",value = "gt",-indiv,-pop)

#making the NeEstimator file
genepop_create(gp_format,output_file = paste0("SNPsets/BMR_Neestimator_110221.txt"),title = "Black Mallard NeEstimator")
