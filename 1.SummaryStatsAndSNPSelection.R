##NbdLamprey - Script 1
#Script goals: 
#1. Calculating the following summary statistics for each population:
#minor allele frequency, heterozygosity, percent genotyped
#2. Select SNPs that have a MAF of 0.05 and >80% genotyped 
#for per population COLONY run

#Written by: Ellie Weise
#Last edited: 10/26/21

#required libraries
library(tidyverse)

#load homebrew functions:
source("Homebrew/match_tags.R")
source("Homebrew/COLONY_filter.R")
source("Homebrew/marker_create.R")
source("Homebrew/vcf_colony.R")
source("Homebrew/colonydat_create.R")

#load RAD-capture reference file
load("Input/rapture_panel_all_SNPs.rda")
#generate pops list and input files
pops <- c("BMR","CHE","OCQ")
Colony_out <- list()
for (i in 3:length(pops)) {
  pop <- pops[i]
  print(pop)
  ##Goal 1####
  print("Calculating Summary Stats")
  #read in genotype table from VCFtools script
  gts <- read.table(paste0("Input_fulldata/",pop,"_GT_8X.GT.FORMAT"),header = T,sep = "\t",stringsAsFactors = F)
  gts$ID <- paste0(gts$CHROM,"-",gts$POS)
  #load in MAF data from VCFtools
  MAF <- read.table(paste0("Input_fulldata/",pop,"_allele.frq"),header = T,sep = "\t",stringsAsFactors = F)
  MAF$ID <- paste(MAF$CHROM,MAF$POS,sep = "-")
  ##calculate het and pGT for SNP selection
  #calculate heterozygosity
  n_indiv <- ncol(gts)-3
  stats <- data.frame(ID = gts$ID,stringsAsFactors = F)
  het_counts <- rowSums(gts == "0/1")
  stats$het <- het_counts/n_indiv
  #calculate percent coverage
  gt_missing <- rowSums(gts == "./.")
  stats$pGT <- (1-(gt_missing/n_indiv))
  ##merge with MAF
  stats <- stats %>% full_join(MAF,by = "ID")
  
  ##match to the RAD-capture panel for reference
  print("Matching to RAD-capture Panel")
  #getting unique tags
  unique_tags <- unique(rapture$rad_tag_name)
  rapture1 <- as.data.frame(unique_tags)
  
  #separating tag IDs to make location ranges on reference genome scaffolds for all tags
  rapture1 <- rapture1 %>% 
    separate(unique_tags, into = c("CHROM","pos"),sep = ":") %>% 
    separate(pos, into = c("min","max"),sep = "-") %>% 
    mutate(ID = paste0(CHROM,":",min,"-",max)) %>% 
    select(ID, CHROM, min, max)
  SNPs <- gts %>% 
    select(CHROM,POS)
  
  target <- match_tags(SNPs = SNPs,tags = rapture1,target = T)
  
  #merging target results with summary statistics
  rapture_class <- target[[1]]
  SNPsumm <- stats %>% full_join(target[[1]],by = c("CHROM","POS"))
  ##save the summary file
  write.table(SNPsumm,file = paste0("Summary_Stats/",pop,"_Summary_bySNP.txt"),append = F,quote = F,sep = "\t")
  ##Goal 2####
  print("Generating SNP set")
  ##filter SNPs for Colony
  #!#Note: for very small populations these criteria may need to change
  prefilter <- stats %>%
    rename(MAF=MinAF) %>% 
    filter(MAF > 0.05 & pGT > 0.8)
  Colony_out[[pop]] <- COLONY_filter(prefilter,window = 1000000,pGT_min = 0.8,MAF_min = 0.05)
  
  #merging SNP sets with the corresponding genotype calls
  gt_Col <- merge(Colony_out[[pop]],gts)
  
  #making a COLONY input file
  #changing gt_Col into a three col data frame
  #only ID, Individual, and gt
  col_input <- gt_Col %>% 
    #making generic loci names for COLONY here
    mutate(locus = paste0("L",1:nrow(gt_Col))) %>% 
    select(-CHROM:-pGT) %>% 
    gather(key = "id",value = "gt",-locus)
  #converting gts from .vcf format to COLONY format
  col_format <- vcf_colony(col_input)
  #making the markers file for COLONY
  SNPs <- colnames(col_format)
  SNPs <- SNPs[-1]
  markers <- marker_create(SNPs,cod = 0,gte = 0.02,ote = 0.001)
  
  tmp_col <- col_format[grep(pattern = pops[i],x = col_format$id),]
  #colonydat_create(moms = NA,dads = NA,kids = tmp_col,markers = markers,update.alfs = 0,spp.type = 2,inbreeding = 0,
  #                 ploidy = 0,fem.gamy = 0,mal.gamy = 0,clone = 0,sib.scale = 0,sib.prior = 0,known.alfs = 0,
  #                 run.number = 1,run.length = 2,monitor = 1,windows.version = 1,full.likelihood = 1,likelihood.precision = 3,
  #                 prob.mom = 0,prob.dad = 0,output_file = paste0("SNPsets/ColonyFiles_062421/",pops[i],"_colony2.dat"))
}

#save COLONY SNP sets for later analysis and reference
save(Colony_out,file = "SNPsets/COLONY_SNPsets_062321.rda")





