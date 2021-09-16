#NbdLamprey - Script 1
##Main Objective: Calculate the following summary statistics:
#1. Calculate percent coverage and heterozygosity for all SNPs
#2. Match SNPs to the Rapture Panel and calculate stats for Rapture success

#load libraries:
library(tidyverse)

#load homebrew functions:
source("Homebrew/match_tags.R")

#load input data:
geno <- read.table("Input_062421/GT_8X.GT.FORMAT",header = T,sep = "\t")
save(geno,file = "Input_fulldata/GTs.all.rda")
af <- read.table("Input_062421/allele.frq",header = T,sep = "\t")
save(af,file = "Input_fulldata/MAF.all.rda")
gdepth <- read.table("Input_fulldata/depth.gdepth",header = T,sep = "\t")
save(gdepth,file = "Input_fulldata/depth.all.rda")
load("Input/rapture_panel_all_SNPs.rda")

#Goal 1####
#calculating %coverage, heterozygosity for all loci
#merging with MAF data
#should have a large table with all genotypes and three summary criteria

#making ID column for genotypes and gdepth
geno <- geno %>% 
  mutate(ID=paste(CHROM,POS,sep="_")) %>% 
  select(ID,everything()) %>% 
  select(-CHROM,-POS)

gdepth <- gdepth %>% 
  mutate(ID=paste(CHROM,POS,sep="_")) %>% 
  select(ID,everything()) %>% 
  select(-CHROM,-POS)

af <- af %>% 
  mutate(ID=paste(CHROM,POS,sep="_")) %>% 
  select(ID,everything()) %>% 
  select(-CHROM,-POS)

#replace -1 values with 0 in gdepth file
gdepth[gdepth==-1]<-0
#mean depth per SNP
rownames(gdepth) <- gdepth$ID
mean_depth <- rowMeans(gdepth[,-c(1)])
hist(mean_depth)
#!#resave gdepth .Rdata object because this can take a long time
save(gdepth,file = "Input_fulldata/depth.all.rda")
##calculate mean depth for each individual

#Checking MAF file: if there are MAF > 0.5, fix so they're < 0.5
af$Min_AF2 <- af$Min_AF
for(i in 1:length(af$Maj_AF)){
  if (af$Min_AF[i] > 0.5) {
    af$Min_AF2[i] <- 1-af$Min_AF[i]
  }
}
#check that all MAF values are below 0.5 and saving as an individual file
max(af$MinAF)
MAF <- data.frame(ID = af$ID, MAF = af$MinAF, stringsAsFactors = F)

#calculate heterozygosity
n_indiv <- ncol(geno)-2
stats <- data.frame(ID = geno$ID,stringsAsFactors = F)
het_counts <- rowSums(geno == "0/1")
stats$het <- het_counts/n_indiv
#calculate percent coverage
gt_missing <- rowSums(geno == "./.")
stats$pGT <- (1-(gt_missing/n_indiv))
#merge with MAF
stats <- merge(stats,MAF)


#Goal 2####
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
SNPs <- geno %>% 
  separate(ID,into = c("scaf","c","POS"),sep = "_") %>% 
  mutate(CHROM = paste0(scaf,"_",c)) %>% 
  select(CHROM,POS)

target <- match_tags(SNPs = SNPs,tags = rapture1,target = T)

#merging target results with summary statistics
rapture_class <- target[[1]]
rapture_class <- rapture_class %>% 
  mutate(ID = paste(CHROM,POS,sep = "_")) %>% 
  select(ID,everything())
SNPsumm <- merge(rapture_class,stats)
#save two main outputs from this script
#summary stats, target classifications, and on-target SNPs

save(SNPsumm,file = "Summary_Stats/SNP_summaries_targets.rda")










