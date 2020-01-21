#Converting .vcf genotypes to COLONY
#Created by: Ellie Weise
#Originally Created on: Oct 25th, 2019
#Last Edited on: Nov 4th, 2019
#Goals:
#1. convert selected SNPs from vcf to COLONY format
#2. subset formatted file into cohorts and create COLONY files
#libraries
library(tidyverse)
library(adegenet)
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
  mutate(locus = paste(CHROM,POS,sep = "_")) %>% 
  select(-CHROM:-target) %>%
  gather(key = "id",value = "gt",-locus)

#converting vcf to pca format
col2 <- vcf_colony(vcf = col1)
loci <- data.frame(id = colnames(col2),stringsAsFactors = F)
loci <- loci %>% 
  separate(id, into = c("scaf","scafnum","pos"),sep = "_") %>% 
  mutate(scaffold = paste0(scaf,scafnum)) %>% 
  select(-scaf:-scafnum)
locs <- col2 %>% select(id) %>% separate(id, into = c("spp","loc","indiv"),sep = "_") %>% select(loc)
rownames(col2) <- col2$id
col2 <- col2 %>% select(-id)

snp <- new("genlight",
           col2,
           chromosome=loci$scaffold,
           position=loci$pos,
           pop=locs$loc)

#running PCA analysis
#getting colors
locs$col <- NA
for(i in 1:length(locs$loc)){
  if(locs$loc[i] == "BMR"){locs$col[i] <- "blue"}
  if(locs$loc[i] == "OCQ"){locs$col[i] <- "dark green"}
  if(locs$loc[i] == "CHE"){locs$col[i] <- "purple"}
}
pca <- glPca(snp,nf = 10)
par(mfrow=c(1,2))
plot(pca$scores[,1], pca$scores[,2],
     col=locs$col,cex=0.5)
#text(pca$scores[,1], pca$scores[,2] + 0.7,
     #labels=rownames(pca$scores), cex= 0.7)
plot(pca$scores[,1], pca$scores[,3],
     col=locs$col,cex=0.5)
legend("bottomleft",legend = c("BMR","CHE","OCQ"),col = c("blue","dark green","purple"),cex = 0.8)
#text(pca$scores[,1], pca$scores[,3] + 0.7,
     #labels=rownames(pca$scores), cex= 0.7)
