#Summary figures for SNP data set
#Created by: Ellie Weise
#Originally Created on: Nov 5th, 2019
#Last Edited on: Nov 5th, 2019

#Goals:
#Create the following figures:

#libraries
library(tidyverse)
library(anchors)

#homebrew functions


#load in data
load(file = "Input/gt_summary_with_targets.rda")
load(file = "Input/selected_loci_ld_filter_summary.rda")
load(file = "Input/tag_selected_SNPs.rda")
load(file = "Input/rapture_panel_all_SNPs.rda")

#getting depth per SNP and combining with all summaries
gdepth1 <- gdepth %>% mutate(ID = paste0(CHROM,"_",POS))

gdepth1 <- gdepth1 %>% 
  select(-CHROM:-POS) %>% 
  gather(key = "indiv",value = "depth",-ID)

gdepth1 <- replace.value(gdepth1,names = "depth",from = -1, to = 0)

depth_summ <- gdepth1 %>% 
  group_by(ID) %>% 
  summarize(depth_mean = mean(depth))

geno1 <- geno1 %>% 
  mutate(ID = paste0(CHROM,"_",POS))
  
geno1 <- geno1 %>% 
  select(-CHROM:-POS) %>% 
  select(ID,everything())

geno1 <- merge(geno1,depth_summ)
loci_select_summ <- loci_select_summ %>% mutate(ID = paste0(CHROM,"_",POS))
loci_select_summ <- merge(loci_select_summ,depth_summ)

SNP_selected <- SNP_selected %>% mutate(ID = paste0(CHROM,"_",POS))
SNP_selected <- merge(SNP_selected,depth_summ)

#calculating proportion of reads that are on target
gdepth2 <- gdepth1 %>% spread(key = indiv,value = depth)
target <- geno1 %>% 
  mutate(ID = paste0(CHROM,"_",POS))
target <- data.frame(ID = target$ID,target = target$target,stringsAsFactors = F)
target_SNPs <- target %>% 
  filter(target != "NonTarget")
gdepth2 <- merge(target_SNPs,gdepth2)
gdepth2 <- gdepth2 %>% 
  select(-target) %>% 
  gather(key = "OffspringID",value = "depth",-ID)
sum(gdepth1$depth)
sum(gdepth2$depth)
#filtering for on-target to compare
geno2 <- geno1 %>% 
  filter(target != "NonTarget")
geno3 <- geno1 %>% 
  filter(target == "NonTarget")
#figuring out how many SNPs per tag 
rapture1 <- unique(rapture$rad_tag_name)

tags <- data.frame(matrix(data = NA ,nrow = length(rapture1),ncol = 2))
colnames(tags) <- c("tag", "count")
i <- 1
for (i in 1:length(rapture1)) {
  tags$tag[i] <- rapture1[i]
  tmp <- which(geno1$target == tags$tag[i])
  tags$count[i] <- length(tmp)
}
#histogram of SNPs per Rapture tag
tiff(file="Output/hist_SNPs_per_tag.tiff",width=6, height=4, units="in", res=200)
hist(tags$count,breaks = 50,
     main = "Distribution of genotyped \nSNPs per Rapture tag",
     xlab = "SNPs per tag")
dev.off()


#boxplots of summary stats for:
  #all SNPs
  #off-target
  #on-target SNPs
  #selected SNPs for COLONY
  #selected SNPs for NeEstimator
tiff(file="Output/boxplot_pGT.tiff",width=6, height=4, units="in", res=200)
boxplot(geno1$pGT,geno2$pGT,geno3$pGT,loci_select_summ$pGT,SNP_selected$pGT,
        main = "Distribution of Percent Genotyped",
        names = c("All","On-\nTarget","Off-\nTarget","SF/ \n PwoP","LD"),
        ylab = "Percent Genotyped")
dev.off()

tiff(file="Output/boxplot_MAF.tiff",width=6, height=4, units="in", res=200)
boxplot(geno1$MAF,geno2$MAF,geno3$MAF,loci_select_summ$MAF,SNP_selected$MAF,
        main = "Distribution of Minor Allele Frequency",
        names = c("All","On-\nTarget","Off-\nTarget","SF/ \n PwoP","LD"),
        ylab = "Minor Allele Frequency")
dev.off()

tiff(file="Output/boxplot_het.tiff",width=6, height=4, units="in", res=200)
boxplot(geno1$het,geno2$het,geno3$het,loci_select_summ$het,SNP_selected$het,
        main = "Distribution of Heterozygosity",
        names = c("All","On-\nTarget","Off-\nTarget","SF/ \n PwoP","LD"),
        ylab = "Heterozygosity")
dev.off()

tiff(file="Output/boxplot_depth.tiff",width=6, height=4, units="in", res=200)
boxplot(geno1$depth_mean,geno2$depth_mean,geno3$depth_mean,loci_select_summ$depth_mean,SNP_selected$depth_mean,
        main = "Distribution of Sequencing Depth",
        names = c("All","On-\nTarget","Off-\nTarget","SF/ \n PwoP","LD"),
        ylab = "Average Depth")
dev.off()

