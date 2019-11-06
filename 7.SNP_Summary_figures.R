#Summary figures for SNP data set
#Created by: Ellie Weise
#Originally Created on: Nov 5th, 2019
#Last Edited on: Nov 5th, 2019

#Goals:
#Create the following figures:

#libraries
library(tidyverse)

#homebrew functions


#load in data
load(file = "Input/gt_summary_with_targets.rda")
load(file = "Input/selected_loci_ld_filter_summary.rda")
load(file = "Input/tag_selected_SNPs.rda")
load(file = "Input/rapture_panel_all_SNPs.rda")
load(file = "Input/gdepth_filtered.rda")

#getting depth per SNP and combining with gt_summary
gdepth <- gdepth %>% mutate(ID = paste0(CHROM,"_",POS))
gdepth %>% 
  select(-CHROM:-POS) %>% 
  gather(key = "indiv",value = "depth",-ID) %>% 
  group_by(ID) %>% 
  summarize(depth_mean = mean(depth))
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
for (i in 1:length(unique_tags)) {
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
        names = c("All","On-Target","Off-Target","COLONY","NeEstimator"),
        ylab = "Percent Genotyped")
dev.off()

tiff(file="Output/boxplot_MAF.tiff",width=6, height=4, units="in", res=200)
boxplot(geno1$MAF,geno2$MAF,geno3$MAF,loci_select_summ$MAF,SNP_selected$MAF,
        main = "Distribution of Minor Allele Frequency",
        names = c("All","On-Target","Off-Target","COLONY","NeEstimator"),
        ylab = "Minor Allele Frequency")
dev.off()

tiff(file="Output/boxplot_het.tiff",width=6, height=4, units="in", res=200)
boxplot(geno1$het,geno2$het,geno3$het,loci_select_summ$het,SNP_selected$het,
        main = "Distribution of Heterozygosity",
        names = c("All","On-Target","Off-Target","COLONY","NeEstimator"),
        ylab = "Heterozygosity")
dev.off()



