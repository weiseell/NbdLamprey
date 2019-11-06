#LD filtering on SNP set
#Created by: Ellie Weise
#Originally Created on: Oct 25th, 2019
#Last Edited on: Nov 4th, 2019
#Goals:

#libraries

#homebrew functions
source("Homebrew/ld_filter.R")

#load in data
load(file = "Input/gt_filtered.rda")
load(file = "Input/gt_summary_targets.rda")

loci_select_all <- ld_filter(summ = geno1, gt = geno)
loci_select <- loci_select_all[[2]]
loci_select_summ <- loci_select_all[[1]]

#output selected loci file
save(loci_select,file = "Input/selected_loci_ld_filter.rda")
