#Nb and Ns Calculations
#Main objectives:
#1.Calculate Nb using the Parentage-without-Parents method
#2.Calculate Ns and extrapolated Ns

#libraries
library(tidyverse)
library(ggrepel)
#homebrew functions
source("Homebrew/PwoP.R")
source("Homebrew/PwoP_uncert.R")
source("Homebrew/Ns_calc.R")
source("Homebrew/multiplot.R")
#load in data
load("Aging_Models/Family_data_all_locations.rda")

#prepping lists for storing results
locs <- unique(all_families$loc)
ca_names <- c("bmr15","bmr16","bmral","chePR","ocq")
names(loc_names) <- locs
Nb_PwoP_all <- data.frame(matrix(data = NA, nrow = length(locs), ncol = 6))
colnames(Nb_PwoP_all) <- c("loc","Nb","kbar","Vk", "CI_Lower", "CI_Upper")
Ns_all <- data.frame(matrix(data = NA, nrow = length(locs), ncol = 6))
colnames(Ns_all) <- c("loc","Ns" ,"Ns_Chao","Chao_uncert","Ns_Jackknife","Jackknife_uncert")

#loop to calculate Nb - PwoP and extrapolated Ns
i <- 1
for (i in 1:length(locs)) {
  #PwoP
  ltmp <- locs[i]
  tmp <- subset(all_families,all_families$loc == ltmp)
  PwoP_tmp <- PwoP(tmp)
  ca_tmp <- readLines(paste0("Software_outputs/",ca_names[i],".ConfigArchive"))
  uncert <- PwoP_uncert(ca,tmp)
  names(uncert) <- c("CI_Lower","CI_Upper")
  tmp1 <- c(PwoP_tmp,uncert)
  Nb_PwoP_all[i,] <- c(ltmp,tmp1)
  
  #Ns
  tmp2 <- tmp[!duplicated(tmp$OffspringID),]
  tmp2 <- na.omit(tmp2)
  Ns_tmp <- Ns_calc(tmp2)
  
  #saving values into a table
  Ns_all[i,] <- c(ltmp,Ns_tmp[[3]],Ns_tmp[[2]]$chao,Ns_tmp[[2]]$chao.se,Ns_tmp[[2]]$jack1,Ns_tmp[[2]]$jack1.se)
  
}


#!# Need to make an Nb/Ns table and save plot (also need to do this for script 4)