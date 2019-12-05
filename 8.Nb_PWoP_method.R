#Calculating Nb using the Parentage without Parents (PwoP) method (Waples and Waples 2011)
#Created by: Ellie Weise
#Originally Created on: Nov 12th, 2019
#Last Edited on: Nov 12th, 2019

#Goals:
#1.calculate Nb using the PwoP method
#2. Calculate estimate uncertainty with a bootstrapping method

#load libraries
library(tidyverse)

#homebrew functions
source("Homebrew/PwoP.R")
source("Homebrew/PwoP_boot.R")

#load in data
bmr_colony <- read.table("Input/colony.bestconfig.bmr.txt",header = T,sep = "\t",stringsAsFactors = F)
che_colony <- read.table("Input/colony.bestconfig.che.txt",header = T,sep = "\t",stringsAsFactors = F)
ocq_colony <- read.table("Input/colony.bestconfig.ocq.txt",header = T,sep = "\t",stringsAsFactors = F)
locs <- read.table(file = "Input/locs_and_sublocations.txt",sep = "\t",header = T,stringsAsFactors = F)
ocq <- read.table(file = "Input/ocq_cohort_info.txt",sep = "\t",header = T,stringsAsFactors = F)
che <- read.table(file = "Input/che_cohort_info.txt",sep = "\t",header = T,stringsAsFactors = F)
bmr <- read.table(file = "Input/bmr_cohort_info.txt",sep = "\t",header = T,stringsAsFactors = F)
#run extrapolated Ns curves (function) for all cohorts
#BMR - Making files to subset cohorts ####
bmr <- bmr %>%  select(ID,cohort) %>% rename(OffspringID = ID)
locs <- locs %>% rename(OffspringID = id)
bmr <- merge(bmr,locs)
colnames(bmr) <- c("OffspringID","cohort","loc","sub.loc")
bmr_colony <- merge(bmr_colony,bmr)
#BMR - above lake #####
bmr_col_al <- bmr_colony %>% 
  filter(sub.loc == "AboveLake")
PwoP_bmr_al <- PwoP(bmr_col_al)
bmr_al_conf <- PwoP_boot(family = bmr_col_al,iter = 1000,alpha = 0.05,real_Nb = PwoP_bmr_al$Nb)
#BMR - below lake ####
bmr_col_bl <- bmr_colony %>% 
  filter(sub.loc == "BelowLake")
PwoP_bmr_bl <- PwoP(bmr_col_bl)
bmr_bl_conf <- PwoP_boot(family = bmr_col_bl,iter = 1000,alpha = 0.05,real_Nb = PwoP_bmr_bl$Nb)
#BMR - 2015 cohort ####
bmr_col_2015 <- bmr_colony %>% 
  filter(sub.loc == "BelowLake" & cohort == "2015")
PwoP_bmr_2015 <- PwoP(bmr_col_2015)
bmr_2015_conf <- PwoP_boot(family = bmr_col_2015,iter = 1000,alpha = 0.05,real_Nb = PwoP_bmr_2015$Nb)
#BMR - 2016 cohort ####
bmr_col_2016 <- bmr_colony %>% 
  filter(sub.loc == "BelowLake" & cohort == "2016")
PwoP_bmr_2016 <- PwoP(bmr_col_2016)
bmr_2016_conf <- PwoP_boot(family = bmr_col_2016,iter = 1000,alpha = 0.05,real_Nb = PwoP_bmr_2016$Nb)
#BMR - 2017 cohort ####
bmr_col_2017 <- bmr_colony %>% 
  filter(sub.loc == "BelowLake" & cohort == "2017")
PwoP_bmr_2017 <- PwoP(bmr_col_2017)
bmr_2017_conf <- PwoP_boot(family = bmr_col_2017,iter = 1000,alpha = 0.05,real_Nb = PwoP_bmr_2017$Nb)

#CHE - Making files to subset cohorts ####
colnames(locs) <- c("OffspringID","loc","sub.loc")
che_colony <- merge(che_colony,locs)
che_colony <- merge(che_colony,che)
#CHE - Pigeon Rd - 2016 cohort####
chePR16 <- che_colony %>% 
  filter(sub.loc == "Pigeon" & cohort == 2016)
PwoP_chePR16 <- PwoP(chePR16)
chePR16_conf <- PwoP_boot(family = chePR16,iter = 1000,alpha = 0.05,real_Nb = PwoP_chePR16$Nb)

#CHE - Pigeon Rd - 2017 cohort####
chePR17 <- che_colony %>% 
  filter(sub.loc == "Pigeon" & cohort == 2017)
PwoP_chePR17 <- PwoP(chePR17)
chePR17_conf <- PwoP_boot(family = chePR17,iter = 1000,alpha = 0.05,real_Nb = PwoP_chePR17$Nb)

#CHE - Cheboygan River - 2016 cohort####
cheCR16 <- che_colony %>%  
  filter(sub.loc == "Cheboygan" & cohort == 2016)
PwoP_cheCR16 <- PwoP(cheCR16)
cheCR16_conf <- PwoP_boot(family = cheCR16,iter = 1000,alpha = 0.05,real_Nb = PwoP_cheCR16$Nb)

#CHE - Cheboygan River - 2017 cohort####
cheCR17 <- che_colony %>%  
  filter(sub.loc == "Cheboygan" & cohort == 2017)
PwoP_cheCR17 <- PwoP(cheCR17)
cheCR17_conf <- PwoP_boot(family = cheCR17,iter = 1000,alpha = 0.05,real_Nb = PwoP_cheCR17$Nb)
#OCQ ####
ocq_colony <- merge(ocq_colony,locs)
ocq_colony <- merge(ocq_colony,ocq)
PwoP_ocq <- PwoP(ocq_colony)
ocq_conf <- PwoP_boot(family = ocq_colony,iter = 1000,alpha = 0.05,real_Nb = PwoP_ocq$Nb)

#OCQ - 2015 cohort ####
ocq15 <- ocq_colony %>% 
  filter(loc == "OCQ" & cohort == 2015)
PwoP_ocq15 <- PwoP(ocq15)
ocq15_conf <- PwoP_boot(family = ocq15,iter = 1000,alpha = 0.05,real_Nb = PwoP_ocq15$Nb)

#OCQ - 2016 cohort ####
ocq16 <- ocq_colony %>% 
  filter(loc == "OCQ" & cohort == 2016)
PwoP_ocq16 <- PwoP(ocq16)
ocq16_conf <- PwoP_boot(family = ocq16,iter = 1000,alpha = 0.05,real_Nb = PwoP_ocq16$Nb)

#making PwoP summary table for all cohorts ####
PwoP_all <- data.frame(matrix(nrow = 8,ncol = 3))
rownames(PwoP_all) <- c("Black Mallard - above lake",
                        "Black Mallard - below lake",
                        "Black Mallard - 2015 cohort",
                        "Black Mallard - 2016 cohort",
                        "Black Mallard - 2017 cohort",
                        "Cheboygan - Pigeon Rd",
                        "Cheboygan - Cheboygan River",
                        "Ocqueoc")
colnames(PwoP_all) <- c("Nb","kbar","Vk")
PwoP_all[1,] <- PwoP_bmr_al
PwoP_all[2,] <- PwoP_bmr_bl
PwoP_all[3,] <- PwoP_bmr_2015
PwoP_all[4,] <- PwoP_bmr_2016
PwoP_all[5,] <- PwoP_bmr_2017
PwoP_all[6,] <- PwoP_chePR
PwoP_all[7,] <- PwoP_cheCR
PwoP_all[8,] <- PwoP_ocq
PwoP_boot_all <- data.frame(matrix(nrow = 8,ncol = 2))
rownames(PwoP_boot_all) <- c("Black Mallard - above lake",
                        "Black Mallard - below lake",
                        "Black Mallard - 2015 cohort",
                        "Black Mallard - 2016 cohort",
                        "Black Mallard - 2017 cohort",
                        "Cheboygan - Pigeon Rd",
                        "Cheboygan - Cheboygan River",
                        "Ocqueoc")
colnames(PwoP_boot_all) <- c("CI95L","CI95U")
PwoP_boot_all[1,] <- bmr_al_conf
PwoP_boot_all[2,] <- bmr_bl_conf
PwoP_boot_all[3,] <- bmr_2015_conf
PwoP_boot_all[4,] <- bmr_2016_conf
PwoP_boot_all[5,] <- bmr_2017_conf
PwoP_boot_all[6,] <- chePR_conf
PwoP_boot_all[7,] <- cheCR_conf
PwoP_boot_all[8,] <- ocq_conf

PwoP_all <- cbind(PwoP_all,PwoP_boot_all)
write.table(PwoP_all,file = "Output/PwoP_summary_table.txt",append = F,quote = F,sep = "\t",row.names = T,col.names = T)

