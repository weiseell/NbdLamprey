#Using COLONY parentage results to generate
#a minimum number of parents curve (Ns)
#Created by: Ellie Weise
#Originally Created on: Nov 12th, 2019
#Last Edited on: Nov 12th, 2019

#Goals:
#1.calculate Ns estimates for all cohorts
#2.plot Ns estimates
#3.estimate asymptote values?

#load libraries
library(vegan)
library(tidyverse)
#homebrew functions
source("Homebrew/Ns_calc.R")

#load in data
bmr_colony <- read.table("Input/colony.bestconfig.bmr.txt",header = T,sep = "\t",stringsAsFactors = F)
che_colony <- read.table("Input/colony.bestconfig.che.txt",header = T,sep = "\t",stringsAsFactors = F)
ocq_colony <- read.table("Input/colony.bestconfig.ocq.txt",header = T,sep = "\t",stringsAsFactors = F)
load(file = "Input/bmr_cohort_data.txt")
locs <- read.table(file = "Input/locs_and_sublocations.txt",sep = "\t",header = T,stringsAsFactors = F)
ocq <- read.table(file = "Input/ocq_cohort_info.txt",sep = "\t",header = T,stringsAsFactors = F)
che <- read.table(file = "Input/che_cohort_info.txt",sep = "\t",header = T,stringsAsFactors = F)
#run extrapolated Ns curves (function) for all cohorts
#BMR - Making files to subset cohorts ####
bmr <- bmr %>% mutate(ID = paste0(species,"_",loc,"_",num)) %>% select(ID,cohort)
colnames(bmr) <- c("OffspringID","cohort")
locs <- locs %>% rename(OffspringID = id)
bmr_colony <- merge(bmr_colony,locs)
bmr_colony <- merge(bmr_colony,bmr)
#BMR - above lake #####
bmr_col_al <- bmr_colony %>% 
  filter(sub.loc == "AboveLake")
bmr_al_Ns <- Ns_calc(family = bmr_col_al)

plot(bmr_al_Ns[[1]]$sites,bmr_al_Ns[[1]]$richness,
     main = "Ns Accumulation Curve - Black Mallard 2019 Collection",
     ylim = c(0,20),
     xlab = "Number of Offspring sampled",
     ylab = "Number of Parents")
abline(h = bmr_al_Ns[[2]]$chao,col = "blue")
abline(h = bmr_al_Ns[[2]]$jack1,col = "red")
#BMR - below lake ####
bmr_col_bl <- bmr_colony %>% 
  filter(sub.loc == "BelowLake")
bmr_bl_Ns <- Ns_calc(family = bmr_col_bl)
plot(bmr_bl_Ns[[1]]$sites,bmr_bl_Ns[[1]]$richness,
     main = "Ns Accumulation Curve - Black Mallard 2015-2017 Cohorts",
     ylim = c(0,180),
     xlab = "Number of Offspring sampled",
     ylab = "Number of Parents")
abline(h = bmr_bl_Ns[[2]]$chao,col = "blue")
abline(h = bmr_bl_Ns[[2]]$jack1,col = "red")

#BMR - 2015 cohort ####
bmr15 <- bmr_colony %>% 
  filter(sub.loc == "BelowLake" & cohort == "2015")
bmr15_Ns <- Ns_calc(family = bmr15)

plot(bmr15_Ns[[1]]$sites,bmr15_Ns[[1]]$richness,
     main = "Ns Accumulation Curve - Black Mallard 2015 Cohort",
     ylim = c(0,180),
     xlab = "Number of Offspring sampled",
     ylab = "Number of Parents")
abline(h = bmr15_Ns[[2]]$chao,col = "blue")
abline(h = bmr15_Ns[[2]]$jack1,col = "red")
#BMR - 2016 cohort ####
bmr16 <- bmr_colony %>% 
  filter(sub.loc == "BelowLake" & cohort == "2016")
bmr16_Ns <- Ns_calc(family = bmr16)

plot(bmr16_Ns[[1]]$sites,bmr16_Ns[[1]]$richness,
     main = "Ns Accumulation Curve - Black Mallard 2016 Cohort",
     ylim = c(0,180),
     xlab = "Number of Offspring sampled",
     ylab = "Number of Parents")
abline(h = bmr16_Ns[[2]]$chao,col = "blue")
abline(h = bmr16_Ns[[2]]$jack1,col = "red")

#BMR - 2017 cohort ####
bmr17 <- bmr_colony %>% 
  filter(sub.loc == "BelowLake" & cohort == "2017")
bmr17_Ns <- Ns_calc(family = bmr17)

plot(bmr17_Ns[[1]]$sites,bmr17_Ns[[1]]$richness,
     main = "Ns Accumulation Curve - Black Mallard 2017 Cohort",
     ylim = c(0,60),
     xlab = "Number of Offspring sampled",
     ylab = "Number of Parents")
abline(h = bmr17_Ns[[2]]$chao,col = "blue")
abline(h = bmr17_Ns[[2]]$jack1,col = "red")

#CHE - Making files to subset cohorts ####
colnames(locs) <- c("OffspringID","loc","sub.loc")
che_colony <- merge(che_colony,locs)
che <- che %>% select(ID,cohort) %>% rename(OffspringID = ID)
che_colony <- merge(che_colony,che)
#CHE - Pigeon Rd - 2016 cohort####
chePR16 <- che_colony %>% 
  filter(sub.loc == "Pigeon" & cohort == 2016)
chePR16_Ns <- Ns_calc(family = chePR16)

plot(chePR16_Ns[[1]]$sites,chePR16_Ns[[1]]$richness,
     main = "Ns Accumulation Curve - Pigeon River 2016 Cohort",
     ylim = c(0,20),
     xlab = "Number of Offspring sampled",
     ylab = "Number of Parents")
abline(h = chePR16_Ns[[2]]$chao,col = "blue")
abline(h = chePR16_Ns[[2]]$jack1,col = "red")
#CHE - Pigeon Rd - 2017 cohort####
che_PR17 <- che_colony %>% 
  filter(sub.loc == "Pigeon" & cohort == 2017)
chePR17_Ns <- Ns_calc(family = che_PR17)

plot(chePR17_Ns[[1]]$sites,chePR17_Ns[[1]]$richness,
     main = "Ns Accumulation Curve - Pigeon River 2017 Cohort",
     ylim = c(0,10),
     xlab = "Number of Offspring sampled",
     ylab = "Number of Parents")
abline(h = chePR17_Ns[[2]]$chao,col = "blue")
abline(h = chePR17_Ns[[2]]$jack1,col = "red")
#CHE - Cheboygan River - 2016 cohort####
cheCR16 <- che_colony %>%
  filter(sub.loc == "Cheboygan" & cohort == 2016)
cheCR16_Ns <- Ns_calc(family = cheCR16)

plot(cheCR16_Ns[[1]]$sites,cheCR16_Ns[[1]]$richness,
     main = "Ns Accumulation Curve - Cheboygan River 2016 Cohort",
     ylim = c(0,35),
     xlab = "Number of Offspring sampled",
     ylab = "Number of Parents")
abline(h = cheCR16_Ns[[2]]$chao,col = "blue")
abline(h = cheCR16_Ns[[2]]$jack1,col = "red")

#CHE - Cheboygan River - 2017 cohort####
cheCR17 <- che_colony %>%
  filter(sub.loc == "Cheboygan" & cohort == 2017)
cheCR17_Ns <- Ns_calc(family = cheCR17)

plot(cheCR17_Ns[[1]]$sites,cheCR17_Ns[[1]]$richness,
     main = "Ns Accumulation Curve - Cheboygan River 2017 Cohort",
     ylim = c(0,30),
     xlab = "Number of Offspring sampled",
     ylab = "Number of Parents")
abline(h = cheCR17_Ns[[2]]$chao,col = "blue")
abline(h = cheCR17_Ns[[2]]$jack1,col = "red")

#OCQ - making files for Ns calcs####
ocq_colony <- merge(ocq_colony,locs)
ocq <- ocq %>% select(ID,cohort) %>% rename(OffspringID = ID)
ocq_colony <- merge(ocq_colony,ocq)
#OCQ - 2015-2016 cohorts
ocq_Ns <- Ns_calc(family = ocq_colony)

plot(ocq_Ns[[1]]$sites,ocq_Ns[[1]]$richness,
     main = "Ns Accumulation Curve - Ocqueoc River all cohorts",
     ylim = c(0,100),
     xlab = "Number of Offspring sampled",
     ylab = "Number of Parents")
abline(h = ocq_Ns[[2]]$chao,col = "blue")
abline(h = ocq_Ns[[2]]$jack1,col = "red")
#OCQ - 2015 cohort ####
ocq15 <- ocq_colony %>% 
  filter(loc == "OCQ" & cohort == 2015)
ocq15_Ns <- Ns_calc(family = ocq15)

plot(ocq15_Ns[[1]]$sites,ocq15_Ns[[1]]$richness,
     main = "Ns Accumulation Curve - Ocqueoc River 2015 Cohort",
     ylim = c(0,100),
     xlab = "Number of Offspring sampled",
     ylab = "Number of Parents")
abline(h = ocq15_Ns[[2]]$chao,col = "blue")
abline(h = ocq15_Ns[[2]]$jack1,col = "red")

#OCQ - 2016 cohort ####
ocq16 <- ocq_colony %>% 
  filter(loc == "OCQ" & cohort == 2016)
ocq16_Ns <- Ns_calc(family = ocq16)

plot(ocq16_Ns[[1]]$sites,ocq16_Ns[[1]]$richness,
     main = "Ns Accumulation Curve - Ocqueoc River 2016 Cohort",
     ylim = c(0,80),
     xlab = "Number of Offspring sampled",
     ylab = "Number of Parents")
abline(h = ocq16_Ns[[2]]$chao,col = "blue")
abline(h = ocq16_Ns[[2]]$jack1,col = "red")



