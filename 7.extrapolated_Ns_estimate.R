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
locs <- read.table(file = "Input/locs_and_sublocations.txt",sep = "\t",header = T,stringsAsFactors = F)
ocq <- read.table(file = "Input/ocq_cohort_info.txt",sep = "\t",header = T,stringsAsFactors = F)
che <- read.table(file = "Input/che_cohort_info.txt",sep = "\t",header = T,stringsAsFactors = F)
bmr <- read.table(file = "Input/bmr_cohort_info.txt",sep = "\t",header = T,stringsAsFactors = F)
#run extrapolated Ns curves (function) for all cohorts
#BMR - Making files to subset cohorts ####
bmr <- bmr %>% select(ID,cohort)
colnames(bmr) <- c("OffspringID","cohort")
locs <- locs %>% rename(OffspringID = id)
bmr_colony <- merge(bmr_colony,locs)
bmr_colony <- merge(bmr_colony,bmr)
#BMR - above lake #####
bmr_col_al <- bmr_colony %>% 
  filter(sub.loc == "AboveLake")
bmr_al_Ns <- Ns_calc(family = bmr_col_al)
plot(bmr_al_Ns[[1]], ci.type = "poly", ci.col = "lightblue",ci.lty = 0,
     main = "Ns Accumulation Curve - Black Mallard 2019 Collection",
     ylim = c(0,20),
     xlab = "Number of Offspring sampled",
     ylab = "Number of Parents")
boxplot(bmr_al_Ns[[1]],add = T,pch = 19,cex = 0.2)
abline(h = bmr_al_Ns[[2]]$chao,col = "blue")
abline(h = bmr_al_Ns[[2]]$jack1,col = "darkgreen")
text(locator(1),"Chao estimate = 15.2",col = "blue",size = 10)
text(locator(1),"Jackknife estimate = 15.91",col = "darkgreen")

#BMR - below lake ####
bmr_col_bl <- bmr_colony %>% 
  filter(sub.loc == "BelowLake")
bmr_bl_Ns <- Ns_calc(family = bmr_col_bl)
plot(bmr_bl_Ns[[1]], ci.type = "poly", ci.col = "lightblue",ci.lty = 0,
     main = "Ns Accumulation Curve - Black Mallard 2015-2017 Cohorts",
     ylim = c(0,170),
     xlab = "Number of Offspring sampled",
     ylab = "Number of Parents")
boxplot(bmr_bl_Ns[[1]],add = T,pch = 19,cex = 0.2)
abline(h = bmr_bl_Ns[[2]]$chao,col = "blue")
abline(h = bmr_bl_Ns[[2]]$jack1,col = "darkgreen")
text(locator(1),"Chao estimate = 160.95",col = "blue")
text(locator(1),"Jackknife estimate = 139.97",col = "darkgreen")

#BMR - 2015 cohort ####
bmr15 <- bmr_colony %>% 
  filter(sub.loc == "BelowLake" & cohort == "2015")
bmr15_Ns <- Ns_calc(family = bmr15)

plot(bmr15_Ns[[1]], ci.type = "poly", ci.col = "lightblue",ci.lty = 0,
     main = "Ns Accumulation Curve - Black Mallard 2015 Cohort",
     ylim = c(0,180),
     xlab = "Number of Offspring sampled",
     ylab = "Number of Parents")
boxplot(bmr15_Ns[[1]],add = T,pch = 19,cex = 0.2)
abline(h = bmr15_Ns[[2]]$chao,col = "blue")
abline(h = bmr15_Ns[[2]]$jack1,col = "darkgreen")
text(locator(1),"Chao estimate = 165.99",col = "blue")
text(locator(1),"Jackknife estimate = 124.93",col = "darkgreen")
#BMR - 2016 cohort ####
bmr16 <- bmr_colony %>% 
  filter(sub.loc == "BelowLake" & cohort == "2016")
bmr16_Ns <- Ns_calc(family = bmr16)

plot(bmr16_Ns[[1]], ci.type = "poly", ci.col = "lightblue",ci.lty = 0,
     main = "Ns Accumulation Curve - Black Mallard 2016 Cohort",
     ylim = c(0,180),
     xlab = "Number of Offspring sampled",
     ylab = "Number of Parents")
boxplot(bmr16_Ns[[1]],add = T,pch = 19,cex = 0.2)
abline(h = bmr16_Ns[[2]]$chao,col = "blue")
abline(h = bmr16_Ns[[2]]$jack1,col = "darkgreen")
text(locator(1),"Chao estimate = 132.71",col = "blue")
text(locator(1),"Jackknife estimate = 117.96",col = "darkgreen")

#BMR - 2017 cohort ####
bmr17 <- bmr_colony %>% 
  filter(sub.loc == "BelowLake" & cohort == "2017")
bmr17_Ns <- Ns_calc(family = bmr17)

plot(bmr17_Ns[[1]], ci.type = "poly", ci.col = "lightblue",ci.lty = 0,
     main = "Ns Accumulation Curve - Black Mallard 2017 Cohort",
     pch = 19,
     ylim = c(0,60),
     xlab = "Number of Offspring sampled",
     ylab = "Number of Parents")
boxplot(bmr17_Ns[[1]],add = T,pch = 19,cex = 0.2)
abline(h = bmr17_Ns[[2]]$chao,col = "blue")
abline(h = bmr17_Ns[[2]]$jack1,col = "darkgreen")
text(locator(1),"Chao estimate = 48.05",col = "blue")
text(locator(1),"Jackknife estimate = 50.53",col = "darkgreen")

#CHE - Making files to subset cohorts ####
colnames(locs) <- c("OffspringID","loc","sub.loc")
che_colony <- merge(che_colony,locs)
che <- che %>% select(OffspringID,cohort)
che_colony <- merge(che_colony,che)

#CHE - Pigeon Rd - 2017 cohort####
che_PR17 <- che_colony %>% 
  filter(sub.loc == "Pigeon" & cohort == 2017)
chePR17_Ns <- Ns_calc(family = che_PR17)

plot(chePR17_Ns[[1]], ci.type = "poly", ci.col = "lightblue",ci.lty = 0,
     main = "Ns Accumulation Curve - Pigeon River 2017 Cohort",
     ylim = c(0,10),
     pch = 19,
     xlab = "Number of Offspring sampled",
     ylab = "Number of Parents")
boxplot(chePR17_Ns[[1]],add = T,pch = 19,cex = 0.2)
abline(h = chePR17_Ns[[2]]$chao,col = "blue")
abline(h = chePR17_Ns[[2]]$jack1,col = "darkgreen")
text(locator(1),"Chao estimate = 5",col = "blue")
text(locator(1),"Jackknife estimate = 5",col = "darkgreen")

#CHE - Cheboygan River - 2017 cohort####
cheCR17 <- che_colony %>%
  filter(sub.loc == "Cheboygan" & cohort == 2017)
cheCR17_Ns <- Ns_calc(family = cheCR17)

plot(cheCR17_Ns[[1]], ci.type = "poly", ci.col = "lightblue",ci.lty = 0,
     main = "Ns Accumulation Curve - Cheboygan River 2017 Cohort",
     ylim = c(0,30),
     pch = 19,
     xlab = "Number of Offspring sampled",
     ylab = "Number of Parents")
boxplot(cheCR17_Ns[[1]],add = T,pch = 19,cex = 0.2)
abline(h = cheCR17_Ns[[2]]$chao,col = "blue")
abline(h = cheCR17_Ns[[2]]$jack1,col = "darkgreen")
text(locator(1),"Chao estimate = 26.4",col = "blue")
text(locator(1),"Jackknife estimate = 19.2",col = "darkgreen")
#OCQ - making files for Ns calcs####
ocq_colony <- merge(ocq_colony,locs)
ocq <- ocq %>% select(OffspringID,cohort)
ocq_colony <- merge(ocq_colony,ocq)
#OCQ - 2015-2016 cohorts
ocq_Ns <- Ns_calc(family = ocq_colony)

plot(ocq_Ns[[1]], ci.type = "poly", ci.col = "lightblue",ci.lty = 0,
     main = "Ns Accumulation Curve - Ocqueoc River all cohorts",
     pch = 19,
     ylim = c(0,100),
     xlab = "Number of Offspring sampled",
     ylab = "Number of Parents")
boxplot(ocq_Ns[[1]],add = T,pch = 19,cex = 0.2)
abline(h = ocq_Ns[[2]]$chao,col = "blue")
abline(h = ocq_Ns[[2]]$jack1,col = "darkgreen")
text(locator(1),"Chao estimate = 73.72",col = "blue")
text(locator(1),"Jackknife estimate = 82.96",col = "darkgreen")

#OCQ - 2015 cohort ####
ocq15 <- ocq_colony %>% 
  filter(loc == "OCQ" & cohort == 2015)
ocq15_Ns <- Ns_calc(family = ocq15)

plot(ocq15_Ns[[1]], ci.type = "poly", ci.col = "lightblue",ci.lty = 0,
     main = "Ns Accumulation Curve - Ocqueoc River 2015 Cohort",
     ylim = c(0,100),
     pch = 19,
     xlab = "Number of Offspring sampled",
     ylab = "Number of Parents")
boxplot(ocq15_Ns[[1]],add = T,pch = 19,cex = 0.2)
abline(h = ocq15_Ns[[2]]$chao,col = "blue")
abline(h = ocq15_Ns[[2]]$jack1,col = "darkgreen")
text(locator(1),"Chao estimate = 95.65",col = "blue")
text(locator(1),"Jackknife estimate = 87.88",col = "darkgreen")
#OCQ - 2016 cohort ####
ocq16 <- ocq_colony %>% 
  filter(loc == "OCQ" & cohort == 2016)
ocq16_Ns <- Ns_calc(family = ocq16)

plot(ocq16_Ns[[1]], ci.type = "poly", ci.col = "lightblue",ci.lty = 0,
     main = "Ns Accumulation Curve - Ocqueoc River 2016 Cohort",
     ylim = c(0,80),
     pch = 19,
     xlab = "Number of Offspring sampled",
     ylab = "Number of Parents")
boxplot(ocq16_Ns[[1]],add = T,pch = 19,cex = 0.2)
abline(h = ocq16_Ns[[2]]$chao,col = "blue")
abline(h = ocq16_Ns[[2]]$jack1,col = "darkgreen")
text(locator(1),"Chao estimate = 62.95",col = "blue")
text(locator(1),"Jackknife estimate = 61.90",col = "darkgreen")


