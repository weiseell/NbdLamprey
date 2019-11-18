#Calculating Nb using the Parentage without Parents (PwoP) method (Waples and Waples 2011)
#Created by: Ellie Weise
#Originally Created on: Nov 12th, 2019
#Last Edited on: Nov 12th, 2019

#Goals:
#1.calculate Ns estimates for all cohorts
#2.plot Ns estimates
#3.estimate asymptote values?

#homebrew functions
source("Homebrew/PwoP.R")
source("Homebrew/PwoP_boot.R")
#load in data
bmr_colony <- read.table("Input/colony.bestconfig.bmr.txt",header = T,sep = "\t",stringsAsFactors = F)
che_colony <- read.table("Input/colony.bestconfig.che.txt",header = T,sep = "\t",stringsAsFactors = F)
ocq_colony <- read.table("Input/colony.bestconfig.ocq.txt",header = T,sep = "\t",stringsAsFactors = F)

#calculating Nb
bmr_pwop <- PwoP(bmr_colony)
che_pwop <- PwoP(che_colony)
ocq_pwop <- PwoP(ocq_colony)
#using a bootstrapping method to construct confidence intervals for Nb estimates
bmr_conf <- PwoP_boot(family = bmr_colony,iter = 1000,alpha = 0.05,real_Nb = bmr_pwop)
che_conf <- PwoP_boot(family = che_colony,iter = 1000,alpha = 0.05,real_Nb = che_pwop)
ocq_conf <- PwoP_boot(family = ocq_colony,iter = 1000,alpha = 0.05,real_Nb = ocq_pwop)

