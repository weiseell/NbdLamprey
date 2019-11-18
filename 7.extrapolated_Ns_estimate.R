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
#homebrew functions
source("Homebrew/Ns_calc.R")

#load in data
bmr_colony <- read.table("Input/colony.bestconfig.bmr.txt",header = T,sep = "\t",stringsAsFactors = F)
che_colony <- read.table("Input/colony.bestconfig.che.txt",header = T,sep = "\t",stringsAsFactors = F)
ocq_colony <- read.table("Input/colony.bestconfig.ocq.txt",header = T,sep = "\t",stringsAsFactors = F)

#run extrapolated Ns curves (function) for all cohorts
bmr_ns <- Ns_calc(family = bmr_colony,step = 10,reps = 100)
che_ns <- Ns_calc(family = che_colony,step = 1,reps = 100)
ocq_ns <- Ns_calc(family = ocq_colony,step = 10,reps = 100)
#plot Ns curves
plot(bmr_ns$noff,bmr_ns$npar,
     main = "Black Mallard River",
     xlab = "Number of Offspring",
     ylab = "Number of Parents")
plot(che_ns$noff,che_ns$npar,
     main = "Cheboygan River",
     xlab = "Number of Offspring",
     ylab = "Number of Parents")

plot(ocq_ns$noff,ocq_ns$npar,
     main = "Ocqueoc River",
     xlab = "Number of Offspring",
     ylab = "Number of Parents")
#fitting best fit lines to determine asymptotes for minimum number of parents
#black mallard
fit_bmr <- lm(bmr_ns$npar~log(bmr_ns$noff))
summary(fit_bmr)
bmr1 <- sortedXyData(bmr_ns$noff,bmr_ns$npar)
NLSstRtAsymptote(bmr1)

#cheboygan
fit_che <- lm(che_ns$npar~log(che_ns$noff))
summary(fit_che)
che1 <- sortedXyData(che_ns$noff,che_ns$npar)
NLSstRtAsymptote(che1)

#ocqueoc
fit_ocq <- lm(ocq_ns$npar~log(ocq_ns$noff))
summary(fit_ocq)
ocq1 <- sortedXyData(ocq_ns$noff,ocq_ns$npar)
NLSstRtAsymptote(ocq1)

#using vegan
specaccum(bmr_colony)



