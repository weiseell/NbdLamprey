#Using COLONY parentage results to generate
#a minimum number of parents curve (Ns)
#Test data and code for Lilian
#Created by: Ellie Weise
#Originally Created on: Nov 12th, 2019
#Last Edited on: Dec 12th, 2019

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
#run extrapolated Ns curves (function)
#BMR - above lake #####
#run extrapolated Ns curves (function)
bmr_al_Ns <- Ns_calc(family = bmr_col_al)
#setting up plot and uncertainty curve
plot(bmr_al_Ns[[1]], ci.type = "poly", ci.col = "lightblue",ci.lty = 0,
     main = "Ns Accumulation Curve - Black Mallard 2019 Collection",
     ylim = c(0,20),
     xlab = "Number of Offspring sampled",
     ylab = "Number of Parents")
#adding boxplot
boxplot(bmr_al_Ns[[1]],add = T,pch = 19,cex = 0.2)
#adding straight lines for my extrapolated estimates
abline(h = bmr_al_Ns[[2]]$chao,col = "blue")
abline(h = bmr_al_Ns[[2]]$jack1,col = "darkgreen")
#adding custom text to indicate the extrapolated estimate values
text(locator(1),"Chao estimate = 15.2",col = "blue",size = 10)
text(locator(1),"Jackknife estimate = 15.91",col = "darkgreen")

