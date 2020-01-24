#Nb and Ns Calculations
#Main objectives:
#1.Calculate Nb using the Parentage-without-Parents method
#2.Calculate Ns and extrapolated Ns

#libraries

#homebrew functions
source("Homebrew/PwoP.R")
source("Homebrew/PwoP_boot.R")
source("Homebrew/Ns_calc.R")
#load in data
load("Software_outputs/test_best_config.rda")

#calculating PwoP
Nb_PwoP <- PwoP(best_config)
uncert <- PwoP_boot(best_config,iter = 1000,alpha = 0.05,real_Nb = Nb_PwoP$Nb)

#calculating Ns and extrapolated Ns
#getting rid of duplicates (not sure why they're there, but they are :/)
best_config2 <- best_config[!duplicated(best_config$OffspringID),]
best_config2 <- na.omit(best_config2)
Ns <- Ns_calc(best_config2)

#Figure for Ns
#adding accumulation curve with uncertainty
plot(Ns[[1]], ci.type = "poly", ci.col = "lightblue",ci.lty = 0,
     main = "Ns Accumulation Curve - Test set",
     ylim = c(0,100),
     xlab = "Number of Offspring sampled",
     ylab = "Number of Parents")
#adding boxplot distribution
boxplot(Ns[[1]],add = T,pch = 19,cex = 0.2)
#adding horizontal lines for the extrapolated estimates
abline(h = Ns[[2]]$chao,col = "blue")
abline(h = Ns[[2]]$jack1,col = "darkgreen")
#adds labels to the extrapolation lines
#!# need to click on the graph to add the label
#running the line will initiate a plus sign on your cursor that you can use to place the label
text(locator(1),"Chao estimate = 85.8",col = "blue")
text(locator(1),"Jackknife estimate = 91.9",col = "darkgreen")

#!# Need to make an Nb/Ns table and save plot (also need to do this for script 4)