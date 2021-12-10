#Nb and Ns estimates
#1. Calculate Nb - PwoP method and both Ns estimates using reconstructed pedigree
#2. Extract Nb - LD method from NeEstimator tabular output
#3. Extract Nb - SF method from Colony Ne output

#load libraries
library(tidyverse)

#load functions
source("Homebrew/PwoP.R")
source("Homebrew/PwoP_uncert.R")
source("Homebrew/Ns_calc.R")

##1. Calculate Nb - PwoP method and both Ns estimates using reconstructed pedigree
pops <- c("BMR15","BMR16","bmrAL","chePR","OCQ")
i <- 1
Nb_PwoP <- data.frame(matrix(nrow = length(pops),ncol = 7))
colnames(Nb_PwoP) <- c("Pop","SampSize","PwoP_Nb","kbar","Vk", "PwoP_LCI", "PwoP_HCI")
Ns_all <- data.frame(matrix(nrow = length(pops),ncol = 6))
colnames(Ns_all) <- c("Pop","Ns" ,"Ns_Chao","Chao_uncert","Ns_Jackknife","Jackknife_uncert")

#calculation loop
for (i in 1:length(pops)) {
  print(i)
  #read in file
  df <- readLines(paste0("Software_Outputs/",pops[i],"_Output.data.BestCluster"))
  #separate file into usable data frame
  df <- strsplit(df,"\\s+")
  df1 <- matrix(unlist(df),byrow = T)
  df1 <- df1[df1 != ""]
  df1 <- matrix(df1,ncol = 5,byrow = T)
  df1 <- as.data.frame(df1)
  colnames(df1) <- df1[1,]
  df1 <- df1[-1,]
  
  ##save bestconfig files for future use in figures
  write.table(df1,paste0("Software_outputs/",pops[i],"BestCluster_formatted.txt"),append = F,quote = F,sep = "\t")
  
  #read in configarchive for PwoP uncertainty
  ca_tmp <- readLines(paste0("Software_Outputs/",pops[i],"_Output.data.ConfigArchive"))
  
  #calculate Nb_PwoP
  PwoP_tmp <- PwoP(df1)
  uncert <- PwoP_uncert(ca = ca_tmp,bc = df1)
  names(uncert) <- c("CI_Lower","CI_Upper")
  
  #combine outputs
  tmp <- c(pops[i],length(df1$OffspringID),PwoP_tmp,uncert)
  Nb_PwoP[i,] <- tmp
  
  #calculate extrapolated Ns
  Ns_tmp <- Ns_calc(df1)
  Ns_all[i,] <- c(pops[i],Ns_tmp[[3]],Ns_tmp[[2]]$chao,Ns_tmp[[2]]$chao.se,Ns_tmp[[2]]$jack1,Ns_tmp[[2]]$jack1.se)
}


##2. Extract Nb - LD method from NeEstimator tabular output
#!# LD table is concatinated into excel and the header was removed
#!# to ease the analysis
Nb_LD <- read.table("Software_Outputs/NeEstimator_results_063021.txt",header = T)
Nb_LD1 <- Nb_LD %>% 
  select(Population,Nb,JackLCI,JackHCI) %>% 
  rename(Pop=Population,LD_Nb=Nb,LD_LCI=JackLCI,LD_HCI=JackHCI)
#3. Extract Nb - SF method from Colony Ne output
i <- 1
Nb_SF <- data.frame(matrix(nrow = length(pops),ncol = 4))
colnames(Nb_SF) <- c("SF_Nb","SF_LCI","SF_HCI","Pop")
for (i in 1:length(pops)) {
  df <- readLines(paste0("Software_Outputs/",pops[i],"_Output.data.Ne"))
  Ne <- grep(pattern = "Ne",x = df,value = T)
  Ci <- grep(pattern = "CI95",x = df, value = T)
  
  tmp <- data.frame(Nb=gsub(" ","",Ne[1]),CILow=gsub(" ","",Ci[1]),CIHigh=gsub(" ","",Ci[2]))
  tmp1 <- tmp %>% 
    gather(key = "stat",value = "value") %>% 
    separate(value,into = c("tmp","value"),sep = "=") %>% 
    select(-tmp) %>% 
    spread(key = stat,value = value) %>% 
    select(Nb,CILow,CIHigh)
  tmp1$Pop <- pops[i]
  Nb_SF[i,] <- tmp1
}


##combining all estimates into a table, along with Vk and kbar
Nb_Ns <- Nb_LD1 %>% 
  full_join(Nb_PwoP,by="Pop") %>% 
  full_join(Nb_SF,by="Pop") %>% 
  full_join(Ns_all,by="Pop")

write.table(Nb_Ns,file = "Summary_Stats/genetic.estimates.txt",append = F,quote = F,sep = "\t")
