#separate config archive in a loop and calculate PwoP on all of them
PwoP_uncert <- function(ca,bc){
  require(tidyverse)
  #read in source
  source("Homebrew/PwoP.R")
  #create Config Archive list
  Configs <- (grep(pattern = "Config#",x = ca))
  i <- 1
  PwoP_ca <- data.frame(matrix(data = NA, nrow = length(Configs), ncol = 3))
  colnames(PwoP_ca) <- c("Nb","kbar","Vk")
  
  for (i in 1:length(Configs)) {
    if(i < length(Configs)){
      tmp <- ca[Configs[i]:(Configs[i+1]-1)]
    }
    else{
      tmp <- ca[Configs[i]:length(ca)]
    }
    likeli <- tmp[1]
    header <- tmp[2]
    header <- strsplit(x = header,split = ",")[[1]]
    tmp <- tmp[-(1:2)]
    tmp <- data.frame(temp = tmp)
    tmp <- tmp %>% 
      separate(col = temp, into = header,sep = ",")
    PwoP_ca[i,] <- PwoP(tmp)
  }
  
  #simulate data and calculate PwoP on those
  PwoP_sim <- data.frame(matrix(data = NA, nrow = 1000, ncol = 3))
  colnames(PwoP_sim) <- c("Nb","kbar","Vk")
  i <- 1
  trueNb <- PwoP(bc)
  
  for (i in 1:1000) {
    npar <- round(trueNb["Nb"]/2)
    family <- data.frame(matrix(data = NA, nrow = length(bc$OffspringID), ncol = 3))
    names(family) <- c("OffspringID","FatherID","MotherID")
    family$OffspringID <- seq(1,length(bc$OffspringID))
    family$FatherID <- sample(seq(1,npar),size = length(bc$OffspringID),replace = T)
    family$MotherID <- sample(seq(1,npar),size = length(bc$OffspringID),replace = T)
    PwoP_sim[i,] <- PwoP(family)
  }
  
  #combine Vc and Vg to get total variance
  Vg <- var(1/(2*PwoP_sim$Nb))
  Vs <- var(1/(2*PwoP_ca$Nb))
  LCIrecip <- 1/(2*trueNb["Nb"]) + 1.96*sqrt(Vg+Vs)
  LCI <- 0.5/LCIrecip
  UCIrecip <- 1/(2*trueNb["Nb"]) - 1.96*sqrt(Vg+Vs) 
  UCI <- 0.5/UCIrecip  
  
  conf_int <- c(LCI,UCI)
  names(conf_int) <- c("LCI","UCI")
  conf_int
}








