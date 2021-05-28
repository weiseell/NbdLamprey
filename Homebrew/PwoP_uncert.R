#separate config archive in a loop and calculate PwoP on all of them
PwoP_uncert <- function(ca,bc){
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
  PwoP_sim <- data.frame(matrix(data = NA, nrow = 100, ncol = 3))
  colnames(PwoP_sim) <- c("Nb","kbar","Vk")
  i <- 1
  trueNb <- PwoP(bc)
  
  for (i in 1:length(Configs)) {
    npar <- round(trueNb["Nb"]/2)
    family <- data.frame(matrix(data = NA, nrow = length(bc$OffspringID), ncol = 3))
    names(family) <- c("OffspringID","FatherID","MotherID")
    family$OffspringID <- seq(1,length(bc$OffspringID))
    family$FatherID <- sample(seq(1,npar),size = length(bc$OffspringID),replace = T)
    family$MotherID <- sample(seq(1,npar),size = length(bc$OffspringID),replace = T)
    PwoP_sim[i,] <- PwoP(family)
  }
  
  #combine Vc and Vg to get total variance
  Vg <- var(PwoP_ca$Nb)
  Vc <- var(PwoP_sim$Nb)
  UCI <- (1/(2*trueNb["Nb"]))+1.96*sqrt(Vg+Vc)
  UCI <- trueNb["Nb"]+(1/UCI)/2
  LCI <- (1/(2*trueNb["Nb"]))-1.96*sqrt(Vg+Vc)
  LCI <- trueNb["Nb"]+(1/LCI)/2
  
  conf_int <- c(LCI,UCI)
  names(conf_int) <- c("LCI","UCI")
  conf_int
}








