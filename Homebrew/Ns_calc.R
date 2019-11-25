#Ns_calc - function
#function takes a colony bestconfig file and does an extrapolated Ns curve

#Inputs:
#family - best config file (has four columns: OffspringID, FatherID, MotherID, and ClusterIndex)
Ns_calc <- function(family){
  require(vegan)
  family$FatherID <- paste0("Dad",family$FatherID)
  family$MotherID <- paste0("Mom",family$MotherID)
  #making matrix for dads
  dads <- data.frame(matrix(0,nrow = length(family$OffspringID),ncol = length(unique(family$FatherID))))
  colnames(dads) <- unique(family$FatherID)
  rownames(dads) <- family$OffspringID
  #making matrix for moms
  moms <- data.frame(matrix(0,nrow = length(family$OffspringID),ncol = length(unique(family$MotherID))))
  colnames(moms) <- unique(family$MotherID)
  rownames(moms) <- family$OffspringID
  
  #loop to fill in matrix with parents
  for (i in 1:length(family$OffspringID)) {
    off <- family[i,]
    dadn <- which(off$FatherID == colnames(dads))
    dads[i,dadn] <- 1
  }
  for (i in 1:length(family$OffspringID)) {
    off <- family[i,]
    momn <- which(off$MotherID == colnames(moms))
    moms[i,momn] <- 1
  }
  parents <- cbind(moms,dads)
  ns_points <- specaccum(parents)
  asymp <- specpool(parents)
  output <- list(ns_points,asymp)
  output
}
