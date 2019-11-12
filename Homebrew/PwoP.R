#PwoP - function
#calculates Nb using a method that uses variance in reproductive success to estimate Nb
#called parentage without parents (PwoP) - Waples and Waples 2011

#Inputs
#df - data frame with four columns: OffspringID, FatherID, MotherID, and ClusterIndex

PwoP <- function(family){
  #turning parent names into character strings
  family$FatherID <- paste0("Dad",family$FatherID)
  family$MotherID <- paste0("Mom",family$MotherID)
  #making k matrix - number of offspring for each parent
  dadcounts <- as.vector(table(family$FatherID))
  momcounts <- as.vector(table(family$MotherID))
  moms <- data.frame(parent = unique(family$MotherID),k = momcounts,stringsAsFactors = F)
  dads <- data.frame(parent = unique(family$FatherID),k = dadcounts,stringsAsFactors = F)
  parents <- rbind(moms,dads)
  
  #calculating Nb from parent counts
  parents$k2 <- parents$k^2
  Nb <- (sum(parents$k)-1)/((sum(parents$k2)/sum(parents$k))-1)
  Nb
}






