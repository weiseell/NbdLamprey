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
  dadcounts <- table(family$FatherID)
  momcounts <- table(family$MotherID)
  moms <- data.frame(parent = names(momcounts),k = as.vector(momcounts),stringsAsFactors = F)
  dads <- data.frame(parent = names(dadcounts),k = as.vector(dadcounts),stringsAsFactors = F)
  parents <- rbind(moms,dads)
  
  #calculating Nb from parent counts
  parents$k2 <- parents$k^2
  N <- length(unique(parents$parent))
  Vk <- (sum(parents$k2)/N)-((sum(parents$k)/N)^2)
  kbar <- mean(parents$k)
  Nb <- (sum(parents$k)-1)/((sum(parents$k2)/sum(parents$k))-1)
  out <- c(Nb,kbar,Vk)
  names(out) <- c("Nb","kbar","Vk")
  out
}






