#PwoP_boot - function
#bootstrapping for Parentage without Parents Nb Method
#!# Not used as of August 2021 - replaced by PwoP_uncert

#Inputs:
#family - data frame with four columns: OffspringID, FatherID, MotherID, and ClusterIndex
#iter - the number of iterations used in the bootstrapping algorithm
#real_Nb - numeric generated from empirical estimate, used to generate CI
PwoP_boot <- function(family,iter,alpha,real_Nb){
  n_fam <- length(family$FatherID)
  boot_Nb <- numeric(iter)
  i <- 1
  for (i in 1:iter) {
    tmp <- family[sample(n_fam,n_fam,replace = T),]
    tmp1 <- PwoP(tmp)
    boot_Nb[i] <- tmp1["Nb"]
  }
  
  #constructing confidence intervals for the bootstrapped values
  se_boot_Nb <- sd(boot_Nb)/sqrt(n_fam)
  cv_boot_Nb <- qnorm(alpha/2,lower.tail=FALSE)
  
  conf_int <- real_Nb+c(-1,1)*se_boot_Nb*cv_boot_Nb
  conf_int
}







