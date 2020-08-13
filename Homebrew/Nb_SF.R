#Nb_SF <- function
#calculates Nb using the sibship method based on Colony output files

Nb_SF <- function(family){
  family$parent_ID <- paste0(family$FatherID,"_",family$MotherID)
  
  #counting the number of full siblings
  fullsib_counts <- 0
  for (i in 1:length(family$OffspringID)) {
    tmp_counts <- length(which(family$parent_ID[i] == family$parent_ID))
    fullsib_counts <- fullsib_counts + tmp_counts
  }
  
  #counting the number of maternal half siblings
  momhalfsib_counts <- 0
  for (i in 1:length(family$OffspringID)) {
    tmp_counts <- length(which(family$MotherID[i] == family$MotherID & family$FatherID[i] != family$FatherID))
    momhalfsib_counts <- momhalfsib_counts + tmp_counts
  }
  
  #counting the number of paternal half siblings
  dadhalfsib_counts <- 0
  for (i in 1:length(family$OffspringID)) {
    tmp_counts <- length(which(family$FatherID[i] == family$FatherID & family$MotherID[i] != family$MotherID))
    dadhalfsib_counts <- dadhalfsib_counts + tmp_counts
  }
  #getting total number of pairwise comparisons
  numcomparisons <- length(family$OffspringID)^2 - length(family$OffspringID)
  
  #calculating probabilities that a pair of offspring are half-sib/full-sib
  Q1 <- dadhalfsib_counts/numcomparisons
  Q2 <- momhalfsib_counts/numcomparisons
  Q3 <- fullsib_counts/numcomparisons
  
  #getting the number of dads and moms
  N1 <- length(unique(family$FatherID))
  N2 <- length(unique(family$MotherID))
  
  alpha <- 0
  
  Nb <- 1/(((1+3*alpha)/4)*(Q1+Q2+Q3) - (alpha/2)*((1/N1)+(1/N2)))
  Nb
}










