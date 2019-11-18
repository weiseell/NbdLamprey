#Ns_calc - function
#function takes a colony bestconfig file and does an extrapolated Ns curve

#Inputs:
#family - best config file (has four columns: OffspringID, FatherID, MotherID, and ClusterIndex)
#step - integer for the step size of the subset size (default is 50)
#reps - integer for the number of replicates at each step (default is 100)
Ns_calc <- function(family,step = 50,reps = 100){
  #selecting the sample size of individuals
  steps <- seq(from = step,to = length(family$OffspringID),by = step)
  #empty matrix to store parents
  min_parents <- data.frame(matrix(NA,nrow = length(steps),ncol = 2))
  colnames(min_parents) <- c("noff","npar")
  i <- 1
  j <- 1
  for (i in 1:length(steps)) {
    size <- steps[i]
    parents <- vector()
    for (j in 1:reps) {
      #randomly choosing offspring
      tmp <- sample.int(length(family$OffspringID),size = size)
      df <- family[tmp,]
      #getting the number of parents for the subsetted offspring
      nparents <- length(unique(df$MotherID))+length(unique(df$FatherID))
      parents <- append(parents,nparents)
    }
    min_parents$noff[i] <- size
    #getting the mean number of parents for each subset size
    min_parents$npar[i] <- mean(parents)
  }
  min_parents
}


