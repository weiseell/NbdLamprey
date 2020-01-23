#marker_create - function
#goal: create a markers file for colony2.dat input file

#Inputs:
#SNPs = factor, names of all loci
#cod = the number of co-dominant markers allowed within COLONY
#gte = level of genotyping error that's accounted for during COLONY run
#ote = level of random error that's accounted for during COLONY run

marker_create <- function(SNPs,cod = 0,gte = 0.02,ote = 0.001){
  #making a matrix of zeros of the correct size and adding colnames
  nmarkers <- (length(SNPs)-1)/2
  markers  <- data.frame(matrix(data = 0,nrow = 3,ncol = nmarkers))
  colnames(markers) <- SNPs[seq(from=2, to = nmarkers*2,by = 2)]
  
  #filling in matrix
  markers[1,] <- cod #for co-dominant markers
  markers[2,] <- gte # genotyping error
  markers[3,] <- ote #other types of error
  
  #output
  markers
}

