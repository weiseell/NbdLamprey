#LD_filter - function
#Goal: select a set of SNPs so that there is one SNP per tag
#Good for SNP sets for NeEstimator/LDNe

#Inputs:
#df: a dataframe containing targeted SNP names and percent genotyped

LD_filter <- function(df){
  require(tidyverse)
  tags <- unique(df$target)
  #making a matrix for the SNPs full of NAs to be filled in
  SNPs  <- data.frame(matrix(data = NA,nrow = length(tags),ncol = ncol(df)))
  rownames(SNPs) <- tags
  colnames(SNPs) <- colnames(df)
  #loop to select SNPs and put them into empty dataframe
  i <- 1
  for (i in 1:length(tags)) {
    tmp_tag <- tags[i]
    tmp <- subset(df,df$target == tmp_tag)
    SNP_pick <- tmp %>% 
      mutate(rvar = runif(n = nrow(tmp))) %>% 
      arrange(desc(pGT),rvar) %>% 
      slice(1)
    SNP_pick$rvar <- NULL
    SNPs[i,] <- SNP_pick
  }
  SNPs
}
