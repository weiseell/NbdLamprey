#ld_filter
#select SNPs based on a window size to prevent excess linked loci
#good for Ne estimates where excess linked loci would create bias
#default window is 1MB, which is a conservative number to ensure independent loci
#can change based on recombination rates for specific species

#Inputs:
#summ = file with SNP summary statistics
#gt = data frame of SNP genotype calls
#mine is already filtered by heterozygosity and sequence depth
#window = integer - size of average separation between SNPs in bases
#pGT_min = integer - minimum number for percent individuals genotyped required for a SNP to be included
#otherwise, another one is selected

ld_filter <- function(summ,gt,window = 1000000,pGT_min = 0.8){
  require(tidyverse)
  summ <- summ %>% 
    filter(pGT > pGT_min)
  #calculating how many loci per scaffold to select
  summ$nloci1 <- summ$POS/window
  summ$nloci2 <- floor(summ$nloci1)
  summ$nloci3 <- abs((summ$nloci1 - summ$nloci2)-0.5)
  hist(summ$nloci3)
  summ %>% 
    group_by(CHROM,nloci2) %>% 
    count() %>% 
    group_by(CHROM) %>% 
    count() %>% 
    pull(n) %>% 
    sum()
  
  #filtering step
  #not sure about the order of these steps, or the number of loci I end up with
  loci_select <- summ %>% 
    mutate(rvar = runif(n = nrow(summ))) %>%
    group_by(CHROM,nloci2) %>%
    arrange(nloci2,desc(MAF),desc(pGT),nloci3,rvar) %>%
    slice(1)
  hist(loci_select$MAF)
  hist(loci_select$pGT)
  
  loci_select1 <- merge(loci_select,geno)
  loci_select1 <- loci_select1 %>% 
    select(-het:-MAF,-nloci1:-rvar)
  
  loci_select1
}
