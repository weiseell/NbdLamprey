#ped_filter - function
#Goal: select SNPs based on a window size to prevent excess linked loci

#ped_filter
#select SNPs based on a window size to prevent excess linked loci
#good for Ne estimates where excess linked loci would create bias
#default window is 1MB, which is a conservative number to ensure independent loci
#can change based on recombination rates for specific species

#Inputs:
#df = file with the following SNP summary statistics:
#data is already filtered by heterozygosity and sequence depth
#window = integer - size of average separation between SNPs in bases
#pGT_min = integer - minimum number for percent individuals genotyped required for a SNP to be included
#otherwise, another one is selected

ped_filter <- function(df,window = 1000000,pGT_min = 0.8){
  require(tidyverse)
  df <- df %>% 
    filter(pGT > pGT_min)
  #calculating how many loci per scaffold to select
  df$nloci1 <- as.numeric(df$POS)/window
  #sorting loci into genome sections
  df$nloci2 <- floor(df$nloci1)
  df$nloci3 <- abs((df$nloci1 - df$nloci2)-0.5)
  #calculating the maximum number of SNPs that will be selected in this filtering
  df %>% 
    group_by(CHROM,nloci2) %>% 
    count() %>% 
    group_by(CHROM) %>% 
    count() %>% 
    pull(n) %>% 
    sum()
  
  #filtering step
  #not sure about the order of these steps, or the number of loci I end up with
  loci_select <- df %>% 
    #creating a random variable as a tie breaker
    mutate(rvar = runif(n = nrow(df))) %>%
    #sorting loci
    group_by(CHROM,nloci2) %>%
            # 1. Genome "section" (determined above)
    arrange(nloci2,
            # 2. Highest MAF
            desc(MAF),
            # 3. highest percent genotyped for all populations
            desc(pGT),
            # 4. Closest to the center of the genome section
            nloci3,
            # 5. random variable to break ties
            rvar) %>%
    #taking the best loci for each section
    slice(1)
  
  
  loci_select
  
}
