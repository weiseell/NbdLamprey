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

COLONY_filter <- function(df,window = 1000000,pGT_min = 0.8,MAF_min = 0.05){
  require(tidyverse)
  df <- df %>% 
    filter(pGT > pGT_min & MAF > MAF_min)
  #filtering step
  loci_select <- df %>% 
    arrange(CHROM,POS)
  chroms <- unique(loci_select$CHROM)
  i <- 1
  COL_SNPs <- as.data.frame(matrix(nrow = 0,ncol = ncol(loci_select)),stringsAsFactors = F)
  colnames(COL_SNPs) <- colnames(loci_select)
  for (i in 1:length(chroms)) {
    #subset to 1 scaffold for selection
    ch <- chroms[i]
    SNPs <- subset(loci_select,loci_select$CHROM == ch)
    #getting the max number of SNPs to select 
    j <- as.numeric(SNPs$POS[1])
    while (j <= max(as.numeric(SNPs$POS))) {
      tmpSNP <- SNPs[min(which(as.numeric(SNPs$POS) >= j)),]
      COL_SNPs <- rbind(COL_SNPs,tmpSNP)
      j <- as.numeric(tmpSNP$POS) + window
    }
  }
  COL_SNPs
}






