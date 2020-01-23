#vcf_genepop
#function to convert a genotyped SNP file from vcf format to GENEPOP format

#Inputs:
#vcf = data frame, genotype file with generic loci names in vcf format
#three columns: Indiv,SNP, and gt
#empty = logical, if T runs additional section that eliminates individuals with no genotypes

vcf_genepop <- function(vcf,empty = T){
  require(tidyverse)
  vcf$gt <- gsub(pattern = "1/1",replacement = "0202",x = vcf$gt)
  vcf$gt <- gsub(pattern = "0/0",replacement = "0101",x = vcf$gt)
  vcf$gt <- gsub(pattern = "0/1",replacement = "0102",x = vcf$gt)
  vcf$gt <- gsub(pattern = "\\./\\.",replacement = "0000",x = vcf$gt)
  
  genepop <- vcf %>%
    spread(key = SNP,value = gt)
  if(empty == T){
    genepop1 <- genepop %>% 
      gather(key = locus,value = gt, -indiv)
    genepop1 <- genepop1 %>% 
      mutate(gtcheck = ifelse(genepop1$gt == "0",0,1)) %>% 
      group_by(indiv) %>% 
      summarise(sums = sum(gtcheck),
                empty = ifelse(sums == 0,T,F))
    
    genepop$empty <- genepop1$empty
    genepop <- genepop %>% 
      filter(empty == F) %>% 
      select(-empty)
  }
}
