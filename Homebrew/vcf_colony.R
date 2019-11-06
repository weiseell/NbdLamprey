#vcf_colony - function
#goal: convert a genotype file in vcf format to a COLONY type format

#Inputs:
#vcf = data frame, genotype file with generic SNP names in vcf format
#empty = logical, if T runs additional section that eliminates individuals with no genotypes
#colony cannot tolerate individuals with no genotypes and will error out

vcf_colony <- function(vcf,empty = T){
  #sorting data for COLONY and reformatting genotypes
  vcf$gt <- gsub(pattern = "1/1",replacement = "2/2",x = vcf$gt)
  vcf$gt <- gsub(pattern = "0/0",replacement = "1/1",x = vcf$gt)
  vcf$gt <- gsub(pattern = "0/1",replacement = "1/2",x = vcf$gt)
  vcf$gt <- gsub(pattern = "\\./\\.",replacement = "0/0",x = vcf$gt)
  
  #separating genotypes into alleles and then merging them together
  #allele one
  a1 <- vcf %>% 
    separate(col = gt,into = c("a","b"),sep = "/") %>%
    select(-b) %>%
    spread(key = "locus",value = "a")
  colnames(a1) <- paste0(colnames(a1),".1")
  names(a1)[1] <- "id"
  
  #allele two
  a2 <- vcf %>% 
    separate(col = gt,into = c("a","b"),sep = "/") %>%
    select(-a) %>%
    spread(key = "locus",value = "b")
  colnames(a2) <- paste0(colnames(a2),".2")
  names(a2)[1] <- "id"
  
  #merging on IDs
  colony <- merge(a1,a2,by="id") 
  colony <- colony %>% select(sort(colnames(colony)))
  
  if(empty == T){
    colony1 <- colony %>% 
      gather(key = locus,value = gt, -id)
    colony1 <- colony1 %>% 
      mutate(gtcheck = ifelse(colony1$gt == "0",0,1)) %>% 
      group_by(id) %>% 
      summarise(sums = sum(gtcheck),
                empty = ifelse(sums == 0,T,F))
    
    colony$empty <- colony1$empty
    colony <- colony %>% 
      filter(empty == F) %>% 
      select(-empty)
    
  }
  colony
}
