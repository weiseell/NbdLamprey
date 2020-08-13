#subsetting full data set to make a test data set for the nbd lamprey package
#test data will include a randomly selected set of 20,000  SNPs and 200 Black Mallard individuals

#randomly selecting rows and columns
rows <- length(geno$CHROM)
cols <- seq(451,1538,by = 1)
samrows <- sample(rows,20000,replace = F)
samcols <- sample(cols,200,replace = F)
ids <- geno[,1:2]

#getting corresponding gdepth and allele frequencies
subgeno <- subgeno %>% mutate(ID = paste(CHROM,POS,sep = "_")) %>% 
  select(ID,everything())
subaf <- subaf %>% mutate(ID = paste(Chrom_Pos,Pos,sep = "_")) %>% 
  select(ID,everything())

#subsampling genotypes
subgeno <- geno[samrows,]
subgeno1 <- subgeno %>% select(ID,samcols)

#subsampling AF
subaf <- af[samrows,]

gdepth <- gdepth1 %>% spread(key = indiv,value = depth)
i <- 1
checks <- matrix(data = NA,nrow = length(subgeno$CHROM),ncol = 1)
for (i in 1:length(subgeno$PM_BMR_120)) {
  checks[i,] <- which(gdepth$ID == subgeno$ID[i])
}
subdepth <- gdepth[checks,]

subcols <- colnames(subgeno1)
subdepth1 <- subdepth %>% select(subcols)

save(subgeno1,file = "Input/test_genotypes.rda")
save(subaf,file = "Input/test_MAF.rda")
save(subdepth1, file = "Input/test_depths.rda")


#selecting mixture analysis and COLONY output data
indivs <- colnames(gt)
indivs <- indivs[-1]

checks <- matrix(data = NA,nrow = length(indivs),ncol = 1)
i <- 1
for (i in 1:length(indivs)) {
  tmp <- which(bmr_colony$OffspringID == indivs[i])
  if(length(tmp)>0){
    checks[i] <- tmp
  }
}
best_config <- bmr_colony[checks,]

checks <- matrix(data = NA,nrow = length(indivs),ncol = 1)
i <- 1
for (i in 1:length(indivs)) {
  tmp <- which(df$ID_indiv == indivs[i])
  if(length(tmp)>0){
    checks[i] <- tmp
  }
}

lw <- df[checks,]
save(best_config,file = "Software_outputs/test_best_config.rda")
save(lw, file = "Input/test_length_weight.rda")
