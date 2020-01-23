#genepop_create - function
#function to create a genepop type output for a SNP data set

#Inputs:
#df - dataframe with pop and ID for all individuals, with all SNP genotype calls in the columns
#output_file - character string that's the desired name of the output file
#title - the character string that goes on the top of the file
genepop_create <- function(df,output_file = "genepop_file.txt",title){
  #making line 1, which is a character string about the file
  cat(title,file = output_file,sep = "\n",append = T)
  #making SNP name section, which is a list of SNPs, one per line
  tmp <- unique(df$SNP)
  cat(tmp,file = output_file,sep = ", ",append = T)
  cat("\n",file = output_file,append = T)
  #adding actual genotypes by population
  tmp <- unique(df$pop)
  i <- 1
  for (i in 1:length(tmp)) {
    pop_tmp <- tmp[i]
    df1 <- df %>% 
      filter(pop == pop_tmp) %>% 
      select(-pop) %>% 
      mutate(name = paste0(indiv,",")) %>% 
      select(-indiv) %>% 
      select(name,everything()) %>% 
      spread(key = SNP,value = gt)
    cat("POP",file = output_file,sep = "\n",append = T)
    write.table(df1,file = output_file,append = T,quote = F,sep = " ",row.names = F,col.names = F)
  }
}
