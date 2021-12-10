#match_tags:
#function to identify target loci from the rapture panel

#Inputs
#SNPs: dataframe of genotyped SNP locations on a reference genome
  #should have two columns: CHROM and POS
#tags: dataframe of tag locations on the same reference genome as genotyped SNPs
  #should have four columns: ID, CHROM, min, and max
#target:logical that determines whether to output all SNPs with locations, or filter out off-target loci
  #!# off-target loci are marked as "NonTarget" in both outputs

#Outputs
#list of what targets each SNP is mapped to or
#list of only targeted SNPS

match_tags <- function(SNPs,tags,target = F){
  i <- 1
  SNPs$target <- NA
  for (i in 1:length(SNPs$CHROM)) {
    tmp <- which(tags$CHROM == SNPs$CHROM[i] & as.numeric(tags$min) <= SNPs$POS[i] & as.numeric(tags$max) >= SNPs$POS[i])
    ifelse(length(tmp) > 0,SNPs$target[i] <- tags$ID[tmp], SNPs$target[i] <- "NonTarget")
  }
  target <- subset(SNPs,SNPs$target != "NonTarget")
  list(SNPs,target)
}
