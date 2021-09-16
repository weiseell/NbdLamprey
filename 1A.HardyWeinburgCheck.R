pops <- c("BAD","BEI","BET","BRL","CAT","CHE","EAG","FOR","MAI","MAN","MIR","MIS","MUS","OCQ","STE","SWN","TAQ","TWO")
i <- 1
##HW check
popdev <- list()
popdevIDs <- list()
for (i in 1:length(pops)) {
  print(i)
  popsumstats <- read.table(file=paste0("Input/",pops[i],".populations.sumstats.tsv"),header = F,sep = "\t")
  colnames(popsumstats) <- c("locus","CHROM","POS","col","pop","Pnt","Qnt","N","P","HetO","HomO","HetE","HomE","Pi","SmoothPi","PiPval","Fis","SmoothFis","FisPval","HWEPval","Private")
  popsumstats$ID <- paste(popsumstats$CHROM,popsumstats$POS,sep = "_")
  
  #bonferroni correction for HWE deviance p-values
  popsumstats$bonferr <- p.adjust(popsumstats$HWEPval,method = "holm")
  popdev[[pops[i]]] <- popsumstats[which(popsumstats$bonferr == 0),]
  popdevIDs[[pops[i]]] <- popdev[[pops[i]]]$ID
}

dev_allpops <- Reduce(intersect,popdevIDs,accumulate = T)









