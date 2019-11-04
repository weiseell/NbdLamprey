colony.dat.create.R <- function(moms = NA, dads = NA, kids, markers, update.alfs = 0,
                                spp.type = 2, inbreeding = 0, ploidy = 0, fem.gamy = 1, mal.gamy = 1,
                                clone = 0, sib.scale = 1, sib.prior = 0, known.alfs = 0, run.number = 1,
                                run.length = 2, monitor = 0, windows.version = 0, full.likelihood = 1,
                                likelihood.precision = 2, prob.mom = 1.0, prob.dad = 1.0){
  
  #getting current working directory and fixing slashes for running on linux
  my.wd <- "'Input.data'"
  my.wd2 <- "'Output.data'" 
  
  #getting the number of kids, moms, and dads
  noff1 <- nrow(kids); nmoms <- nrow(moms); ndads <- nrow(dads)
  noff <- paste(noff1,"! Number of offspring in the sample",sep = "\t")
  
  #getting the number of loci
  nloci <- ncol(markers)
  nloci <- paste(nloci, "! Number of loci",sep = "\t") 
  
  #setting a random number seed
  random.seed <- round(runif(n = 1,min = 1,max = 9999))
  random.seed <- paste(random.seed, "! Seed for random number generator",sep = "\t")
  
  #adding comments to input parameters
  update.alfs <- paste(update.alfs,"! Not upate/update allele frequency",sep = "\t")
  spp.type <- paste(spp.type, "! 2/1=Dioecious/Monoecious",sep = "\t")
  inbreeding <- paste(inbreeding,"! 0/1=No inbreeding/inbreeding", sep = "\t")
  ploidy <- paste(ploidy, "! 0/1=Diploid species/HaploDiploid species",sep="\t")
  gamy <- paste(paste(mal.gamy,fem.gamy),"! 0/1=Polygamy/Monogamy for males & females",sep = "\t")
  clone <- paste(clone,"! 0/1=Clone inference =No/Yes",sep="\t")
  sib.scale <- paste(sib.scale,"! 0/1=Scale full sibship=No/Yes",sep = "\t")
  sib.prior <- paste(sib.prior,"! 0/1/2/3=No sibship prior/Weak sibship prior/Medium sibship prior/Strong sibship prior",sep = "\t")
  known.alfs <- paste(known.alfs,"! 0/1=Unknown/Known population allele frequency",sep = "\t")
  run.number <- paste(run.number,"! Number of runs",sep = "\t")
  run.length <- paste(run.length,"! Length of run",sep = "\t")
  monitor <- paste(monitor,"! 0/1=Monitor method by Iterate#/Time in second",sep = "\t")
  monitor.interval <- paste(1000000,"! Monitor interval in Iterate# / in seconds",sep = "\t")
  windows.version <- paste(windows.version,"! Windows version",sep = "\t")
  full.likelihood <- paste(full.likelihood,"! Fulllikelihood",sep = "\t")
  likelihood.precision <- paste(likelihood.precision,"! 1/2/3=low/medium/high Precision for Fulllikelihood",sep = "\t")
  
  #collating info for parents
  if(prob.dad == 0 & prob.mom == 0){
    probs <- c("0.0","0.0")
  } else {
    probs <- c(prob.dad,prob.mom)
  }

  if(is.na(dads) & is.na(moms)){
    npars <- c(0,0)
  } else {
    npars <- c(nrow(dads),nrow(moms))
  }

  my.value <- paste0(0,"\n")
  my.exc.value <- paste(0,0,"\n")
  #making the actual file
  cat(my.wd,my.wd2,noff,nloci,random.seed, update.alfs,spp.type,
      inbreeding,ploidy,gamy,clone,sib.scale,sib.prior,
      known.alfs,run.number,run.length,monitor,monitor.interval,
      windows.version,full.likelihood,likelihood.precision,file = "colony2.dat",sep = "\n",append = T)
  cat("\n",file = "colony2.dat",append = T)
  write.table(x = markers,file = "colony2.dat",append = T,quote = F,sep = ",",row.names = F,col.names = T)
  cat("\n",file = "colony2.dat",append = T)
  write.table(x = kids,file = "colony2.dat",append = T,quote = F,sep = " ",row.names = F,col.names = F)
  cat("\n",file = "colony2.dat",append = T)
  cat(probs,file = "colony2.dat",append = T)
  cat("\n",file = "colony2.dat",append = T)
  cat(npars,file = "colony2.dat",append = T)
  cat("\n",file = "colony2.dat",append = T)
  if(!is.na(dads) & !is.na(moms)){
    cat("\n",file = "colony2.dat",append = T)
    write.table(x = dads,file = "colony2.dat",append = T,quote = F,sep = " ",row.names = F,col.names = F)
    cat("\n",file = "colony2.dat",append = T)
    write.table(x = moms,file = "colony2.dat",append = T,quote = F,sep = " ",row.names = F,col.names = F)
  }
  cat("\n",file = "colony2.dat",append = T)
  cat(my.exc.value,file = "colony2.dat",append = T)
  cat("\n",file = "colony2.dat",append = T)
  cat(my.exc.value,file = "colony2.dat",append = T)
  cat("\n",file = "colony2.dat",append = T)
  cat(my.value,file = "colony2.dat",append = T)
  cat("\n",file = "colony2.dat",append = T)
  cat(my.value,file = "colony2.dat",append = T)
  cat("\n",file = "colony2.dat",append = T)
  cat(my.value,file = "colony2.dat",append = T)
  cat("\n",file = "colony2.dat",append = T)
  cat(my.value,file = "colony2.dat",append = T)
  cat("\n",file = "colony2.dat",append = T)
  cat(my.value,file = "colony2.dat",append = T)
  cat("\n",file = "colony2.dat",append = T)
  cat(my.value,file = "colony2.dat",append = T)
}

