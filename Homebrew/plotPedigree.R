df <- che1 %>% select(ID,Mother,Father,cohort)
#colnames(df) <- c("offspringID", "momID", "dadID", "cohort")
df1 <- ocq1 %>% select(ID,Mother,Father,cohort)
#colnames(df1) <- c("offspringID", "momID", "dadID", "cohort")
df2 <- bmr2 %>% select(ID,Mother,Father,cohort)
#colnames(df2) <- c("offspringID", "momID", "dadID", "cohort")
par(mfrow=c(1,3))
text(0.5,0.5,"Visualization of Pedigree Reconstruction",cex=2,font=2)
plotPedigree(df2,title = "Black Mallard River")
plotPedigree(df,title = "Cheboygan River")
plotPedigree(df1,title = "Ocqueoc River")
plotPedigree <- function(df,title){
  df <- df %>% arrange(cohort,Mother,Father)
  cohort <- unique(df$cohort)
  offspringID <- seq(from=1,to=length(df$ID),length.out=length(unique(df$ID)))
  momlocs <- data.frame(ID = unique(df$Mother), x = rep(1,length(unique(df$Mother))), y = seq(from=1,to=length(offspringID),length.out=length(unique(df$Mother))))
  dadlocs <- data.frame(ID = unique(df$Father), x = rep(3,length(unique(df$Father))), y = seq(from=1,to=length(offspringID),length.out=length(unique(df$Father))))
  
  #now make the plot
  plot(x = c(1,2,3), 
       y = c(0, length(offspringID), length(offspringID)), 
       pch = "", axes = FALSE, xlab = "", ylab = "",main = title)
  points(x = rep(2,length(offspringID)), y = offspringID, pch = "-", cex = 0.4)
  points(x = momlocs$x, y = momlocs$y, pch = 19, cex = 0.5)
  points(x = dadlocs$x, y = dadlocs$y, pch = 19, cex = 0.5)
  r <- 1
  for(r in 1:length(df$ID)) {
    #mom line
    lines(x=c(2,1), y = c(offspringID[r], momlocs$y[which(momlocs$ID == df$Mother[r])]), lwd = 0.5, col = "grey")
    
    #dad line
    lines(x =c(2,3), y = c(offspringID[r], dadlocs$y[which(dadlocs$ID == df$Father[r])]), lwd = 0.5, col = "grey")
  }
  
  mtext("Inferred \nParent 1", side = 3, line = 0, at = 1, cex = 0.7)
  mtext("Inferred \nParent 2", side = 3, line = 0, at = 3, cex = 0.7)
  mtext("Offspring", side = 3, line = 0, at = 2, cex = 0.7)
  
  r <- 1
  for (r in 1:length(cohort)) {
    tmp <- cohort[r]
    rect(xleft=1.75,
         xright = 2.25,
         ybottom = max(which(df$cohort == tmp)), 
         ytop = min(which(df$cohort == tmp)),
         border = "black", lwd = 2)
    text(x = 1.7, y = mean(which(df$cohort == tmp)), paste0("Inferred \n cohort ",tmp), pos =2, col = "black",cex = 1.1)
  }
}
