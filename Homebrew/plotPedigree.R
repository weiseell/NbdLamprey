df <- che1 %>% select(ID,Mother,Father,cohort)
colnames(df) <- c("offspringID", "momID", "dadID", "cohort")
df1 <- ocq1 %>% select(ID,Mother,Father,cohort)
colnames(df1) <- c("offspringID", "momID", "dadID", "cohort")
df2 <- bmr1 %>% select(ID,Mother,Father,cohort)
colnames(df2) <- c("offspringID", "momID", "dadID", "cohort")
par(mfrow=c(1,3))
text(0.5,0.5,"Visualization of Pedigree Reconstruction",cex=2,font=2)
plotPedigree(df2,title = "Black Mallard River")
plotPedigree(df,title = "Cheboygan River")
plotPedigree(df1,title = "Ocqueoc River")
plotPedigree <- function(df,title){
  df <- df[order(df$cohort),]
  cohort <- unique(df$cohort)
  offspringID <- seq(from=1,to=length(df$offspringID),length.out=length(unique(df$offspringID)))
  momlocs <- data.frame(ID = unique(df$momID), x = rep(1,length(unique(df$momID))), y = seq(from=1,to=length(offspringID),length.out=length(unique(df$momID))))
  dadlocs <- data.frame(ID = unique(df$dadID), x = rep(3,length(unique(df$dadID))), y = seq(from=1,to=length(offspringID),length.out=length(unique(df$dadID))))
  
  #now make the plot
  plot(x = c(1,2,3), y = c(0, length(offspringID), length(offspringID)), pch = "", axes = FALSE, xlab = "", ylab = "",main = title)
  points(x = rep(2,length(offspringID)), y = offspringID, pch = "-", cex = 0.25)
  points(x = momlocs$x, y = momlocs$y, pch = 19, cex = 0.5)
  points(x = dadlocs$x, y = dadlocs$y, pch = 19, cex = 0.5)
  r <- 1
  for(r in 1:length(df$offspringID)) {
    #mom line
    lines(x=c(2,1), y = c(offspringID[r], momlocs$y[which(momlocs$ID == df$momID[r])]), lwd = 0.5, col = "grey")
    
    #dad line
    lines(x =c(2,3), y = c(offspringID[r], dadlocs$y[which(dadlocs$ID == df$dadID[r])]), lwd = 0.5, col = "grey")
  }
  
  mtext("Parent 1", side = 3, line = 0, at = 1, cex = 0.9)
  mtext("Parent 2", side = 3, line = 0, at = 3, cex = 0.9)
  mtext("Offspring", side = 3, line = 0, at = 2, cex = 0.9)
  
  for (r in 1:length(cohort)) {
    tmp <- cohort[r]
    rect(xleft=1.75,
         xright = 2.25,
         ybottom = max(which(df$cohort == tmp)), 
         ytop = min(which(df$cohort == tmp)),
         border = "black", lwd = 2)
    text(x = 1.7, y = mean(which(df$cohort == tmp)), tmp, pos =2, col = "black")
  }
}
