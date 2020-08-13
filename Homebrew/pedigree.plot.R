#pedigree.plot

#inputs a COLONY file with a column for offspring, moms, and dads
#also includes a cohort column (created earlier)

pedigree.plot <- function(family,title = "Pedigree Plot"){
  family <- family[order(family$cohort),]
  moms <- unique(family$Mother)
  dads <- unique(family$Father)
  
  momdots <- data.frame(ID = c(1:length(moms)), y = seq(from=1,to=length(family$ID),length.out=length(moms)))
  daddots <- data.frame(ID = c(1:length(dads)), y = seq(from=1,to=length(family$ID),length.out=length(dads)))
  
  #make the plot points
  plot(x = c(1,2,3), 
       y = c(0, length(family$ID), length(family$ID)), 
       pch = "", axes = FALSE, 
       xlab = "", ylab = "",
       main = title)
  points(x = rep(2,length(family$ID)), y = c(1:length(family$ID)), pch = "-", cex = 0.25)
  points(x = rep(1,length(moms)), y = momdots$y, pch = 19, cex = 0.5)
  points(x = rep(3,length(dads)), y = daddots$y, pch = 19, cex = 0.5)
  
  #make the lines in the plot
  for(r in 1:length(family$ID)) {
    #mom line
    lines(x=c(2,1), y = c(r, momdots$y[family$Mother[r]]), lwd = 0.5,col = rgb(0,0,0,alpha = 0.5))
    
    #dad line
    lines(x =c(2,3), y = c(r, daddots$y[family$Father[r]]), lwd = 0.5,col = rgb(0,0,0,alpha = 0.5))
  }
  
  #adding labels
  mtext("Parent 1", side = 3, line = 0, at = 1, cex = 0.9)
  mtext("Parent 2", side = 3, line = 0, at = 3, cex = 0.9)
  mtext("Offspring", side = 3, line = 0, at = 2, cex = 0.9)
  #adding rectangles
  cohorts <- unique(family$cohort)
  for (r in 1:length(cohorts)) {
    c <- cohorts[r]
    rect(xleft=1.75,
         xright = 2.25,
         ybottom = max(which(family$cohort == c)), 
         ytop = min(which(family$cohort == c)),
         border = "black", lwd = 2)
    text(x = 1.7, y = mean(which(family$cohort == c)), c, pos =2, col = "black")
  }
}
