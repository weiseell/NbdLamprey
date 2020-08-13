##Mixture analysis functions:
#Written by: Ellie Weise
#Last edited: 8/1/19

#Inputs:
#data table with at least three columns: ID, V1 and V2 (need these labels)
#Note: data table should already be sorted by location and year collected

#Output: a list with 3 elements
#1. complete mclust output
#2. Individual summary table with length, weight, class, and uncertainty for each individual
#3. Location summary table with Length/Weight ranges for each cluster, size of each class, model type,
#BIC value for best fit model, and means of each cluster for length and weight

mixture <- function(df, maxclust, pop = "pop1"){
  #running mclust
  #loading required packages
  require(tidyverse,mclust)
  #grabbing V1 and V2 from df and turning it into a matrix
  df1 <- as.matrix(data.frame(V1 = df$V1, V2 = df$V2))
  #mclust.BIC function:
  #df.BIC <- mclustBIC(df1, G = 1:maxclust, modelNames = c("VVV"))
  #running mclust
  df.mclust <- Mclust(df1, x = maxclust, modelNames = c("VVV"))
  
  ##interpreting mclust output
  #adding individual classification and uncertainty to df
  df <- df %>% mutate(class = df.mclust$classification, uncert = df.mclust$uncertainty)
  #subsetting data to only include data with low uncertainty
  nclass <- length(unique(df$class))
  
  #if data has more than one cluster, subset for low uncertainty
  #if data has one cluster, then take all data
  ifelse(nclass > 1, 
         df2 <- subset(df,df$uncert < ((1-1/nclass)/2)), 
         df2 <- df)
  
  #creating location summary table using tidyverse
  df2 <- df2 %>% 
    #group by indicates that all summary stats should be calculated for each class
    group_by(class) %>% 
    #summarize creates a table of custom summary statistics
    summarize(pop = pop,
              #length/weight ranges - calculated to include 95% of observed data
              low_Length = quantile(V1,probs = 0.025),
              high_Length = quantile(V1,probs = 0.975),
              low_Weight = quantile(V2,probs = 0.025),
              high_Weight = quantile(V2,probs = 0.975),
              #sample size
              n = length(V1),
              #model type
              model = df.mclust$modelName,
              #BIC of best fit model
              BIC = df.mclust$bic) %>% 
    #resorts the data table so each statistic is a row
    gather(key = stats,value = value,-class,-pop) %>% 
    #selects location, class, stat name and value columns and puts them in that order
    select(loc = pop, class, stats, value)
  
  #creates an object out of cluster means
  means <- data.frame(t(df.mclust$parameters$mean))
  colnames(means) <- c("Cluster_Mean_Length","Cluster_Mean_Weight")
  
  #adding cluster means to the summary table
  summ <- rbind(df2, means %>% 
                  mutate(class = seq(nrow(means))) %>% 
                  gather(key = stats, value = value,-class) %>% 
                  mutate(loc = pop) %>% 
                  select(loc,class,stats,value))
  #reordering data table so each statistic is in its own column
  #separate rows for each class
  summ <- summ %>% spread(key = stats,value = value)
  
  #output (described above)
  list(df.mclust,df,summ)
}





