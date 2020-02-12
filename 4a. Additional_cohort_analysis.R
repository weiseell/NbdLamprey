#Supplmental analysis for the SupCon meeting
#may add some of these figures to the publication?
#Created on: Feb 10th, 2020

#Cohort determination
#Main goals:
#1.Perform Gaussian mixture analysis on length and weight data
#2.Comparing the results of the mixture analysis with genetic data

#libraries
library(tidyverse)
library(mclust)
#homebrew functions
source("Homebrew/mixture_function.R")
source("Homebrew/multiplot.R")

#load in data
bmr_colony <- read.table("Input_fulldata/colony.bestconfig.bmr.txt",header = T,sep = "\t",stringsAsFactors = F)
lw <- read.table("Input_fulldata/exp_lengths_weights.txt",header = T,sep = "\t",stringsAsFactors = F)

#Goal 1####
#separating data by location and age collected
lw <- lw %>% separate(ID_indiv, into = c("species","loc","num"),sep = "_") %>% 
  rename(V1 = Length, V2 = Weight)
bmr17 <- subset(lw, lw$loc == "BMR" & lw$Year_collect == "2017")
bmr18 <- subset(lw, lw$loc == "BMR" & lw$Year_collect == "2018")
bmr19 <- subset(lw, lw$loc == "BMR" & lw$Year_collect == "2019")
#running mixture analysis
bmr17_mix <- mixture(bmr17, pop = "BMR_17")
bmr18_mix <- mixture(bmr18, pop = "BMR_18")
bmr19_mix <- mixture(bmr19, pop = "BMR_19")

#use outputs to get cutoffs for each year
cutoff_all <- rbind(bmr17_mix[[3]],bmr18_mix[[3]],bmr19_mix[[3]])
#!# need a coding way to do this, right now I'm going to use the cutoffs from the actual project
#which I got from the outputs manually because the clusters have overlap that I can't fix ://///
#hopefully developing the new analysis will fix this problem anyways
cutoffs <- data.frame(loc = "BMR",
                      year = c(2017,2017,2017,2018,2018,2018),
                      age = c(0,1,2,1,2,3),
                      min = as.integer(c(0,33,76,0,61,89)),
                      max = as.integer(c(32,75,125,60,88,134)))
cutoffs$rep_year <- cutoffs$year - cutoffs$age

#adding cohorts to the mixture analysis
i <- 1
j <- 1
df1 <- lw
df1$cohort <- NA
for (i in 1:length(df1$species)) {
  tmp <- df1[i,]
  for (j in 1:length(cutoffs$loc)) {
    if(tmp$loc == cutoffs$loc[j] & tmp$Year_collect == cutoffs$year[j] & tmp$V1 >= cutoffs$min[j] & tmp$V1 <= cutoffs$max[j]){
      df1$cohort[i] <- cutoffs$rep_year[j]
    }
  }
  if(tmp$loc == "BMR" & tmp$Year_collect == 2019){
    df1$cohort[i] <- 2016
  }
}

#merging mixture analysis and genetic data
df1 <- df1 %>% 
  mutate(ID = paste(species,loc,num,sep = "_")) %>% 
  select(-Sample_number:-num) %>% 
  rename(Length = V1,Weight = V2) %>% 
  select(ID,Year_collect:cohort)
colnames(bmr_colony) <-  c("ID", "Father","Mother","Cluster")
df2 <-  merge(df1,bmr_colony)
df2 <- df2 %>% 
  arrange(Year_collect,desc(Length),Father,Mother) %>% 
  mutate(family = paste(Mother,Father,sep = "_"))

#creating half-sibling family clusters
df2$halfsib <- NA
moms <- unique(df2$Mother)

#find a way to eliminate moms that are already in groups from the mom vector
#have a while loop that escapes when the moms are empty
#keep track of the # of moms and # of dads in each group
i <- 1
counts <- data.frame(GroupID = character(0),
                     nmoms = character(0),
                     ndads = character(0),
                     noff = character(0),stringsAsFactors = F)
while (length(moms)>0) {
  #define a half-sib group
  group <- i
  #getting a unique group of mothers
  mom <- moms[1]
  tmp <- which(df2$Mother == mom)
  #getting the corresponding group of dads
  tmp2 <- unique(df2$Father[tmp])
  tmp3 <- which(df2$Father%in%tmp2)
  #getting the moms from this group (this should fix the NA problem)
  tmp4 <- unique(df2$Mother[tmp3])
  tmp5 <- which(df2$Mother%in%tmp4)
  #defining the group of offspring in the group
  df2[unique(c(tmp5,tmp3)),]$halfsib <- group
  #finding which moms are in the full group and eliminating them from the mom vector
  allmom <- unique(df2[unique(c(tmp,tmp3)),]$Mother)
  moms <- moms[-which(moms%in%allmom)]
  
  #counting #of moms/dads, #of offspring per group
  counts[nrow(counts)+1,] <- list(GroupID = i,
                                  nmoms = length(allmom),
                                  ndads = length(unique(df2[unique(c(tmp,tmp3)),]$Father)),
                                  noff = length(df2[unique(c(tmp,tmp3)),]$ID))
  #group counter
  i <- i+1
}

#calculating the average number of moms/dads, 
counts$ndads <- as.integer(counts$ndads)
counts$nmoms <- as.integer(counts$nmoms)
counts$noff <- as.integer(counts$noff)

counts <- counts %>% 
  mutate(npar = ndads+nmoms)

counts %>% 
  summarise(mean_moms = mean(nmoms),
            mean_dads = mean(ndads),
            mean_off = mean(noff),
            mean_par = mean(npar))
#looking for isolated families (npar = 2)
table(counts$npar)

#and the average number of mates per parent for each half-sibling group
momcount <- matrix(data = NA,ncol = 1,nrow = length(unique(df2$Mother)))
for (i in 1:length(unique(df2$Mother))){
  tmp <- which(df2$Mother == i)
  momcount[i] <- length(unique(df2$Father[tmp]))
}
dadcount <- matrix(data = NA,ncol = 1,nrow = length(unique(df2$Father)))
for (i in 1:length(unique(df2$Father))){
  tmp <- which(df2$Father == i)
  dadcount[i] <- length(unique(df2$Mother[tmp]))
}
mean(momcount)
mean(dadcount)
mean(momcount+dadcount)

#making a counts file for the cluster
counts2 <- data.frame(GroupID = character(0),
                      nmoms = character(0),
                      ndads = character(0),
                      noff = character(0),stringsAsFactors = F)
Clusters <- unique(df2$Cluster)
for(i in 1:length(Clusters)){
  tmp <- which(df2$Cluster == Clusters[i])
  momcount <- length(unique(df2$Mother[tmp]))
  dadcount <- length(unique(df2$Father[tmp]))
  counts2[nrow(counts2)+1,] <- list(GroupID = Clusters[i],
                                  nmoms = momcount,
                                  ndads = dadcount,
                                  noff = length(tmp))
}

#doing the same for the cluster analysis
Clusters <- unique(df2$Cluster)
parents <- list(character(0))
for(i in 1:length(Clusters)){
  tmp <- df2[which(df2$Cluster == Clusters[i]),]
  mom <- unique(tmp$Mother)
  dad <- unique(tmp$Father)
  name <- paste0("Cluster",Clusters[i])
  parents[[name]] <- matrix(data = c("Moms",mom,"Dads",dad))
}
#Goal 2####
#plots to look at how well the mixture analysis and genetic data overlap
#looking at length distribution (families w/ multiple offspring only)
#2017 group
df17 <- subset(df2,df2$Year_collect == 2017)
df17 <- df17[df17$family %in% names(which(table(df17$family)>=3)),]
#boxplot for full-sibling groups
boxplot(df17$Length~df17$family, las = 2,
        main = "Black Mallard length distribution by full-sibling family \n2017 Collection",
        xlab = "family",
        ylab = "length (mm)",
        ylim = c(0,130),
        cex.axis = 0.5)
#boxplot for half-sibling groups
boxplot(df17$Length~df17$halfsib,
        main = "Black Mallard length distribution by half-sibling group \n2017 Collection",
        xlab = "half-sibling group",
        ylab = "length (mm)",
        ylim = c(0,130),
        cex.axis = 0.7)
#boxplot for cluster
boxplot(df17$Length~df17$Cluster,
        main = "Black Mallard length distribution by half-sibling group \n2017 Collection",
        xlab = "half-sibling group",
        ylab = "length (mm)",
        ylim = c(0,130),
        cex.axis = 0.9)
#pie chart
df3 <- df2[df2$family %in% names(which(table(df2$family)>=3)),]
cp <- coord_polar(theta = "y")
cp$is_free <- function() TRUE
counts1 <- df3 %>% 
  group_by(family,cohort) %>% 
  tally()
ggplot(counts1,aes(x = family,y = n,fill = factor(cohort)))+
  geom_bar(stat="identity")+
  cp+
  facet_wrap(~family,scales = "free")+
  theme_bw()+ 
  scale_x_discrete(NULL, expand = c(0,0)) +
  scale_y_continuous(NULL, expand = c(0,0)) +
  scale_fill_grey(name = "Cohort")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x=element_blank(),
        axis.text.y = element_blank(),
        aspect.ratio = 1) +
  ggtitle("Cohort identification for Black Mallard full-sibling families")

#2018 group
#2017 group
df18 <- subset(df2,df2$Year_collect == 2018)
df18 <- df18[df18$family %in% names(which(table(df18$family)>=3)),]
#boxplot for full-sibling groups
boxplot(df18$Length~df18$family, las = 2,
        main = "Black Mallard length distribution by full-sibling family \n2018 Collection",
        xlab = "family",
        ylab = "length (mm)",
        ylim = c(0,130),
        cex.axis = 0.5)
#boxplot for half-sibling groups
boxplot(df18$Length~df18$halfsib,
        main = "Black Mallard length distribution by half-sibling group \n2018 Collection",
        xlab = "half-sibling group",
        ylab = "length (mm)",
        ylim = c(0,130),
        cex.axis = 0.7)
#boxplot for colony cluster
boxplot(df18$Length~df18$Cluster,
        main = "Black Mallard length distribution by half-sibling group \n2018 Collection",
        xlab = "half-sibling group",
        ylab = "length (mm)",
        ylim = c(0,130),
        cex.axis = 0.9)
#length-frequency histogram with half-sibling groups color coded
#2017 group
ggplot(df17,aes(x = Length,fill = Cluster))+
  geom_histogram(aes(fill = factor(halfsib)),bins = 50)+
  labs(x = "Length (mm)",fill = "Half-sibling \n group")+
  ggtitle("Length histogram - Black Mallard 2017 cohort")

#2018 group
ggplot(df18,aes(x = Length,fill = Cluster))+
  geom_histogram(aes(fill = factor(halfsib)),bins = 50)+
  labs(x = "Length (mm)",fill = "Half-sibling \n group")+
  ggtitle("Length histogram - Black Mallard 2018 cohort")+
  geom_vline(xintercept = 60)+
  geom_vline(xintercept = 88)+
  geom_vline(xintercept = 113)

