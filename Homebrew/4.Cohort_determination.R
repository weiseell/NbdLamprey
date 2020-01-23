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
load("Software_outputs/test_best_config.rda")
load("Input/test_length_weight.rda")

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
#which I got from the outputs manually because the clusters have overlap that I can't fix
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
colnames(best_config) <-  c("ID", "Father","Mother","Cluster")
df2 <-  merge(df1,best_config)
df2 <- df2 %>% 
  arrange(Year_collect,desc(Length),Father,Mother) %>% 
  mutate(family = paste(Mother,Father,sep = "_"))

#Goal 2####
#plots to look at how well the mixture analysis and genetic data overlap
#looking at length distribution (families w/ multiple offspring only)
#2017 group
df17 <- subset(df2,df2$Year_collect == 2017)
df17 <- df17[df17$family %in% names(which(table(df17$family)>=3)),]
#boxplot
boxplot(df17$Length~df17$family, las = 2,
        main = "Black Mallard length distribution by family \n2017 Collection",
        xlab = "family",
        ylab = "length (mm)",
        ylim = c(0,130),
        cex.axis = 0.6)
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




