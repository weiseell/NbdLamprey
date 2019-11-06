#separating expansion locations into cohorts

#setting working directory
setwd("/Users/ellenweise/OneDrive - Michigan State University/Documents/Sea_Lamprey_MS_project/mixture_analysis")

#loading libraries
#Note: the requirements are built into the function so you don't need to load this before running,
#but it's good to have the requirements in case you don't have these packages in R
library(tidyverse)
library(mclust)

#loading in my functions
source("mixture_function.R")
#loading in data
df <- read.table("input/exp_lengths_weights_081219.txt",header = T, stringsAsFactors = F)

#subsetting for locations and year collected
#separating ID into species, location, and individual number
df <- df %>% separate(ID_indiv, into = c("species","loc","num"),sep = "_") %>% 
  rename(V1 = Length, V2 = Weight)
#look at the table of locations and collection year to figure out which locations
#and collection years that need to be analyzed
bmr17 <- subset(df, df$loc == "BMR" & df$Year_collect == "2017")
bmr18 <- subset(df, df$loc == "BMR" & df$Year_collect == "2018")
bmr19 <- subset(df, df$loc == "BMR" & df$Year_collect == "2019")
che18 <- subset(df, df$loc == "CHE" & df$Year_collect == "2018")
ocq18 <- subset(df, df$loc == "OCQ" & df$Year_collect == "2018")
ocq19 <- subset(df, df$loc == "OCQ" & df$Year_collect == "2019")

#running the mixture function for all locations and collections
bmr17_mix <- mixture(bmr17, pop = "BMR_17")
bmr18_mix <- mixture(bmr18, pop = "BMR_18")
bmr19_mix <- mixture(bmr19, pop = "BMR_19")
che18_mix <- mixture(che18, pop = "CHE_18")
ocq18_mix <- mixture(ocq18, pop = "OCQ_18")
ocq19_mix <- mixture(ocq19, pop = "OCQ_19")

#make graphs and determine length ranges for each cohort for each group
ggplot(ocq18_mix[[2]], aes(x=V1, y=V2,color = class)) +
  geom_point(aes(color = factor(class))) +
  #scale_color_manual(values = c("#2ca25f","#2c7fb8","#810f7c"),
   #                  name = "Age",
    #                 labels = c("3", "2","1"),
     #                guide = guide_legend(reverse = TRUE))+
  labs(x="Length (mm)", y="Weight (g)")+
  theme_bw(base_size = 12)+
  ggtitle("Age Classifications for Ocqueoc River Individuals - 2018 cohort")
#age 1 cutoff - 40mm, age 2 cutoff - 60mm

ggplot(ocq19_mix[[2]], aes(x=V1, y=V2,color = class)) +
  geom_point(aes(color = factor(class))) +
  scale_color_manual(values = c("#2ca25f","#2c7fb8","#810f7c"),
                     name = "Age",
                     labels = c("3", "2","1"),
                     guide = guide_legend(reverse = TRUE))+
  labs(x="Length (mm)", y="Weight (g)")+
  theme_bw(base_size = 12)+
  ggtitle("Age Classifications for Ocqueoc River Individuals - 2019 cohort")+
  scale_x_continuous(limits = c(10,80))+
  scale_y_continuous(limits = c(0,1))

ggplot(bmr19_mix[[2]], aes(x=V1, y=V2,color = class)) +
  geom_point(aes(color = factor(class))) +
  scale_color_manual(values = c("#2ca25f","#2c7fb8","#810f7c"),
                     name = "Age",
                     labels = c("3", "2","1"),
                     guide = guide_legend(reverse = TRUE))+
  labs(x="Length (mm)", y="Weight (g)")+
  theme_bw(base_size = 12)+
  ggtitle("Age Classifications for Black Mallard River Individuals - 2019 cohort")
#only one cohort  -  age 3+?

ggplot(bmr18_mix[[2]], aes(x=V1, fill = class)) +
  geom_histogram(aes(fill = factor(class)),bins = 100) +
  labs(x="Length (mm)", y="Weight (g)")+
  theme_bw(base_size = 12)+
  ggtitle("Age Classifications for Black Mallard River Individuals - 2018 cohort")
#age 1 - 60, age 2 - 88mm, age 3 - 113mm

ggplot(bmr17_mix[[2]], aes(x=V1, fill = class)) +
  geom_histogram(aes(fill = factor(class)),bins = 100) +
  labs(x="Length (mm)")+
  theme_bw(base_size = 12)+
  ggtitle("Age Classifications for Black Mallard River Individuals - 2017 cohort")

#age 0 - 32mm, age 1 - 75mm, age 2 - 100mm

ggplot(che18_mix[[2]], aes(x=V1, y=V2,color = class)) +
  geom_point(aes(color = factor(class))) +
  scale_color_manual(values = c("#2ca25f","#2c7fb8","#810f7c"),
                     name = "Age",
                     labels = c("3", "2","1"),
                     guide = guide_legend(reverse = TRUE))+
  labs(x="Length (mm)", y="Weight (g)")+
  theme_bw(base_size = 12)+
  ggtitle("Age Classifications for Cheboygan River Individuals - 2018 cohort")
#age 1 - 69mm, age 2 - 119mm













