#separating expansion locations into cohorts

#loading libraries
#Note: the requirements are built into the function so you don't need to load this before running,
#but it's good to have the requirements in case you don't have these packages in R
library(tidyverse)
library(mclust)

#loading in my functions
source("Homebrew/mixture_function.R")
#loading in data
df <- read.table("Input/exp_lengths_weights_081219.txt",header = T, stringsAsFactors = F)
#subsetting for locations and year collected
#separating ID into species, location, and individual number
df <- df %>% separate(ID_indiv, into = c("species","loc","num"),sep = "_") %>% 
  rename(V1 = Length, V2 = Weight)
#look at the table of locations and collection year to figure out which locations
#and collection years that need to be analyzed
bmr17 <- subset(df, df$loc == "BMR" & df$Year_collect == "2017")

#running the mixture function
bmr17_mix <- mixture(bmr17, pop = "BMR_17")

#Plots
bmr17_plot1 <- ggplot(bmr17_mix[[2]], aes(x=V1, y=V2,color = class)) +
  geom_point(aes(color = factor(class))) +
  scale_color_manual(values = c("#c7e9b4","#7fcdbb","#41b6c4","#2c7fb8"),
                     name = "Age",
                     labels = c("0", "1","2","3"),
                     guide = F)+
  labs(x="Length (mm)", y="Weight (g)")+
  theme_bw(base_size = 8) +
  ggtitle("Age Classifications for Black Mallard River Individuals - \n2017 collection")
bmr17_plot2 <- ggplot(bmr17_mix[[2]], aes(x=V1, fill = class)) +
  geom_histogram(aes(fill = factor(class)),bins = 100) +
  scale_fill_manual(values = c("#c7e9b4","#7fcdbb","#41b6c4","#2c7fb8"),
                    name = "Age",
                    labels = c("0", "1","2","3"),
                    guide = guide_legend(reverse = TRUE))+
  labs(x="Length (mm)")+
  theme_bw(base_size = 8)+
  ggtitle("Age Classifications for Black Mallard River Individuals - \n2017 collection")