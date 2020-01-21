#separating expansion locations into cohorts

#loading libraries
#Note: the requirements are built into the function so you don't need to load this before running,
#but it's good to have the requirements in case you don't have these packages in R
library(tidyverse)
library(mclust)

#loading in my functions
source("Homebrew/mixture_function.R")
source("Homebrew/multiplot.R")
#loading in data
df <- read.table("Input/exp_lengths_weights_081219.txt",header = T, stringsAsFactors = F)
che_colony <- read.table("Input/colony.bestconfig.che.txt",header = T,sep = "\t",stringsAsFactors = F)
bmr_colony <- read.table("Input/colony.bestconfig.bmr.txt",header = T,stringsAsFactors = F)
ocq_colony <- read.table("Input/colony.bestconfig.ocq.txt",header = T,stringsAsFactors = F)
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

#running the mixture function for all locations and collections
bmr17_mix <- mixture(bmr17, pop = "BMR_17")
bmr18_mix <- mixture(bmr18, pop = "BMR_18")
bmr19_mix <- mixture(bmr19, pop = "BMR_19")
che18_mix <- mixture(che18, pop = "CHE_18")
ocq18_mix <- mixture(ocq18, pop = "OCQ_18")

#make graphs and determine length ranges for each cohort for each group ####
#OCQ
ocq_plot1 <- ggplot(ocq18_mix[[2]], aes(x=V1, y=V2,color = class)) +
  geom_point(aes(color = factor(class))) +
  scale_color_manual(values = c("#2c7fb8","#41b6c4","#253494"),
                     name = "Age",
                     labels = c("3", "2","4"),
                     guide = F)+
  labs(x="Length (mm)", y="Weight (g)")+
  theme_bw(base_size = 8)+
  ggtitle("Age Classifications for Ocqueoc River Individuals")
ocq_plot2 <- ggplot(ocq18_mix[[2]], aes(x=V1, fill = class)) +
  geom_histogram(aes(fill = factor(class)),bins = 100) +
  scale_fill_manual(values = c("#2c7fb8","#41b6c4","#253494"),
                     name = "Age",
                     labels = c("3", "2","4"),
                     guide = guide_legend(reverse = FALSE))+
  labs(x="Length (mm)", y="counts")+
  theme_bw(base_size = 8)+
  ggtitle("Age Classifications for Ocqueoc River Individuals")

#BMR - 2019 collection
bmr19_plot1 <- ggplot(bmr19_mix[[2]], aes(x=V1, y=V2,color = class)) +
  geom_point(aes(color = factor(class))) +
  scale_color_manual(values = c("#2c7fb8"),
                     name = "Age",
                     labels = c("3"),
                     guide = F) +
  labs(x="Length (mm)", y="Weight (g)") +
  theme_bw(base_size = 8) +
  ggtitle("Age Classifications for Black Mallard River Individuals - \n2019 collection")
bmr19_plot2 <- ggplot(bmr19_mix[[2]], aes(x=V1, fill = class)) +
  geom_histogram(aes(fill = factor(class)),bins = 20) +
  scale_fill_manual(values = c("#2c7fb8"),
                     name = "Age",
                     labels = c("3"))+
  labs(x="Length (mm)", y="counts")+
  theme_bw(base_size = 8)+
  ggtitle("Age Classifications for Black Mallard River Individuals - \n2019 collection")

#BMR - 2018 collection
bmr18_plot1 <- ggplot(bmr18_mix[[2]], aes(x=V1, y=V2,color = class)) +
  geom_point(aes(color = factor(class))) +
  scale_color_manual(values = c("#2c7fb8","#41b6c4","#a1dab4"),
                     name = "Age",
                     labels = c("3", "2","1"),
                     guide = F)+
  labs(x="Length (mm)", y="Weight (g)")+
  theme_bw(base_size = 8)+
  ggtitle("Age Classifications for Black Mallard River Individuals - \n2018 collection")
bmr18_plot2 <- ggplot(bmr18_mix[[2]], aes(x=V1, fill = class)) +
  geom_histogram(aes(fill = factor(class)),bins = 100) +
  scale_fill_manual(values = c("#2c7fb8","#41b6c4","#a1dab4"),
                     name = "Age",
                     labels = c("3", "2","1"),
                     guide = guide_legend(reverse = TRUE))+
  labs(x="Length (mm)", y="counts")+
  theme_bw(base_size = 8)+
  ggtitle("Age Classifications for Black Mallard River Individuals - \n2018 collection")
#age 1 - 60mm, age 2 - 88mm, age 3 - 113mm

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
#age 0 - 32mm, age 1 - 75mm, age 2 - 100mm

che_plot1 <- ggplot(che18_mix[[2]], aes(x=V1, y=V2,color = class)) +
  geom_point(aes(color = factor(class))) +
  scale_color_manual(values = c("#7fcdbb","#41b6c4","#2c7fb8"),
                     name = "Age",
                     labels = c("1", "2","3"),
                     guide = F)+
  labs(x="Length (mm)", y="Weight (g)")+
  theme_bw(base_size = 8)+
  ggtitle("Age Classifications for Cheboygan River")
che_plot2 <- ggplot(che18_mix[[2]], aes(x=V1, fill = class)) +
  geom_histogram(aes(fill = factor(class)),bins = 100) +
  scale_fill_manual(values = c("#7fcdbb","#41b6c4","#2c7fb8"),
                     name = "Age",
                     labels = c("1", "2","3"),
                     guide = guide_legend(reverse = FALSE))+
  labs(x="Length (mm)", y="counts")+
  theme_bw(base_size = 8)+
  ggtitle("Age Classifications for Cheboygan River")
#age 1 - 69mm, age 2 - 119mm
#putting all 10 plots together
multiplot(che_plot1,bmr17_plot1,bmr18_plot1,bmr19_plot1,ocq_plot1,che_plot2,bmr17_plot2,bmr18_plot2,bmr19_plot2,ocq_plot2,cols = 2)

#making cutoffs and splitting into cohorts####
cutoffs <- data.frame(loc = as.character(c("BMR","BMR","BMR","BMR","BMR","BMR","CHE","CHE","CHE","OCQ","OCQ")),
           year = c(2017,2017,2017,2018,2018,2018,2018,2018,2018,2018,2018),
           age = c(0,1,2,1,2,3,1,2,3,2,3),
           min = as.integer(c(0,33,76,0,61,89,0,70,119,0,100)),
           max = as.integer(c(32,75,125,60,88,134,69,119,146,100,160)))
cutoffs$rep_year <- cutoffs$year - cutoffs$age

i <- 1
j <- 1
df1 <- df
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

bmr <- subset(df1,df1$loc == "BMR")
che <- subset(df1,df1$loc == "CHE")
ocq <- subset(df1,df1$loc == "OCQ")
che <- che %>% 
  mutate(ID = paste(species,loc,num,sep = "_")) %>% 
  select(-Sample_number:-num) %>% 
  rename(Length = V1,Weight = V2) %>% 
  select(ID,Year_collect:cohort)
colnames(che_colony) <-  c("ID", "Father","Mother","Cluster")
che1 <- merge(che,che_colony)

che1 %>% 
  arrange(desc(Length))

bmr <- bmr %>% 
  mutate(ID = paste(species,loc,num,sep = "_")) %>% 
  select(-Sample_number:-num) %>% 
  rename(Length = V1,Weight = V2) %>% 
  select(ID,Year_collect:cohort)

colnames(bmr_colony) <-  c("ID", "Father","Mother","Cluster")
bmr1 <-  merge(bmr,bmr_colony)
bmr2 <- bmr1 %>% 
  arrange(Year_collect,desc(Length),Father,Mother)

ocq <- ocq %>% 
  mutate(ID = paste(species,loc,num,sep = "_")) %>% 
  select(-Sample_number:-num) %>% 
  rename(Length = V1,Weight = V2) %>% 
  select(ID,Year_collect:cohort)
colnames(ocq_colony) <-  c("ID", "Father","Mother","Cluster")
ocq1 <- merge(ocq,ocq_colony)
ocq2 <- ocq1 %>% 
  arrange(desc(Length))
#separating out reaches
bmr2 <- bmr1 %>% 
  mutate(ID2 = ID) %>% 
  separate(col = ID2,into = c("spp","loc","num"),sep = "_") %>% 
  select(spp:num,ID:Cluster) %>% 
  arrange(num)

bmr2$reach <- NA
for (i in 1:length(bmr2$num)) {
  if(as.numeric(bmr2$num[i]) > 0 & as.numeric(bmr2$num[i]) <= 128){
    bmr2$reach[i] <- 1
  }
  if(as.numeric(bmr2$num[i]) > 128 & as.numeric(bmr2$num[i]) <= 387){
    bmr2$reach[i] <- 2
  }
  if(as.numeric(bmr2$num[i]) > 387 & as.numeric(bmr2$num[i]) <= 1054){
    bmr2$reach[i] <- 3
  }
  if(as.numeric(bmr2$num[i]) > 1054){
    bmr2$reach[i] <- 4
  }
  
}




