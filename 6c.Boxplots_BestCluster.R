#6c.Length Boxplots - separated by COLONY cluster 
#and highlighted by cluster likelihood

#load libraries
library(tidyverse)

#load homebrew functions
source("Homebrew/multiplot.R")
#load best cluster data
bmr_clust <- read.table(file = "Software_outputs/BMR_Output.data.BestCluster",header = T, sep = "\t",stringsAsFactors = F)
ocq_clust <- read.table(file = "Software_outputs/OCQ_Output.data.BestCluster",header = T, sep = "\t",stringsAsFactors = F)
load("Aging_Models/Family_data_all_locations.rda")
all_families <- all_families %>% select(-FatherID:-ClusterIndex)
clust_all <- rbind(ocq_clust,bmr_clust)
clust_all1 <- all_families %>% full_join(clust_all,by = c("OffspringID"))

#Black Mallard
bmr1 <- clust_all1 %>% 
  filter(samp == "BMR_2017" | samp == "BMR_2018" | samp == "BMR_2019")

give.n <- function(x){
  return(c(y = max(x)*1.08, label = length(x)))
  # experiment with the multiplier to find the perfect position
}

collect.labs <- c("2017 Collection","2018 Collection","2019 Collection")
names(collect.labs) <- c("BMR_2017", "BMR_2018", "BMR_2019")
bmr_plot <- ggplot(bmr1,aes(x=ClusterIndex,group=ClusterIndex,y=Length,fill=Probability))+
  facet_wrap(vars(samp),labeller =labeller(samp=collect.labs))+
  stat_summary(fun.data = give.n, geom = "text", fun = median,
               position = position_dodge(width = 0.75),size = 2)+
  geom_boxplot(alpha=0.3)+
  theme_bw()+
  scale_fill_gradient(low="red", high="white",name = "Cluster \nProbability")+
  xlab("Family Cluster")+
  ylab("Length (mm)")+
  ggtitle("A) Black Mallard River")

#Ocqueoc
ocq1 <- clust_all1 %>% 
  filter(loc == "OCQ")

collect.labs <- c("2018 Collection")
names(collect.labs) <- c("OCQ_2018")
ocq_plot <- ggplot(ocq1,aes(x=ClusterIndex,group=ClusterIndex,y=Length,fill=Probability))+
  facet_wrap(vars(samp),labeller =labeller(samp=collect.labs))+
  stat_summary(fun.data = give.n, geom = "text", fun = median,
               position = position_dodge(width = 0.75),size = 2)+
  geom_boxplot(alpha=0.3)+
  theme_bw()+
  scale_fill_gradient(low="red", high="white",name = "Cluster \nProbability")+
  xlab("Family Cluster")+
  ylab("Length (mm)")+
  ggtitle("B) Ocqueoc River")

#putting plots together
tiff(filename = "Figures/Length_Boxplots_prob.tiff",height = 8,width = 8,units = "in",res = 200)
multiplot(cols = 1,bmr_plot,ocq_plot)
dev.off()
