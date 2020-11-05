#6c.Length Boxplots - separated by COLONY cluster 
#and highlighted by cluster likelihood

#load best cluster data
bmr_clust <- read.table(file = "Software_outputs/bmr_BestCluster_091820.txt",header = T, sep = "\t",stringsAsFactors = F)
ocq_clust <- read.table(file = "Software_outputs/ocq_BestCluster_091820.txt",header = T, sep = "\t",stringsAsFactors = F)
load("Aging_Models/Family_data_all_locations.rda")

clust_all <- rbind(ocq_clust,bmr_clust)
clust_all1 <- merge(clust_all,all_families2)

#Black Mallard
bmr1 <- clust_all1 %>% 
  filter(loc != "bmrAL") %>% 
  filter(loc != "OCQ")

collect.labs <- c("2017 Collection","2018 Collection")
names(collect.labs) <- c("BMR_2017", "BMR_2018")
bmr_plot <- ggplot(bmr1,aes(x=ClusterIndex,group=ClusterIndex,y=Length,fill=Probability))+
  facet_wrap(vars(samp),labeller =labeller(collect=collect.labs))+
  geom_boxplot(alpha=0.3)+
  theme_bw()+
  scale_fill_gradient(low="red", high="white",name = "Cluster \nLikelihood")+
  xlab("Family Cluster")+
  ylab("Length (mm)")+
  ggtitle("A) Downstream Black Mallard River")

#Ocqueoc
ocq1 <- clust_all1 %>% 
  filter(loc == "OCQ")

collect.labs <- c("2018 Collection")
names(collect.labs) <- c("OCQ_2018")
ocq_plot <- ggplot(ocq1,aes(x=ClusterIndex,group=ClusterIndex,y=Length,fill=Probability))+
  facet_wrap(vars(samp),labeller =labeller(collect=collect.labs))+
  geom_boxplot(alpha=0.3)+
  theme_bw()+
  scale_fill_gradient(low="red", high="white",name = "Cluster \n Likelihood")+
  xlab("Cluster")+
  ylab("Length (mm)")+
  ggtitle("B) Ocqueoc River")

#putting plots together
tiff(filename = "Figures/Length_Boxplots_prob.tiff",height = 10,width = 10,units = "in",res = 200)
multiplot(cols = 1,bmr_plot,ocq_plot)
dev.off()
