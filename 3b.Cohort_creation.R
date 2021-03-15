#NbdLamprey - Script 3b: Cohort creation
#Objective: Comparing length-based aging models and reconstructed pedigrees

#libraries
library(tidyverse)

#functions
source("Homebrew/pedigree.plot.R")
##reading in data
#load in data
all_locs <- read.table("AgingModels/lw_Bayes_assignments.txt",header = T,sep = "\t",stringsAsFactors = F)
all_locs <- all_locs %>% 
  rename(OffspringID=ID)
#note - BestConfig files were reformatted to be tab delimited and 
#special characters in the file were removed prior to load
#identify locations with multiple inferred cohorts
all_locs %>% 
  group_by(samp) %>% 
  summarise(nclust=length(unique(clust)),ss=n(),max_len = max(Length))
locs <- c("CAT","MAN","MIR","TWO")
best_config <- data.frame(matrix(ncol=5,nrow = 0))
#read in pedigree data for locations with multiple inferred cohorts
for (i in 1:length(locs)) {
  print(i)
  df <- readLines(paste0("SoftwareOutput/",locs[i],".Output.data.BestCluster"))
  #separate file into usable data frame
  df <- strsplit(df,"\\s+")
  df1 <- matrix(unlist(df),byrow = T)
  df1 <- df1[df1 != ""]
  df1 <- matrix(df1,ncol = 5,byrow = T)
  df1 <- as.data.frame(df1)
  colnames(df1) <- df1[1,]
  df1 <- df1[-1,]

  best_config <- rbind(best_config,df1)
}
df <- merge(best_config,all_locs)
#adding ages to model results
df %>% 
  group_by(loc,clust) %>% 
  summarise(ss=n(),max_len = max(Length))

df <- merge(all_locs,best_config)
##splitting pedigrees by collection and length cluster
#MIR
mir.1 <- subset(df,df$samp == "MIR_2017" & df$clust == "clust1")
mir.1$cohort <- "2016"
mir.2 <- subset(df,df$samp == "MIR_2017" & df$clust == "clust2")
mir.2$cohort <- "2015"
mir.3 <- subset(df,df$samp == "MIR_2017" & df$clust == "clust3")
mir.3$cohort <- "2014"
mir <- rbind(mir.1,mir.2,mir.3)
#MAN
man.1 <- subset(df,df$samp == "MAN_2019" & df$clust == "clust1")
man.1$cohort <- "2018"
man.2 <- subset(df,df$samp == "MAN_2019" & df$clust == "clust2")
man.2$cohort <- "2017"
man <- rbind(man.1,man.2)
#TWO
two.1 <- subset(df,df$samp == "TWO_2019" & df$clust == "clust2")
two.1$cohort <- "2018"
two.2 <- subset(df,df$samp == "TWO_2019" & df$clust == "clust3")
two.2$cohort <- "2017"
two <- rbind(two.1,two.2)

##quantifying family relationships across clusters
#MIR
#testing overlap
table(mir$ClusterIndex)
table(mir.1$ClusterIndex%in%mir.2$ClusterIndex)
table(mir.1$ClusterIndex%in%mir.3$ClusterIndex)
table(mir.2$ClusterIndex%in%mir.3$ClusterIndex)
bmr_sing <- mir.1[!(mir.1$ClusterIndex %in% mir.2$ClusterIndex),]
bmr_sing <- rbind(bmr_sing,bmr18[bmr18$ClusterIndex %in% bmr_sing$ClusterIndex,])

#comparing families
table(mir.1$ClusterIndex)
table(mir.2$ClusterIndex)
table(mir.3$ClusterIndex)

ggplot(mir,aes(x=ClusterIndex,y=Length,color=cohort))+
  geom_point()+theme_bw()
##man
#testing overlap
table(man$ClusterIndex)
table(man.1$ClusterIndex%in%man.2$ClusterIndex)
#comparing families
table(man.1$ClusterIndex)
table(man.2$ClusterIndex)

ggplot(man,aes(x=ClusterIndex,y=Length,color=cohort))+
  geom_point()+theme_bw()

##two
#testing overlap
table(two$ClusterIndex)
table(two.2$ClusterIndex%in%two.1$ClusterIndex)
#comparing families
table(two.1$ClusterIndex)
table(two.2$ClusterIndex)

ggplot(two,aes(x=ClusterIndex,y=Length,color=cohort))+
  geom_point()+theme_bw()

##family diagrams
all_families1 <- all_families %>% 
  select(OffspringID,MotherID,FatherID,ClusterIndex,clust,samp,cohort)
#saving family data
save(all_families,file = "Aging_Models/Family_data_all_locations.rda")
#pedigree visualization plots
bmr.plot <- subset(all_families1,all_families1$samp=="BMR_2017"|all_families1$samp=="BMR_2018")
ocq.plot <- subset(all_families1,all_families1$samp == "OCQ_2018")
che.plot <- subset(all_families1,all_families1$samp == "PR_2018")
tiff(filename = "Figures/Pedigree_plots.tiff",width = 10,height = 8,units = "in",res = 200)
par(mfrow=c(1,2))
pedigree.plot(bmr.plot,title = "Lower Black Mallard River")
pedigree.plot(ocq.plot,title = "Ocqueoc River")
pedigree.plot(che.plot)
dev.off()

##bayesmix figure
all_families$samp <- factor(all_families$samp,
                            levels = c("BMR_2017","BMR_2018","BMR_2019","OCQ_2018","PR_2018"),
                            labels = c("Lower Black Mallard - 2017 Collection","Lower Black Mallard - 2018 Collection","Upper Black Mallard - 2019 Collection","Ocqueoc River - 2018 Collection","Pigeon River - 2018 Collection"))
tiff(filename = "Figures/Age_classification_plot.tiff",width = 4,height = 5,units = "in",res = 400)
ggplot(all_families, aes(x=Length, fill = clust)) +
  facet_wrap(vars(samp),scales = "free_y",ncol = 1)+
  geom_histogram(aes(fill = factor(clust)),binwidth = 2) +
  scale_fill_manual(values = c("#225ea8","#02818a","#a8ddb5"),
                    guide = F)+
  labs(x="Length (mm)", y="counts")+
  theme_bw(base_size = 10)+
  theme(strip.background = element_rect(color="#ffffff", fill="#ffffff", size=1.5, linetype="solid"),
        plot.title = element_text(size = 12),
        panel.spacing=unit(1, "lines"))
dev.off()
