#NbdLamprey - Script 3b: Cohort creation
#Objective: Comparing length-based aging models and reconstructed pedigrees

#libraries
library(tidyverse)

#functions
source("Homebrew/pedigree.plot.R")
##reading in data
#load in data
all_locs <- read.table("Aging_Models/lw_Bayes_assignments.txt",header = T,sep = "\t",stringsAsFactors = F)
#note - BestConfig files were reformatted to be tab delimited and 
#special characters in the file were removed prior to load
#identify locations with multiple inferred cohorts
all_locs %>% 
  group_by(samp) %>% 
  summarise(nclust=length(unique(clust)),ss=n(),max_len = max(Length))
locs <- c("BMR","CHE")
best_config <- data.frame(matrix(ncol=5,nrow = 0))
#read in pedigree data for locations with multiple inferred cohorts
for (i in 1:length(locs)) {
  print(i)
  df <- readLines(paste0("Software_outputs/",locs[i],"_Output.data.BestCluster"))
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

#adding ages to model results
df <- merge(best_config,all_locs)
df$full_sib <- paste(df$MotherID,df$FatherID,sep = "_")
##splitting pedigrees by collection and length cluster
#BMR 2018
bmr18.1 <- subset(df,df$samp == "BMR_2018" & df$clust == "clust1")
bmr18.1$cohort <- "2016"
bmr18.2 <- subset(df,df$samp == "BMR_2018" & df$clust == "clust2")
bmr18.2$cohort <- "2015"
bmr18.3 <- subset(df,df$samp == "BMR_2018" & df$clust == "clust3")
bmr18.3$cohort <- "2014"
bmr18 <- rbind(bmr18.1,bmr18.2,bmr18.3)
#BMR 2017
bmr17.1 <- subset(df,df$samp == "BMR_2017" & df$clust == "clust1")
bmr17.1$cohort <- "2016"
bmr17.2 <- subset(df,df$samp == "BMR_2017" & df$clust == "clust2")
bmr17.2$cohort <- "2015"
bmr17 <- rbind(bmr17.1,bmr17.2)
bmr <- rbind(bmr17,bmr18)
bmr19 <- subset(df,df$samp == "BMR_2019")
ocq <- subset(df,df$samp == "OCQ_2018")
che <- subset(df,df$samp == "CHE_2018")
chePR <- che[che$Sample_number < 542 | che$Sample_number == 544,]
##quantifying family relationships across clusters
#BMR 18
#testing overlap
table(bmr18$ClusterIndex)
table(bmr18.1$ClusterIndex%in%bmr18.2$ClusterIndex)
table(bmr18.1$ClusterIndex%in%bmr18.3$ClusterIndex)
table(bmr18.2$ClusterIndex%in%bmr18.3$ClusterIndex)
bmr_sing <- bmr18.1[!(bmr18.1$ClusterIndex %in% bmr18.2$ClusterIndex),]
bmr_sing <- rbind(bmr_sing,bmr18[bmr18$ClusterIndex %in% bmr_sing$ClusterIndex,])

#comparing families
table(bmr18.1$ClusterIndex)
table(bmr18.2$ClusterIndex)
table(bmr18.3$ClusterIndex)

ggplot(bmr18,aes(x=ClusterIndex,y=Length,color=cohort))+
  geom_point()+theme_bw()

ggplot(bmr18,aes(x=full_sib,y=Length,color=cohort))+
  geom_point()+theme_bw()

#compare full-sibling groups
table(bmr18$full_sib)
table(bmr18.1$full_sib)
table(bmr18.2$full_sib)
table(bmr18.3$full_sib)
table(bmr18.1$full_sib%in%bmr18.2$full_sib)
table(bmr18.1$full_sib%in%bmr18.3$full_sib)
table(bmr18.3$full_sib%in%bmr18.2$full_sib)

#BMR 17
#testing overlap
table(bmr17$ClusterIndex)
table(bmr17.1$ClusterIndex%in%bmr17.2$ClusterIndex)
table(bmr17.2$ClusterIndex%in%bmr17.1$ClusterIndex)
bmr_sing <- bmr17.1[!(bmr17.1$ClusterIndex %in% bmr17.2$ClusterIndex),]
bmr_sing <- rbind(bmr_sing,bmr17[bmr17$ClusterIndex %in% bmr_sing$ClusterIndex,])

#comparing families
table(bmr17.1$ClusterIndex)
table(bmr17.2$ClusterIndex)

ggplot(bmr17,aes(x=ClusterIndex,y=Length,color=cohort))+
  geom_point()+theme_bw()

ggplot(bmr17,aes(x=full_sib,y=Length,color=cohort))+
  geom_point()+theme_bw()

#compare full-sibling groups
table(bmr17$full_sib)
table(bmr17.1$full_sib)
table(bmr17.2$full_sib)
table(bmr17.1$full_sib%in%bmr17.2$full_sib)
table(bmr17.2$full_sib%in%bmr17.1$full_sib)

#generate cohort sets
bmr_cohort16 <- bmr17.1
bmr_cohort16$cohort <- "BMR_2016"
bmr_cohort15 <- rbind(bmr17.2,bmr18)
bmr_cohort15$cohort <- "BMR_2015"
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
