#NbdLamprey - Script 3b: Cohort creation
#Objective: Comparing length-based aging models and reconstructed pedigrees

#libraries
library(tidyverse)

#functions
source("Homebrew/plotPedigree.R")
##reading in data
#load in data
all_locs <- read.table("Aging_Models/lw_Bayes_assignments.txt",header = T,sep = "\t",stringsAsFactors = F)
#note - BestConfig files were reformatted to be tab delimited and 
#special characters in the file were removed prior to load
bmr <- read.table("Software_outputs/bmr_BestConfig_091820.txt",header = T,sep = "\t",stringsAsFactors = F)
che <- read.table("Software_outputs/che_BestConfig_091820.txt",header = T,sep = "\t",stringsAsFactors = F)
ocq <- read.table("Software_outputs/ocq_BestConfig_091820.txt", header = T, sep = "\t",stringsAsFactors = F)
best_config <- rbind(bmr,ocq,che)

#adding ages to model results
all_locs %>% 
  group_by(samp,clust) %>% 
  summarise(maxi <- max(Length))

#converting clusters into inferred cohorts
all_locs$cohort <- NA
all_locs[which(all_locs$samp == "BMR_2017" & all_locs$clust == "clust1"),]$cohort <- "2016"
all_locs[which(all_locs$samp == "BMR_2017" & all_locs$clust == "clust2"),]$cohort <- "2015"

all_locs[which(all_locs$samp == "BMR_2018" & all_locs$clust == "clust1"),]$cohort <- "2016"
all_locs[which(all_locs$samp == "BMR_2018" & all_locs$clust == "clust2"),]$cohort <- "2015"
all_locs[which(all_locs$samp == "BMR_2018" & all_locs$clust == "clust3"),]$cohort <- "2014"

all_locs[which(all_locs$samp == "CHE_2018" & all_locs$clust == "clust1"),]$cohort <- "2016"
all_locs[which(all_locs$samp == "CHE_2018" & all_locs$clust == "clust2"),]$cohort <- "2015"

all_locs[which(all_locs$samp == "OCQ_2018" & all_locs$clust == "clust1"),]$cohort <- "2016"

all_locs[which(all_locs$samp == "BMR_2019" & all_locs$clust == "clust2"),]$cohort <- "2016"

all_locs1 <- all_locs[which(all_locs$OffspringID %in% best_config$OffspringID),]
all_locs1 <- all_locs1 %>% select(-loc)

df <- merge(all_locs,best_config)

##splitting pedigrees by collection and length cluster
#BMR - 2017
bmr17 <- subset(df,df$samp == "BMR_2017")
bmr17.1 <- subset(df,df$samp == "BMR_2017" & df$clust == "clust1")
bmr17.2 <- subset(df,df$samp == "BMR_2017" & df$clust == "clust2")
#BMR - 2018
bmr18 <- subset(df,df$samp == "BMR_2018")
bmr18.1 <- subset(df,df$samp == "BMR_2018" & df$clust == "clust1")
bmr18.2 <- subset(df,df$samp == "BMR_2018" & df$clust == "clust2")
bmr18.3 <- subset(df,df$samp == "BMR_2018" & df$clust == "clust3")
#OCQ
ocq <- subset(df,df$samp == "OCQ_2018")
ocq.1 <- subset(df,df$samp == "OCQ_2018" & df$clust == "clust1")

#CHE - 2018
che <- subset(df,df$samp == "CHE_2018")
che.1 <- subset(df,df$samp == "CHE_2018" & df$clust == "clust1")
che.2 <- subset(df,df$samp == "CHE_2018" & df$clust == "clust2")
##quantifying family relationships across clusters
#OCQ
table(ocq.1$ClusterIndex)

#BMR
#BMR17 small group
table(bmr17.1$ClusterIndex%in%bmr17.2$ClusterIndex)
table(bmr17.1$ClusterIndex%in%bmr18.1$ClusterIndex)
table(bmr17.1$ClusterIndex%in%bmr18.2$ClusterIndex)
table(bmr17.1$ClusterIndex%in%bmr18.3$ClusterIndex)
bmr_sing <- bmr17.1[!(bmr17.1$ClusterIndex %in% bmr17.2$ClusterIndex),]
bmr_sing <- rbind(bmr_sing,bmr18[bmr18$ClusterIndex %in% bmr_sing$ClusterIndex,])

#BMR18 small group
table(bmr18.1$ClusterIndex)
table(bmr18.2$ClusterIndex)
table(bmr18.3$ClusterIndex)

##separating outlier groups
bmr_cohort16 <- bmr_sing
bmr <- subset(df,df$samp == "BMR_2017" | df$samp == "BMR_2018")
bmr_cohort15 <- bmr[!(bmr$OffspringID %in% bmr_cohort16$OffspringID),]

##separating out Pigeon River CHE samples
chePR <- che %>% 
  mutate(ID1 = OffspringID) %>% 
  separate(ID1,into = c("spp","loc1","num"),sep = "_") %>% 
  mutate(num1 = as.numeric(num)) %>% 
  filter(num1 <= 29) %>% 
  filter(num1 != 25) %>% 
  select(-spp:-num1)
chePR$samp <- "PR_2018"
##putting all locs together
bmr_cohort15$loc <- "bmrBL15"
bmr_sing$loc <- "bmrBL16"
bmral <- subset(df,df$samp == "BMR_2019")
bmral$loc <- "bmrAL"
table(bmral$clust)
bmral$clust <- "clust1"
ocq <- subset(df,df$loc == "OCQ")
all_families <- rbind(bmr_sing, bmr_cohort15, ocq, bmral, chePR)

##family diagrams
all_families1 <- all_families %>% 
  select(OffspringID,MotherID,FatherID,ClusterIndex,cohort,samp)
bmr.plot <- subset(all_families1,all_families1$samp=="BMR_2017"|all_families1$samp=="BMR_2018")
ocq.plot <- subset(all_families1,all_families1$samp == "OCQ_2018")

tiff(filename = "Figures/Pedigree_plots.tiff",width = 6,height = 4,units = "in",res = 400)
par(mfrow=c(1,2))
pedigree.plot(bmr.plot,title = "Lower Black Mallard River")
pedigree.plot(ocq.plot,title = "Ocqueoc River")
dev.off()

##bayesmix figure
all_families$samp <- factor(all_families$samp,
                            levels = c("BMR_2017","BMR_2018","BMR_2019","OCQ_2018","PR_2018"),
                            labels = c("Lower Black Mallard - 2017 Collection","Lower Black Mallard - 2018 Collection","Upper Black Mallard","Ocqueoc River","Pigeon River"))
tiff(filename = "Figures/Age_classification_plot.tiff",width = 4,height = 5,units = "in",res = 400)
ggplot(all_families, aes(x=Length, fill = clust)) +
  facet_wrap(vars(samp),scales = "free_y",ncol = 1)+
  geom_histogram(aes(fill = factor(clust)),binwidth = 2) +
  scale_fill_manual(values = c("#7fcdbb","#2c7fb8","#253494"),
                    guide = F)+
  labs(x="Length (mm)", y="counts")+
  theme_bw(base_size = 10)+
  theme(strip.background = element_rect(color="#ffffff", fill="#ffffff", size=1.5, linetype="solid"),
        plot.title = element_text(size = 12),
        panel.spacing=unit(1, "lines"))
dev.off()
