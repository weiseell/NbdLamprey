#NbdLamprey - Script 3b: Cohort creation
#Objective: Comparing length-based aging models and reconstructed pedigrees

#libraries
library(tidyverse)

#functions
source("Homebrew/pedigree.plot.R")
##reading in data
#load in data
all_locs <- read.table("Aging_Models/lw_Bayes_assignments.txt",header = T,sep = "\t",stringsAsFactors = F)
load("Software_outputs/test_best_config.rda")

#adding ages to model results
#all_families <- rbind(ocq_old,chePR,bmrBL,bmrAL)
all_locs <- all_locs %>% 
  rename(OffspringID=ID)
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

all_locs1 <- all_locs[which(all_locs$OffspringID %in% best_config$OffspringID),]
all_locs1 <- all_locs1 %>% select(-loc)

df <- merge(all_locs1,best_config)

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
ocq.2 <- subset(df,df$samp == "OCQ_2018" & df$clust == "clust2")

##quantifying family relationships across clusters
#OCQ
table(ocq.1$ClusterIndex)
table(ocq.2$ClusterIndex)
table(ocq.2$ClusterIndex %in% ocq.1$ClusterIndex)

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

##putting all locs together
bmr_cohort15$loc <- "bmrBL15"
bmr_sing$loc <- "bmrBL16"
bmral <- subset(df,df$loc == "bmrAL")
che <- subset(df,df$loc == "chePR")
all_families <- rbind(bmr_sing,bmr_cohort15,ocq,bmral,che)

all_families[which(all_families$samp == "CHE_2018"),]$samp <- "PR_2018"
##bayesmix figure
all_families$samp <- factor(all_families$samp,
                            levels = c("BMR_2017","BMR_2018","BMR_2019","OCQ_2018","PR_2018"),
                            labels = c("Lower Black Mallard - 2017 Collection","Lower Black Mallard - 2018 Collection","Upper Black Mallard","Ocqueoc River","Pigeon River"))
mix_Bayes <- ggplot(all_families, aes(x=Length, fill = clust)) +
  facet_wrap(vars(samp),scales = "free_y",ncol = 1)+
  geom_histogram(aes(fill = factor(clust)),binwidth = 2) +
  scale_fill_manual(values = c("#7fcdbb","#2c7fb8","#253494"),
                    guide = F)+
  labs(x="Length (mm)", y="counts")+
  theme_bw(base_size = 8)+
  ggtitle("A) Age Classifications - Bayesian")+
  theme(strip.background = element_rect(color="#ffffff", fill="#ffffff", size=1.5, linetype="solid"),
        plot.title = element_text(size = 12))

all_mix <- rbind(bmr17_mix[[2]],bmr18_mix[[2]],bmr19_mix[[2]],ocq_mix[[2]],che_mix[[2]])
all_mix$OffspringID <- paste(all_mix$species,all_mix$loc,all_mix$num,sep = "_")
all_mix <- all_mix %>% 
  select(-Sample_number:-V2) %>% 
  select(OffspringID,everything())
all_mix <- merge(all_mix,all_families)

mixML <- ggplot(all_mix, aes(x=Length, fill = class)) +
  facet_wrap(vars(samp),scales = "free_y",ncol = 1)+
  geom_histogram(aes(fill = factor(class)),binwidth = 2) +
  scale_fill_manual(values = c("#7fcdbb","#2c7fb8","#253494","#addd8e"),
                    guide = F)+
  labs(x="Length (mm)", y="counts")+
  theme_bw(base_size = 8)+
  ggtitle("B) Age Classifications - Maximum Likelihood")+
  theme(strip.background = element_rect(color="#ffffff", fill="#ffffff", size=1.5, linetype="solid"),
        plot.title = element_text(size = 12))
multiplot(mix_Bayes,mixML,cols = 2)

##family diagrams
all_families1 <- all_families %>% 
  select(OffspringID,MotherID,FatherID,ClusterIndex,cohort,loc)
bmr <- subset(all_families1,all_families1$loc=="bmrBL16"|all_families1$loc=="bmrBL15")

par(mfrow=c(1,2))
pedigree.plot(bmr,title = "Lower Black Mallard River")
