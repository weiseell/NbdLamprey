#Ns accumulation curves

#load libraries
library(tidyverse)
#read in functions
source("Homebrew/Ns_calc.R")
source("Homebrew/multiplot.R")
#load in data
load("Aging_Models/Family_data_all_locations.rda")
#calculate Ns and making accumulation curves
unique(all_families$loc)
##OCQ
#calculate Ns
Ns_tmp <- Ns_calc(subset(all_families,all_families$loc == "OCQ"))
#isolating data to plot
df <- data.frame(sites=Ns_tmp[[1]]$sites,richness=Ns_tmp[[1]]$richness,sd=Ns_tmp[[1]]$sd,stringsAsFactors = F)
reps <- as.data.frame(Ns_tmp[[1]]$perm)
reps1 <- as.data.frame(reps[seq(1, nrow(reps), 10),])
reps1 <- reps1 %>% 
  mutate(sites=seq(1, nrow(reps), 10)) %>% 
  select(sites,everything()) %>% 
  gather(key="rep",value = "richness",-sites)
#plotting data
ocq_plot <- ggplot(df,aes(x=sites,y=richness))+
  theme_classic()+
  geom_ribbon(aes(ymin=richness-sd,ymax=richness+sd),fill="lightgrey")+
  geom_line()+
  geom_boxplot(data = reps1,aes(group=sites),outlier.shape = NA)+
  geom_hline(yintercept=Ns_tmp[[2]]$chao)+
  geom_text(data=data.frame(x=0,y=Ns_tmp[[2]]$chao), aes(x, y), label=paste("Chao estimate = ",round(Ns_tmp[[2]]$chao,digits = 2)),hjust=0, vjust=1.5)+
  geom_hline(yintercept=Ns_tmp[[2]]$jack1)+
  geom_text(data=data.frame(x=0,y=Ns_tmp[[2]]$jack1), aes(x, y), label=paste("Jackknife estimate = ",round(Ns_tmp[[2]]$jack1,digits = 2)),hjust=0, vjust=-1)+
  ylim(0,110)+
  xlab("Number of offspring sampled")+
  ylab("Number of parent genotypes")+
  ggtitle("D) Ocqueoc River")
##BMR
#calculate Ns
Ns_tmp <- Ns_calc(subset(all_families,all_families$loc == "2015"))
#isolating data to plot
df <- data.frame(sites=Ns_tmp[[1]]$sites,richness=Ns_tmp[[1]]$richness,sd=Ns_tmp[[1]]$sd,stringsAsFactors = F)
reps <- as.data.frame(Ns_tmp[[1]]$perm)
reps1 <- as.data.frame(reps[seq(1, nrow(reps), 10),])
reps1 <- reps1 %>% 
  mutate(sites=seq(1, nrow(reps), 10)) %>% 
  select(sites,everything()) %>% 
  gather(key="rep",value = "richness",-sites)
#plotting data
bmrbl_plot1 <- ggplot(df,aes(x=sites,y=richness))+
  theme_classic()+
  geom_ribbon(aes(ymin=richness-sd,ymax=richness+sd),fill="lightgrey")+
  geom_line()+
  geom_boxplot(data = reps1,aes(group=sites),outlier.shape = NA)+
  geom_hline(yintercept=Ns_tmp[[2]]$chao)+
  geom_text(data=data.frame(x=0,y=Ns_tmp[[2]]$chao), aes(x, y), label=paste("Chao estimate = ",round(Ns_tmp[[2]]$chao,digits = 2)),hjust=0, vjust=-1)+
  geom_hline(yintercept=Ns_tmp[[2]]$jack1)+
  geom_text(data=data.frame(x=0,y=Ns_tmp[[2]]$jack1), aes(x, y), label=paste("Jackknife estimate = ",round(Ns_tmp[[2]]$jack1,digits = 2)),hjust=0, vjust=1.5)+
  ylim(0,175)+
  xlab("Number of offspring sampled")+
  ylab("Number of parent genotypes")+
  ggtitle("A) Lower Black Mallard - 2015")
#calculate Ns
Ns_tmp <- Ns_calc(subset(all_families,all_families$loc == "2016"))
#isolating data to plot
df <- data.frame(sites=Ns_tmp[[1]]$sites,richness=Ns_tmp[[1]]$richness,sd=Ns_tmp[[1]]$sd,stringsAsFactors = F)
reps <- as.data.frame(Ns_tmp[[1]]$perm)
reps1 <- as.data.frame(reps[seq(1, nrow(reps), 1),])
reps1 <- reps1 %>% 
  mutate(sites=seq(1, nrow(reps), 1)) %>% 
  select(sites,everything()) %>% 
  gather(key="rep",value = "richness",-sites)
#plotting data
bmrbl_plot2 <- ggplot(df,aes(x=sites,y=richness))+
  theme_classic()+
  geom_ribbon(aes(ymin=richness-sd,ymax=richness+sd),fill="lightgrey")+
  geom_line()+
  geom_boxplot(data = reps1,aes(group=sites),outlier.shape = NA)+
  geom_hline(yintercept=Ns_tmp[[2]]$chao)+
  geom_text(data=data.frame(x=0,y=Ns_tmp[[2]]$chao), aes(x, y), label=paste("Chao estimate = ",round(Ns_tmp[[2]]$chao,digits = 2)),hjust=0, vjust=-1.5)+
  geom_hline(yintercept=Ns_tmp[[2]]$jack1)+
  geom_text(data=data.frame(x=0,y=Ns_tmp[[2]]$jack1), aes(x, y), label=paste("Jackknife estimate = ",round(Ns_tmp[[2]]$jack1,digits = 2)),hjust=0, vjust=1.5)+
  ylim(0,30)+
  xlab("Number of offspring sampled")+
  ylab("Number of parent genotypes")+
  ggtitle("B) Lower Black Mallard - 2016")
##chePR
Ns_tmp <- Ns_calc(subset(all_families,all_families$loc == "CHE"))
#isolating data to plot
df <- data.frame(sites=Ns_tmp[[1]]$sites,richness=Ns_tmp[[1]]$richness,sd=Ns_tmp[[1]]$sd,stringsAsFactors = F)
reps <- as.data.frame(Ns_tmp[[1]]$perm)
reps1 <- as.data.frame(reps[seq(1, nrow(reps), 1),])
reps1 <- reps1 %>% 
  mutate(sites=seq(1, nrow(reps), 1)) %>% 
  select(sites,everything()) %>% 
  gather(key="rep",value = "richness",-sites)
#plotting data
che_plot <- ggplot(df,aes(x=sites,y=richness))+
  theme_classic()+
  geom_ribbon(aes(ymin=richness-sd,ymax=richness+sd),fill="lightgrey")+
  geom_line()+
  geom_boxplot(data = reps1,aes(group=sites),outlier.shape = NA)+
  geom_hline(yintercept=Ns_tmp[[2]]$chao)+
  geom_text(data=data.frame(x=0,y=Ns_tmp[[2]]$chao), aes(x, y), label=paste("Chao estimate = ",round(Ns_tmp[[2]]$chao,digits = 2)),hjust=0, vjust=-1.5)+
  geom_hline(yintercept=Ns_tmp[[2]]$jack1)+
  geom_text(data=data.frame(x=0,y=Ns_tmp[[2]]$jack1), aes(x, y), label=paste("Jackknife estimate = ",round(Ns_tmp[[2]]$jack1,digits = 2)),hjust=0, vjust=1.25)+
  ylim(0,20)+
  xlab("Number of offspring sampled")+
  ylab("Number of parent genotypes")+
  ggtitle("E) Pigeon River")

##bmrAL
Ns_tmp <- Ns_calc(subset(all_families,all_families$loc == "bmrAL"))
#isolating data to plot
df <- data.frame(sites=Ns_tmp[[1]]$sites,richness=Ns_tmp[[1]]$richness,sd=Ns_tmp[[1]]$sd,stringsAsFactors = F)
reps <- as.data.frame(Ns_tmp[[1]]$perm)
reps1 <- as.data.frame(reps[seq(1, nrow(reps), 1),])
reps1 <- reps1 %>% 
  mutate(sites=seq(1, nrow(reps), 1)) %>% 
  select(sites,everything()) %>% 
  gather(key="rep",value = "richness",-sites)
#plotting data
bmral_plot <- ggplot(df,aes(x=sites,y=richness))+
  theme_classic()+
  geom_ribbon(aes(ymin=richness-sd,ymax=richness+sd),fill="lightgrey")+
  geom_line()+
  geom_boxplot(data = reps1,aes(group=sites),outlier.shape = NA)+
  geom_hline(yintercept=Ns_tmp[[2]]$chao)+
  geom_text(data=data.frame(x=0,y=Ns_tmp[[2]]$chao), aes(x, y), label=paste("Chao estimate = ",round(Ns_tmp[[2]]$chao,digits = 2)),hjust=0, vjust=1.5)+
  geom_hline(yintercept=Ns_tmp[[2]]$jack1)+
  geom_text(data=data.frame(x=0,y=Ns_tmp[[2]]$jack1), aes(x, y), label=paste("Jackknife estimate = ",round(Ns_tmp[[2]]$jack1,digits = 2)),hjust=0, vjust=-1)+
  ylim(0,20)+
  xlab("Number of offspring sampled")+
  ylab("Number of parent genotypes")+
  ggtitle("C) Upper Black Mallard River")

tiff(filename = "Figures/Ns_extrapolation_plots.tiff",height = 10,width = 10,units = "in",res = 200)
multiplot(bmrbl_plot1,bmrbl_plot2,bmral_plot,ocq_plot,che_plot,cols = 2)
dev.off()
