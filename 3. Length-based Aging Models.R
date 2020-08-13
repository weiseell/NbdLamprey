#Length-Based Age Models
#1. EM models
#2. Bayes Cluster determining models
#3. Bayes models
#4. Plots
#Written by: Ellie Weise
#Last edited: 06/25/20

#loading libraries
library(tidyverse)
library(mclust)
library(bayesmix)
library(bmixture)

#setwd
setwd("~/OneDrive - Michigan State University/Documents/Sea_Lamprey_MS_project/Bioinformatics/BayesMix/")

#loading homebrew functions
source("Homebrew/mixture_function.R")

#load in length and weight data
load("Input/test_length_weight.rda")
df <- lw

#manipulating data frame
df1 <- df %>% 
  mutate(ID2 = ID_indiv) %>% 
  rename(ID=ID_indiv) %>% 
  separate(ID2,into = c("spp","loc","num"),sep = "_") %>% 
  select(ID,loc,everything()) %>% 
  select(-spp:-num) %>% 
  mutate(samp = paste(loc,Year_collect,sep = "_"))
samples <- unique(df1$samp)

##1. EM models for all locations ####
EMclusters <- vector(mode = "list", length = length(samples))
EMmodels <- vector(mode = "list", length = length(samples))
names(EMclusters) <- samples
names(EMmodels) <- samples
df2 <- df1 %>% rename(V1=Length,V2=Weight)

#loop to determine the number of clusters
i <- 1
for (i in 1:length(samples)) {
  samp_i <- samples[i]
  tmp <- subset(df2, df2$samp == samp_i)
  lw <- as.matrix(data.frame(tmp$V1, tmp$V2))
  
  #calculating BIC for 1-4 clusters
  model <- mclustBIC(lw, G = 1:4, modelNames = c("VVV"))
  EMclusters[[samp_i]] <- model
  
}

#looking at BIC values
EMclusters[["BMR_2017"]] #best model: 4 clusters
EMclusters[["BMR_2018"]] #best model: 3 clusters
EMclusters[["BMR_2019"]] #best model: 1 cluster

#running EM models
#loop to determine the number of clusters
i <- 1
for (i in 1:length(samples)) {
  samp_i <- samples[i]
  tmp <- subset(df2, df2$samp == samp_i)
  lw <- data.frame(V1=tmp$V1, V2=tmp$V2,stringsAsFactors = F)
  
  #models
  model <- mixture(lw,maxclust = EMclusters[[samp_i]]) 
  EMmodels[[samp_i]] <- model
  tmp1 <- data.frame(model[[2]],stringsAsFactors = F)
  tmp1$pop <- samp_i
  model[[2]] <- tmp1
  
}


##2. Bayes Cluster determining models ####
RMclust <- vector(mode = "list", length = length(samples))
Bclust <- vector(mode = "list", length = length(samples))
BDclust <- vector(mode = "list", length = length(samples))
names(RMclust) <- samples
names(Bclust) <- samples
names(BDclust) <- samples
i <- 1
for(i in 1:length(samples)){
  samp_i <- samples[i]
  tmp <- subset(df1, df1$samp == samp_i)
  len <- tmp$Length
  
  #running BayesMix models with RM criteria
  model_tmp <- BMMmodel(len, k = 10, 
                        priors = list(kind = "independence",
                                      parameter = "priorsUncertain"))
  control <- JAGScontrol(variables = c("mu", "tau", "eta", "S"),
                         burn.in = 10000, n.iter = 50000)
  z <- JAGSrun(len, model = model_tmp, control = control)
  #calculating probabilities
  eta.thresh <- 0.035
  post.k <- rowSums(z$results[,grep("eta", colnames(z$results))] > eta.thresh)
  post.k <- table(post.k)/sum(table(post.k))
  post.k
  
  RMclust[[samp_i]] <- post.k
  
  #BayesMix models for BMMdiag plots
  model_tmp <- BMMmodel(len,k = 4, initialValues = list(S0 = 2),
                        priors = list(kind = "independence",
                                      parameter = "priorsFish", 
                                      hierarchical = "tau"))
  control <- JAGScontrol(variables = c("mu", "tau", "eta", "S"),
                         burn.in = 10000, n.iter = 50000) 
  z <- JAGSrun(len, model = model_tmp, control = control)
  Bclust[[samp_i]] <- z
  
  #BD-MCMC with bmixture
  bmixt.model <- bmixt(len, k = "unknown", iter = 500000, burnin = 100000, k_max = 4)
  BDclust[[samp_i]] <- bmixt.model
}

RMclust[["BMR_2017"]] #best: 2 clusters
RMclust[["BMR_2018"]] #best: 3 clusters
RMclust[["BMR_2019"]] #best: 2 clusters

BMMdiag(Bclust[["BMR_2017"]]) #best cluster: 2
BMMdiag(Bclust[["BMR_2018"]]) #best cluster: 2
BMMdiag(Bclust[["BMR_2019"]]) #best cluster: 1

#BMR_2017
bmix_vals <- data.frame(kval = BDclust[["BMR_2017"]]$all_k, weight = BDclust[["BMR_2017"]]$all_weights,stringsAsFactors = F)
bmix_vals %>% group_by(kval) %>% summarise(weightsum=sum(weight),prop = sum(weight)/sum(bmix_vals$weight))
#BMR_2018
bmix_vals <- data.frame(kval = BDclust[["BMR_2018"]]$all_k, weight = BDclust[["BMR_2018"]]$all_weights,stringsAsFactors = F)
bmix_vals %>% group_by(kval) %>% summarise(weightsum=sum(weight),prop = sum(weight)/sum(bmix_vals$weight))
#BMR_2019
bmix_vals <- data.frame(kval = BDclust[["BMR_2019"]]$all_k, weight = BDclust[["BMR_2019"]]$all_weights,stringsAsFactors = F)
bmix_vals %>% group_by(kval) %>% summarise(weightsum=sum(weight),prop = sum(weight)/sum(bmix_vals$weight))


##3. Bayes models for all locations ####
i <- 1
bestk <- c(2,3,1)
Bmodels <- vector(mode = "list", length = length(samples))
names(Bmodels) <- samples

for (i in 1:length(samples)) {
  samp_i <- samples[i]
  tmp <- subset(df1, df1$samp == samp_i)
  len <- tmp$Length
  
  if(bestk[i] > 1) {
    #bayesmix model
    model_tmp <- BMMmodel(len,k = bestk[i], 
                          priors = list(kind = "independence",
                                        parameter = "priorsUncertain"),
                          restrict = "tau")
    control <- JAGScontrol(variables = c("mu", "tau", "eta", "S"),
                           burn.in = 10000, n.iter = 50000)
    
    z <- JAGSrun(len, model = model_tmp, control = control)
    Bmodels[[samp_i]] <- z
  }
}
##diagnostic plots to determine the best sorting method
BMMdiag(Bmodels[["BMR_2017"]])
BMMdiag(Bmodels[["BMR_2018"]])
BMMdiag(Bmodels[["BMR_2019"]])

#sorting for models with multiple cohorts
Bmodels[["BMR_2017"]] <- Sort(Bmodels[["BMR_2017"]],by = "mu")
Bmodels[["BMR_2018"]] <- Sort(Bmodels[["BMR_2018"]],by = "mu")

i <- 1
for (i in 1:length(samples)) {
  samp_i <- samples[i]
  tmp <- subset(df1, df1$samp == samp_i)
  len <- tmp$Length
  zSort <- Bmodels[[samp_i]]
  
  if(bestk[i] > 1) {    
    #Extracting probabilities and assigning classes
    probs <- BMMposteriori(zSort, plot = F)
    probs <- as.data.frame(cbind(probs$data,t(probs$post)),stringsAsFactors = F)
    probs$clust <- apply(probs[,-1], 1, FUN = function(x){which(x == max(x))})
    probs$clust <- paste0("clust",probs$clust)
    
    #assignments from probabilities
    probs <- probs %>% arrange(V1)
    j <- 1
    tmp$clust <- NA
    for (j in 1:length(probs$V1)) {
      tmp[which(tmp$Length == probs$V1[j]),]$clust <- probs[j,]$clust
    }
  }
  if(bestk[i] == 1){
    tmp$clust <- "clust1"
  }
  Bmodels[[samp_i]] <- tmp
}

#saving individual assignments
all_locs <- rbind(Bmodels[[1]],Bmodels[[2]],Bmodels[[3]])
write.table(all_locs,file = "Aging_Models/lw_Bayes_assignments.txt",sep = "\t",row.names = F,col.names = T,quote = F)
all_locs_EM <- rbind(EMmodels[[1]][[2]],EMmodels[[2]][[2]],EMmodels[[3]][[2]])
write.table(all_locs,file = "Aging_Models/lw_EM_assignments.txt",sep = "\t",row.names = F,col.names = T,quote = F)

##4. Plotting all models ####
ggplot(all_locs_EM, aes(x=V1, fill = class)) +
  facet_wrap(~pop, scales = "free_y") +
  geom_histogram(aes(fill = factor(class)),bins = 100) +
  scale_fill_manual(values = c("#000000","#cccccc","#969696","#636363"),
                    name = "Age",
                    labels = c("0", "1","2","3"),
                    guide = guide_legend(reverse = TRUE))+
  labs(x="Length (mm)",y="counts")+
  theme_bw(base_size = 8)

ggplot(all_locs_Bayes, aes(x=Length, fill = clust)) +
  facet_wrap(~samp, scales = "free_y") +
  geom_histogram(aes(fill = factor(clust)),bins = 100) +
  scale_fill_manual(values = c("#000000","#cccccc","#969696","#636363"),
                    name = "Age",
                    labels = c("0", "1","2","3"),
                    guide = guide_legend(reverse = TRUE))+
  labs(x="Length (mm)",y="counts")+
  theme_bw(base_size = 8)











