#Length-Based Age Models
#1. EM models
#2. Bayes Cluster determining models
#3. Bayes models
#4. Plots
#Written by: Ellie Weise
#Last edited: 06/25/20

#loading libraries
library(tidyverse)
library(bayesmix)
library(bmixture)

#loading homebrew functions

#load in length and weight data
df <- read.table("Input_fulldata/exp_lengths_weights.txt",header = T,sep = "\t",stringsAsFactors = F)

#manipulating data frame
df1 <- df %>% 
  mutate(ID2 = ID_indiv) %>% 
  rename(ID=ID_indiv) %>% 
  separate(ID2,into = c("spp","loc","num"),sep = "_") %>% 
  select(ID,loc,everything()) %>% 
  select(-spp:-num) %>% 
  mutate(samp = paste(loc,Year_collect,sep = "_")) %>% 
  filter(samp != "OCQ_2019")
samples <- unique(df1$samp)

##2. Bayes Cluster determining models ####
RMclust <- vector(mode = "list", length = length(samples))
BDclust <- vector(mode = "list", length = length(samples))
names(RMclust) <- samples
names(BDclust) <- samples
i <- 4
for(i in 1:length(samples)){
  samp_i <- samples[i]
  tmp <- subset(df1, df1$samp == samp_i)
  len <- tmp$Length
  
  #running BayesMix models with RM criteria
  #check thinning parameter to not save all iterations
  #write models as an .Rdata file instead of a regular R file?
  model_tmp <- BMMmodel(len, k = 10, 
                        priors = list(kind = "independence",
                                      parameter = "priorsUncertain"))
  control <- JAGScontrol(variables = c("mu", "tau", "eta", "S"),
                         burn.in = 10000, n.iter = 500000)
  z <- JAGSrun(len, model = model_tmp, control = control)
  RMclust[[samp_i]] <- z
  
  #BD-MCMC with bmixture
  bmixt.model <- bmixt(len, k = "unknown", iter = 500000, burnin = 10000, k_max = 4)
  BDclust[[samp_i]] <- bmixt.model
}

#calculating probabilities
#BMR - 2017 (2 clusters)
eta.thresh <- 0.035
post.k <- rowSums(RMclust[["BMR_2017"]]$results[,grep("eta", colnames(RMclust[["BMR_2017"]]$results))] > eta.thresh)
post.k <- table(post.k)/sum(table(post.k))
post.k

#BMR - 2018 (2 clusters)
eta.thresh <- 0.035
post.k <- rowSums(RMclust[["BMR_2018"]]$results[,grep("eta", colnames(RMclust[["BMR_2018"]]$results))] > eta.thresh)
post.k <- table(post.k)/sum(table(post.k))
post.k

#BMR - 2019 (2 clusters (not converged))
eta.thresh <- 0.035
post.k <- rowSums(RMclust[["BMR_2019"]]$results[,grep("eta", colnames(RMclust[["BMR_2019"]]$results))] > eta.thresh)
post.k <- table(post.k)/sum(table(post.k))
post.k

#OCQ - 2018 (1 cluster)
eta.thresh <- 0.035
post.k <- rowSums(RMclust[["OCQ_2018"]]$results[,grep("eta", colnames(RMclust[["OCQ_2018"]]$results))] > eta.thresh)
post.k <- table(post.k)/sum(table(post.k))
post.k

#CHE - 2018 (2 clusters (not converged))
eta.thresh <- 0.035
post.k <- rowSums(RMclust[["CHE_2018"]]$results[,grep("eta", colnames(RMclust[["CHE_2018"]]$results))] > eta.thresh)
post.k <- table(post.k)/sum(table(post.k))
post.k

#OCQ_2018: 2 clusters (not good convergence)
bmix_vals <- data.frame(kval = BDclust[["OCQ_2018"]]$all_k, weight = BDclust[["OCQ_2018"]]$all_weights,stringsAsFactors = F)
bmix_vals %>% group_by(kval) %>% summarise(weightsum=sum(weight),prop = sum(weight)/sum(bmix_vals$weight))

#CHE_2018: 3 clusters (not good convergence)
bmix_vals <- data.frame(kval = BDclust[["CHE_2018"]]$all_k, weight = BDclust[["CHE_2018"]]$all_weights,stringsAsFactors = F)
bmix_vals %>% group_by(kval) %>% summarise(weightsum=sum(weight),prop = sum(weight)/sum(bmix_vals$weight))

#BMR_2017: 2 clusters (not good convergence)
bmix_vals <- data.frame(kval = BDclust[["BMR_2017"]]$all_k, weight = BDclust[["BMR_2017"]]$all_weights,stringsAsFactors = F)
bmix_vals %>% group_by(kval) %>% summarise(weightsum=sum(weight),prop = sum(weight)/sum(bmix_vals$weight))

#BMR_2018: 4 clusters (not good convergence)
bmix_vals <- data.frame(kval = BDclust[["BMR_2018"]]$all_k, weight = BDclust[["BMR_2018"]]$all_weights,stringsAsFactors = F)
bmix_vals %>% group_by(kval) %>% summarise(weightsum=sum(weight),prop = sum(weight)/sum(bmix_vals$weight))

#BMR_2019: 3 clusters (not good convergence)
bmix_vals <- data.frame(kval = BDclust[["BMR_2019"]]$all_k, weight = BDclust[["BMR_2019"]]$all_weights,stringsAsFactors = F)
bmix_vals %>% group_by(kval) %>% summarise(weightsum=sum(weight),prop = sum(weight)/sum(bmix_vals$weight))


##3. Bayes models for all locations ####
i <- 4
bestk <- c(1,1,2,2,1)
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
                                        parameter = "priorsUncertain"))
    control <- JAGScontrol(variables = c("mu", "tau", "eta", "S"),
                           burn.in = 10000, n.iter = 500000)
    
    z <- JAGSrun(len, model = model_tmp, control = control)
    Bmodels[[samp_i]] <- z
  }
}

#sorting for models with multiple cohorts
Bmodels[["BMR_2017"]] <- Sort(Bmodels[["BMR_2017"]],by = "mu")
Bmodels[["BMR_2018"]] <- Sort(Bmodels[["BMR_2018"]],by = "mu")

BayesAssign <- vector(mode = "list", length = length(samples))
names(BayesAssign) <- samples
i <- 4
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
    tmp$clust <- "clust0"
  }
  BayesAssign[[samp_i]] <- tmp
}

#saving individual assignments
all_locs_Bayes <- rbind(BayesAssign[[1]],BayesAssign[[2]],BayesAssign[[3]],BayesAssign[[4]],BayesAssign[[5]])
write.table(all_locs,file = "Aging_Models/lw_Bayes_assignments.txt",sep = "\t",row.names = F,col.names = T,quote = F)

##4. Comparing goodness of fit for chosen k ####
#BMR - 2018
norm1 <- data.frame(set="mu1",len=rnorm(n=667/2,mean = 88.45,sd = 2.0647))
norm2 <- data.frame(set="mu2",len=rnorm(n=667/2,mean = 87.47,sd = 0.5331))
normdata <- rbind(norm1,norm2)
ggplot(normdata,aes(x=len,color=set))+
  geom_histogram()

#getting sum of squares for each cluster designation
tmp <- BayesAssign[["BMR_2018"]]
tmp.1 <- subset(tmp,tmp$clust == "clust1")
tmp.2 <- subset(tmp,tmp$clust == "clust2")


##4. Plotting all models ####
all_locs_Bayes <- read.table("Aging_Models/lw_Bayes_assignments.txt",header = T)
all_locs_Bayes <- subset(all_locs_Bayes,all_locs_Bayes$samp != "BMR_2018")
#BayesAssign[["BMR_2018"]] <- BayesAssign[["BMR_2018"]] %>% rename(OffspringID=ID)
all_locs_Bayes <- rbind(all_locs_Bayes,BayesAssign[["BMR_2018"]])
all_locs_Bayes$clust[which(all_locs_Bayes$samp == "CHE_2018")] <- "clust0"
all_locs_Bayes$clust[which(all_locs_Bayes$samp == "BMR_2019")] <- "clust0"
#labels for plot
all_locs_Bayes$samp <- factor(all_locs_Bayes$samp,
                            levels = c("BMR_2017","BMR_2018","BMR_2019","OCQ_2018","CHE_2018"),
                            labels = c("Lower Black Mallard - 2017 Collection","Lower Black Mallard - 2018 Collection","Upper Black Mallard - 2019 Collection","Ocqueoc River - 2018 Collection","Pigeon River - 2018 Collection"))
tiff(filename = "Figures/Length_Histogram_plot.tiff",width = 4,height = 6,units = "in",res = 400)
ggplot(all_locs_Bayes, aes(x=Length, fill = clust)) +
  facet_wrap(~samp, scales = "free_y",ncol = 1) +
  geom_histogram(aes(fill = factor(clust)),bins = 50) +
  scale_fill_manual(values = c("grey","#74a9cf","#034e7b"),
                    guide = F)+
  labs(x="Length (mm)",y="counts")+
  theme_bw(base_size = 8)
dev.off()

ggplot(BayesAssign[["BMR_2018"]], aes(x=Length, fill = clust)) +
  facet_wrap(~samp, scales = "free_y") +
  geom_histogram(aes(fill = factor(clust)),bins = 50) +
  scale_fill_manual(values = c("#000000","#cccccc","#969696"),
                    guide = F)+
  labs(x="Length (mm)",y="counts")+
  theme_bw(base_size = 8)








