library(tidyverse)
library(adegenet)
library(ape)

nl <- read.table("Input_fulldata/GTs_8X_051721.GT.FORMAT",header = T,sep = "\t",stringsAsFactors = F)
load("Input_fulldata/rapture_panel_all_SNPs.rda")
load("Summary_Stats/SNP_summaries_targets.rda")
source("Homebrew/match_tags.R")
source("Homebrew/vcf_colony.R")
rapture <- rapture %>% 
  rename(ID=rad_tag_name) %>% 
  mutate(ID1=ID) %>% 
  separate(ID1,into = c("CHROM","POS"),sep = ":") %>% 
  separate(POS,into = c("min","max"),sep = "-") %>% 
  select(ID, CHROM, min, max)

##Test PCA ####
#subset native lamprey and calculate percent genotyped
la_if <- nl[,grep(pattern = "LA",x=colnames(nl),ignore.case = T)]
la_if <- cbind(la_if,nl[,grep(pattern = "IF",x=colnames(nl),ignore.case = T)])
la_if <- data.frame(CHROM=nl$CHROM,POS=nl$POS,la_if)
colSums(la_if == "./.")
la_if <- la_if %>% select(-LA_06,-IF_01)
gt_missing <- rowSums(la_if == "./.")
la_if$pGT <- (1-(gt_missing/23))

#select SNPs with minimal missing data
nl_SNPs <- la_if %>% filter(pGT > 0.9) %>% select(CHROM,POS)
target <- match_tags(SNPs = nl_SNPs,tags = rapture)
target1 <- target[[2]]
nl_targets <- merge(target1,la_if)

#select test genotypes
test_gts <- nl[,sample(28:ncol(nl),30,replace = F)]
test_gts <- cbind(nl[,1:2],test_gts)

#get SNPs for sea lamrpey and combine with native lamprey
pm_gts <- merge(nl_SNPs,test_gts)
all_SNPs <- merge(pm_gts,nl_targets)

#convert to genlight formatting
col1 <- all_SNPs %>% 
  mutate(ID=paste(CHROM,POS,sep = "-")) %>% 
  select(ID,everything()) %>% 
  select(-CHROM:-POS) %>% 
  select(-target,-pGT)
row.names(col1) <- col1$ID
SNPs <- rownames(col1)
indiv <- colnames(col1)
indiv <- indiv[-1]
col2 <- col1 %>% select(-ID)

loci <- data.frame(ID=SNPs,stringsAsFactors = F)
loci <- loci %>% 
  separate(ID, into = c("scaffold","pos"),sep = "-")
loc <- data.frame(ID=indiv,stringsAsFactors = F)
locs <- loc %>% separate(ID, into = c("spp","loc","indiv"),sep = "_") %>% select(spp)
PM_subpops <- loc %>% separate(ID, into = c("spp","loc","indiv"),sep = "_") %>% select(loc)
#transpose matrix
col2 <- data.frame(t(col2),stringsAsFactors = F)
rownames(col2) <- indiv
colnames(col2) <- SNPs

#check amount of missing data per individual
col2 <- data.frame(lapply(col2, gsub, pattern = "./.", replacement = NA, fixed = TRUE))
col2 <- data.frame(lapply(col2, gsub, pattern = "0/0", replacement = 0, fixed = TRUE))
col2 <- data.frame(lapply(col2, gsub, pattern = "0/1", replacement = 1, fixed = TRUE))
col2 <- data.frame(lapply(col2, gsub, pattern = "1/1", replacement = 2, fixed = TRUE))
save(col2,file = "SNPsets/gentype.for.PCA.test.rda")
#convert to genlight object for PCA
snp <- new("genlight",
           col2,
           chromosome=loci$scaffold,
           position=loci$pos,
           pop=locs$spp)
indNames(snp) <- indiv
save(snp,file = "genlight_pca_test.RData")
#running PCA analysis
#getting colors
locs$col <- NA
for(i in 1:length(locs$spp)){
  if(locs$spp[i] == "LA"){locs$col[i] <- "blue"}
  if(locs$spp[i] == "IF"){locs$col[i] <- "dark green"}
  if(locs$spp[i] == "PM"){locs$col[i] <- "purple"}
}

#run pca
pca <- glPca(snp,nf = 10)
save(pca,file = "pca_test.RData")
#plot pca
plot(pca$scores[,1], pca$scores[,2],
     col=locs$col,cex=0.5)
plot(pca$scores[,1], pca$scores[,3],
     col=locs$col,cex=0.5)
#make a tree
#construct tree
tre <- nj(dist(as.matrix(snp)))
#plot tree
plot(tre, typ="fan", cex=0.7)

##run PCA for each population ####
##Black Mallard
bmr <- nl[,grep(pattern = "BMR",x=colnames(nl),ignore.case = T)]
bmr <- data.frame(CHROM=nl$CHROM,POS=nl$POS,bmr)

#get SNPs for sea lamrpey and combine with native lamprey
pm_gts <- merge(nl_SNPs,bmr)
bmrsums <- colSums(bmr == "./.")
gt_missing <- rowSums(bmr == "./.")
bmr$pGT <- (1-(gt_missing/1088))

#select SNPs with minimal missing data
bmr_SNPs <- bmr %>% filter(pGT > 0.95) %>% select(CHROM,POS)
target <- match_tags(SNPs = bmr_SNPs,tags = rapture)
target1 <- target[[2]]
bmr_targets <- merge(target1,bmr)

all_SNPs <- merge(bmr_targets,nl_targets)
all_SNPs <- bmr_targets %>% 
  full_join(nl_targets, by = c("CHROM","POS","target"))

all_SNPs[is.na(all_SNPs)] <- "0/0"

table(colSums(all_SNPs))
#remove individuals with > 0.8 proportion of missing data
all_SNPs <- all_SNPs[,-which(colSums(all_SNPs == "./.") > 30000)]
#convert to genlight formatting
col1 <- all_SNPs %>% 
  mutate(ID=paste(CHROM,POS,sep = "-")) %>% 
  select(ID,everything()) %>% 
  select(-CHROM:-POS) %>% 
  select(-target,-pGT.x,-pGT.y)
row.names(col1) <- col1$ID
SNPs <- rownames(col1)
indiv <- colnames(col1)
indiv <- indiv[-1]
col2 <- col1 %>% select(-ID)

loci <- data.frame(ID=SNPs,stringsAsFactors = F)
loci <- loci %>% 
  separate(ID, into = c("scaffold","pos"),sep = "-")
loc <- data.frame(ID=indiv,stringsAsFactors = F)
locs <- loc %>% separate(ID, into = c("spp","loc","indiv"),sep = "_") %>% select(spp)
PM_subpops <- loc %>% separate(ID, into = c("spp","loc","indiv"),sep = "_") %>% select(loc)
#transpose matrix
col2 <- data.frame(t(col2),stringsAsFactors = F)
rownames(col2) <- indiv
colnames(col2) <- SNPs

#check amount of missing data per individual
col2 <- data.frame(lapply(col2, gsub, pattern = "./.", replacement = NA, fixed = TRUE))
col2 <- data.frame(lapply(col2, gsub, pattern = "0/0", replacement = 0, fixed = TRUE))
col2 <- data.frame(lapply(col2, gsub, pattern = "0/1", replacement = 1, fixed = TRUE))
col2 <- data.frame(lapply(col2, gsub, pattern = "1/1", replacement = 2, fixed = TRUE))
save(col2,file = "SNPsets/gentype.for.PCA.BMR.rda")
#convert to genlight object for PCA
snp <- new("genlight",
           col2,
           chromosome=loci$scaffold,
           position=loci$pos,
           pop=locs$spp)
indNames(snp) <- indiv
save(snp,file = "genlight_pca_BMR.RData")
#running PCA analysis
#getting colors
locs <- data.frame(spp=snp@pop,stringsAsFactors = F)
locs <- locs %>% 
  select(spp)
locs$col <- NA
for(i in 1:length(locs$spp)){
  if(locs$spp[i] == "LA"){locs$col[i] <- "blue"}
  if(locs$spp[i] == "IF"){locs$col[i] <- "dark green"}
  if(locs$spp[i] == "PM"){locs$col[i] <- "purple"}
}

#plot pca
tiff(filename = "Figures/BMR_PCA_plot.tiff",width = 5,height = 5,units = "in",res = 400)
plot(pca$scores[,1], pca$scores[,2],
     col=locs$col,cex=0.5)
text(x = 0, y = 0, "P.marinus", pos =4, col = "purple",cex = 0.8)
text(x = 90, y = -40, "I.fosser", pos =1, col = "dark green",cex = 0.8)
text(x = 90, y = 55, "L.appendix", pos =3, col = "blue",cex = 0.8)
title("Black Mallard River")
dev.off()

#make a tree
#construct tree
tre <- nj(dist(as.matrix(snp)))
#plot tree
plot(tre, typ="fan", cex=0.7)

##Ocqueoc
ocq <- nl[,grep(pattern = "OCQ",x=colnames(nl),ignore.case = T)]
ocq <- data.frame(CHROM=nl$CHROM,POS=nl$POS,ocq)

#get SNPs for sea lamrpey and combine with native lamprey
pm_gts <- merge(nl_SNPs,ocq)
ocqsums <- colSums(ocq == "./.")
gt_missing <- rowSums(ocq == "./.")
ocq$pGT <- (1-(gt_missing/396))

#select SNPs with minimal missing data
ocq_SNPs <- ocq %>% filter(pGT > 0.95) %>% select(CHROM,POS)
target <- match_tags(SNPs = ocq_SNPs,tags = rapture)
target1 <- target[[2]]
ocq_targets <- merge(target1,ocq)

all_SNPs <- merge(ocq_targets,nl_targets)
all_SNPs <- ocq_targets %>% 
  full_join(nl_targets, by = c("CHROM","POS","target"))

all_SNPs[is.na(all_SNPs)] <- "0/0"
#remove individuals with > 0.8 proportion of missing data
all_SNPs <- all_SNPs[,-which(colSums(all_SNPs == "./.") > 30000)]
#convert to genlight formatting
col1 <- all_SNPs %>% 
  mutate(ID=paste(CHROM,POS,sep = "-")) %>% 
  select(ID,everything()) %>% 
  select(-CHROM:-POS) %>% 
  select(-target,-pGT.x,-pGT.y)
row.names(col1) <- col1$ID
SNPs <- rownames(col1)
indiv <- colnames(col1)
indiv <- indiv[-1]
col2 <- col1 %>% select(-ID)

loci <- data.frame(ID=SNPs,stringsAsFactors = F)
loci <- loci %>% 
  separate(ID, into = c("scaffold","pos"),sep = "-")
loc <- data.frame(ID=indiv,stringsAsFactors = F)
locs <- loc %>% separate(ID, into = c("spp","loc","indiv"),sep = "_") %>% select(spp)
PM_subpops <- loc %>% separate(ID, into = c("spp","loc","indiv"),sep = "_") %>% select(loc)
#transpose matrix
col2 <- data.frame(t(col2),stringsAsFactors = F)
rownames(col2) <- indiv
colnames(col2) <- SNPs

#check amount of missing data per individual
col2 <- data.frame(lapply(col2, gsub, pattern = "./.", replacement = NA, fixed = TRUE))
col2 <- data.frame(lapply(col2, gsub, pattern = "0/0", replacement = 0, fixed = TRUE))
col2 <- data.frame(lapply(col2, gsub, pattern = "0/1", replacement = 1, fixed = TRUE))
col2 <- data.frame(lapply(col2, gsub, pattern = "1/1", replacement = 2, fixed = TRUE))
save(col2,file = "SNPsets/gentype.for.PCA.ocq.rda")
#convert to genlight object for PCA
snp <- new("genlight",
           col2,
           chromosome=loci$scaffold,
           position=loci$pos,
           pop=locs$spp)
indNames(snp) <- indiv
save(snp,file = "genlight_pca_ocq.RData")
#running PCA analysis
#getting colors
locs <- data.frame(spp=snp@pop,stringsAsFactors = F)
locs <- locs %>% 
  select(spp)
locs$col <- NA
for(i in 1:length(locs$spp)){
  if(locs$spp[i] == "LA"){locs$col[i] <- "blue"}
  if(locs$spp[i] == "IF"){locs$col[i] <- "dark green"}
  if(locs$spp[i] == "PM"){locs$col[i] <- "purple"}
}

#plot pca
tiff(filename = "Figures/OCQ_PCA_plot.tiff",width = 5,height = 5,units = "in",res = 400)
plot(pca$scores[,1], pca$scores[,2],
     col=locs$col,cex=0.5)
text(x = 0, y = 0, "P.marinus", pos =4, col = "purple",cex = 0.8)
text(x = 90, y = 50, "I.fosser", pos =1, col = "dark green",cex = 0.8)
text(x = 90, y = -60, "L.appendix", pos =3, col = "blue",cex = 0.8)
title("Ocqueoc River")
dev.off()

#make a tree
#construct tree
tre <- nj(dist(as.matrix(snp)))
#plot tree
plot(tre, typ="fan", cex=0.7)

##Cheboygan
che <- nl[,grep(pattern = "CHE",x=colnames(nl),ignore.case = T)]
che <- data.frame(CHROM=nl$CHROM,POS=nl$POS,che)

#subset to only include Pigeon River samples
che <- che[,c(1:26,28:29)]
#get SNPs for sea lamrpey and combine with native lamprey
pm_gts <- merge(nl_SNPs,che)
chesums <- colSums(che == "./.")
gt_missing <- rowSums(che == "./.")
che$pGT <- (1-(gt_missing/26))

#select SNPs with minimal missing data
che_SNPs <- che %>% filter(pGT > 0.7) %>% select(CHROM,POS)
target <- match_tags(SNPs = che_SNPs,tags = rapture)
target1 <- target[[2]]
che_targets <- merge(target1,che)

all_SNPs <- che_targets %>% 
  full_join(nl_targets, by = c("CHROM","POS","target"))

all_SNPs[is.na(all_SNPs)] <- "0/0"
#remove individuals with > 0.8 proportion of missing data
all_SNPs <- all_SNPs[,-which(colSums(all_SNPs == "./.") > 30000)]
#convert to genlight formatting
col1 <- all_SNPs %>% 
  mutate(ID=paste(CHROM,POS,sep = "-")) %>% 
  select(ID,everything()) %>% 
  select(-CHROM:-POS) %>% 
  select(-target,-pGT.x,-pGT.y)
row.names(col1) <- col1$ID
SNPs <- rownames(col1)
indiv <- colnames(col1)
indiv <- indiv[-1]
col2 <- col1 %>% select(-ID)

loci <- data.frame(ID=SNPs,stringsAsFactors = F)
loci <- loci %>% 
  separate(ID, into = c("scaffold","pos"),sep = "-")
loc <- data.frame(ID=indiv,stringsAsFactors = F)
locs <- loc %>% separate(ID, into = c("spp","loc","indiv"),sep = "_") %>% select(spp)
PM_subpops <- loc %>% separate(ID, into = c("spp","loc","indiv"),sep = "_") %>% select(loc)
#transpose matrix
col2 <- data.frame(t(col2),stringsAsFactors = F)
rownames(col2) <- indiv
colnames(col2) <- SNPs

#check amount of missing data per individual
col2 <- data.frame(lapply(col2, gsub, pattern = "./.", replacement = NA, fixed = TRUE))
col2 <- data.frame(lapply(col2, gsub, pattern = "0/0", replacement = 0, fixed = TRUE))
col2 <- data.frame(lapply(col2, gsub, pattern = "0/1", replacement = 1, fixed = TRUE))
col2 <- data.frame(lapply(col2, gsub, pattern = "1/1", replacement = 2, fixed = TRUE))
save(col2,file = "SNPsets/gentype.for.PCA.che.rda")
#convert to genlight object for PCA
snp <- new("genlight",
           col2,
           chromosome=loci$scaffold,
           position=loci$pos,
           pop=locs$spp)
indNames(snp) <- indiv
save(snp,file = "genlight_pca_che.RData")
#running PCA analysis
#run pca
pca <- glPca(snp,nf = 10)
save(pca,file = "SNPsets/pca_che.RData")

#getting colors
locs <- data.frame(spp=snp@pop,stringsAsFactors = F)
locs <- locs %>% 
  select(spp)
locs$col <- NA
for(i in 1:length(locs$spp)){
  if(locs$spp[i] == "LA"){locs$col[i] <- "blue"}
  if(locs$spp[i] == "IF"){locs$col[i] <- "dark green"}
  if(locs$spp[i] == "PM"){locs$col[i] <- "purple"}
}


#plot pca
tiff(filename = "Figures/CHE_PCA_plot.tiff",width = 5,height = 5,units = "in",res = 400)
plot(pca$scores[,1], pca$scores[,2],
     col=locs$col,cex=0.5)
text(x = 0, y = -40, "P.marinus", pos = 2, col = "purple",cex = 0.8)
text(x = -50, y = 25, "I.fosser", pos = 3, col = "dark green",cex = 0.8)
text(x = 100, y = 20, "L.appendix", pos = 2, col = "blue",cex = 0.8)
title("Cheboyan River")
dev.off()
#make a tree
#construct tree
tre <- nj(dist(as.matrix(snp)))
#plot tree
plot(tre, typ="fan", cex=0.7)
