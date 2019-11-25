library(car)

bmr2 <- bmr2 %>% 
  filter(reach != 4) %>% 
  mutate(family = paste(Mother,Father,sep = "_"))


bmr3 <- bmr2[bmr2$family %in% names(which(table(bmr2$family)>=5)),]
bmr3_2017 <- subset(bmr3,bmr3$Year_collect == 2017)
bmr3_2018 <- subset(bmr3,bmr3$Year_collect == 2018)

bmr3_2017_counts <- bmr3_2017 %>% group_by(family,reach) %>% count() %>% filter(n >= 3)
i <- 1
bmr4 <- data.frame(matrix(data  = NA,nrow = 1,ncol = 1:ncol(bmr3_2017)))
for (i in 1:length(bmr3_2017_counts$family)) {
  tmp <- bmr3_2017 %>% filter(family == bmr3_2017_counts$family[i])
  append(x = bmr4, values = tmp) 
}
summary(lm(data = bmr3,Length~as.factor(reach)))


an_model <- Anova(mod = lm(data = bmr3_2017,Length~as.factor(family)+as.factor(reach)),type = "III",singular.ok = T)

plot.design(data = bmr3_2018,Length~as.factor(family)+as.factor(reach))

boxplot(bmr3_2017$Length~bmr3_2017$family*bmr3_2017$reach,col = bmr3_2017$reach)

tiff(file="Output/bmr_2018_length_family.tiff",width=6, height=4, units="in", res=200)
boxplot(bmr3_2018$Length~bmr3_2018$family, las = 2,
        main = "Black Mallard length distribution by family \n2018 Collection",
        xlab = "family",
        ylab = "length (mm)",
        ylim = c(0,130),
        cex.axis = 0.6)
dev.off()

tiff(file="Output/bmr_2017_length_family.tiff",width=6, height=4, units="in", res=400)
boxplot(bmr3_2017$Length~bmr3_2017$family, las = 2,
        main = "Black Mallard length distribution by family \n2017 Collection",
        xlab = "family",
        ylab = "length (mm)",
        ylim = c(0,130),
        cex.axis = 0.6)
dev.off()

tiff(file="Output/bmr_length_year.tiff",width=6, height=4, units="in", res=200)
boxplot(bmr3$Length~bmr3$Year_collect, col = "darkblue",
        main = "Black Mallard length distributions \n  by collection year",
        xlab = "Year",
        ylab = "length (mm)",
        ylim = c(0,130),
        cex.axis = 0.5)
dev.off()
par(mfrow=c(3,1))
boxplot(bmr3_2018$Length~bmr3_2018$Father, col = "darkorange")
boxplot(bmr3_2017$Length~bmr3_2017$Father, col = "purple")
boxplot(bmr3$Length~bmr3$reach)


#Family analysis from the ocqueoc
ocq2 <- ocq2 %>% mutate(family = paste(Mother,Father,sep = "_"))
ocq3 <- ocq2 %>% filter(Mother == 1)
ocq4 <- ocq2 %>% filter(Mother != 1)
boxplot(ocq3$Length~ocq3$family, las = 2,
        main = "Ocqueoc length distribution for large half-sibling groups",
        xlab = "family",
        ylab = "length (mm)",
        ylim = c(60,160),
        cex.axis = 0.75)
boxplot(ocq4$Length~ocq4$family,las = 2,
        main = "Ocqueoc length distribution for outlier families",
        xlab = "family",
        ylab = "length (mm)",
        ylim = c(60,160))

