data2<-read.csv(file.choose(),header=T)
head(data2)
data3 = list(as.double(data2$WT[!is.na(data2$WT)]),as.double(data2$KO[!is.na(data2$KO)]))
data3
library(ggplot2)

#1-notBoxPlot: mean, 95% confidence interval & 1SD
#function
notboxplot <- function(x, ...){
  means <- sapply(x, mean)
  sds <- sapply(x, sd)
  cis.d <- sapply(x,function(x) mean(x)-1.96*sd(x)/sqrt(length(x)))
  cis.u <- sapply(x,function(x) mean(x)+1.96*sd(x)/sqrt(length(x)))
  n<- length(x)
  a <- 1:n - 0.25
  b <- 1:n + 0.25
  stripchart(x, xlim=c(0.5,n+.5), vertical = TRUE, col="white",...)
  rect(a,means-sds,b,means+sds,col="mediumpurple2",border=NA)
  rect(a,cis.d,b,cis.u,col="pink",border=NA)
  for(i in 1:n) lines(c(i-.25,i+.25),c(means[i],means[i]),lwd=8,col='red')
  stripchart(x,method = 'jitter',add=TRUE,pch=16,vertical = TRUE, col='grey40',cex=1,lwd=2)
}
#notBoxplot
par(mar=c(3,5,3,5)+.1)
notboxplot(data3,ylab="Mean IEI per active cell (sec)",group.names= c("WT","KO"),cex.lab = 2.2, cex.axis = 2)






data<-read.csv(file.choose(),header=T)
head(data)

library(ggplot2)
#1-boxplot
#从下往上: “minimum”, first quartile (Q1=25%), median, third quartile (Q3=75%), and “maximum”
IEI_pac<- ggplot(data,aes(x= factor(grp, level = c('WT','KO')),y=dv))+geom_boxplot(aes(color=grp),size=2,position=position_dodge(width=0.5), outlier.shape = NA)+
  geom_jitter(colour = "grey60",size=1) #to prevent points overlay
print(IEI_pac + labs(y="Mean IEI/active cell (sec)", x = NULL,colour = "Genotype") +
        theme(
          axis.title.x = element_text(size = 35),
          axis.title.y = element_text(size = 35),
          legend.title = element_text(size=25),
          legend.text = element_text(size = 25),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20)))



#2-t-test
t.test(dv~grp,data=data,alternative="two.sided")
#p-value =  0.9155, not significant difference
# the smaller the p-value, the stronger the evidence is that the two populations have different means
wilcox.test(dv~grp,data=data,alternative="two.sided")
#p-value = 0.4515, not significant difference
# the p-value is the probability that observed difference in mean happens by chance


#change C & A name
length(data$dv)
C<-subset(data,data$grp=="KO")
head(C)
length(C$dv)
sapply(C, na.rm=TRUE)
summary(C)
A<-subset(data,data$grp=="WT")
length(A$dv)
sapply(A, na.rm=TRUE)
summary(A)

#2-sample K-S test
ks.test(A$dv,C$dv,data=data,alternative="two.sided")
#p-value = 0.9325, not significant difference
#You reject the null hypothesis that the two samples were drawn from the same distribution if the p-value is less than your significance level.

#3-normal distribution
#Normal Q-Q plot for A
qqnorm(A$dv,main = 'Normal Q-Q Plot for WT',cex.lab = 1.4)
qqline(A$dv)
#Normal Q-Q plot for C
qqnorm(C$dv,main = 'Normal Q-Q Plot for KO',cex.lab = 1.4)
qqline(C$dv)
#histogram for A=WT
par(mar=c(5,1,4,3)+.1)
a<-ggplot(A,aes(x=dv)) +
  geom_histogram(color=hue_pal()(4)[1], fill= hue_pal()(4)[1], alpha=0.5, breaks = seq(0,500,5))
print(a + labs(y="Number of active cells", x = "Mean IEI/active cell (sec)",fill = "Genotype", color = "Genotype") +
        theme(
          axis.title.x = element_text(size = 35),
          axis.title.y = element_text(size = 35),
          axis.text.x = element_text(size = 25),
          axis.text.y = element_text(size = 25)))
#histogram for C=KO
c<-ggplot(C,aes(x=dv)) +
  geom_histogram(color=hue_pal()(4)[3], fill= hue_pal()(4)[3], alpha=0.5, breaks = seq(0,500,5))
print(c + labs(y="Number of active cells", x = "Mean IEI/active cell (sec)",fill = "Genotype", color = "Genotype") +
        theme(
          axis.title.x = element_text(size = 35),
          axis.title.y = element_text(size = 35),
          axis.text.x = element_text(size = 25),
          axis.text.y = element_text(size = 25)))

#4-correlation between Mean IEI/active cell (sec) VS Number of bursts/active cells
# Enhanced Scatterplot of `AMean IEI/active cell (sec)` on the y-axis vs. `Number of bursts/active cells` on the x-axis
# by Genotype
data4<-read.csv(file.choose(),header=T)
head(data4)
library(car)
library(reshape2)
library(ggpubr)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(broom)

basic <- ggplot(data4, aes(No_bursts, mean_IEI, colour = Genotype, shape = Genotype)) +
  geom_point(size=1.3)+
  geom_smooth(method=lm, se=FALSE,size=1.5)+
  stat_cor(method = "pearson", size = 5, label.x = 10, show.legend = FALSE)

print(basic+ labs(x="Number of bursts per active cell",y="Mean IEI per active cell (sec)",color="Genotype")+
        theme(
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          legend.title = element_text(size=15),
          legend.text = element_text(size = 15), 
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15)))

#-1 indicates a strong negative correlation : this means that every time x increases, y decreases (left panel figure)
#0 means that there is no association between the two variables (x and y) (middle panel figure)
#1 indicates a strong positive correlation : this means that y increases with x (right panel figure)

library(car)
library(reshape2)
library(ggpubr)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(broom)
# ANCOVA1
data4 %>%
  anova_test(
    mean_IEI ~ No_bursts + Genotype + No_bursts*Genotype
  )

# Normality of Residuals
# Fit the model, the covariate goes first
model <- lm(mean_IEI ~ No_bursts + Genotype + No_bursts*Genotype, data = data4)
# Inspect the model diagnostic metrics
model.metrics <- augment(model) %>%
  select(-.hat, -.sigma, -.fitted, -.se.fit) # Remove details
head(model.metrics, 3)

# when using ANCOVA, need to check the variance of the residuals is equal for all groups by Levene’s test
levene_test(.resid ~ Genotype, data = model.metrics)
# p=0.625, p>0.05, The Levene’s test was not significant (p > 0.05), so we can assume homogeneity of the residual variances for all groups.

# ANCOVA2
# Or use 
Anova(lm(mean_IEI ~ No_bursts * Genotype, data=data4, type=2))
# Anova Table (Type II tests)
# Response: mean_IEI
#                     Sum Sq  Df  F value           Pr(>F)    
# No_bursts           439199   1 110.2261           <2e-16 ***
# Genotype              1461   1   0.3668           0.5452    
# No_bursts:Genotype     115   1   0.0289           0.8651    
# Residuals          1227235 308                        
#---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

