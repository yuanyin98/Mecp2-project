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
notboxplot(data3,ylab="Mean burst duration per active cell (sec)",group.names= c("WT","KO"),cex.lab = 2.2, cex.axis = 2)





data<-read.csv(file.choose(),header=T)
head(data)

#2-t-test
t.test(dv~grp,data=data,alternative="two.sided")
#p-value =  0.1206, not significant difference
# the smaller the p-value, the stronger the evidence is that the two populations have different means
wilcox.test(dv~grp,data=data,alternative="two.sided")
#p-value =  0.01215, significant difference
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
  geom_histogram(color=hue_pal()(4)[1], fill= hue_pal()(4)[1], alpha=0.5, breaks = seq(0,20,1))
print(a + labs(y="Number of active cells", x = "Mean burst duration/active cell (sec)",fill = "Genotype", color = "Genotype") +
        theme(
          axis.title.x = element_text(size = 35),
          axis.title.y = element_text(size = 35),
          axis.text.x = element_text(size = 25),
          axis.text.y = element_text(size = 25)))
#histogram for C=KO
c<-ggplot(C,aes(x=dv)) +
  geom_histogram(color=hue_pal()(4)[3], fill= hue_pal()(4)[3], alpha=0.5, breaks = seq(0,20,1))
print(c + labs(y="Number of active cells", x = "Mean burst duration/active cell (sec)",fill = "Genotype", color = "Genotype") +
        theme(
          axis.title.x = element_text(size = 35),
          axis.title.y = element_text(size = 35),
          axis.text.x = element_text(size = 25),
          axis.text.y = element_text(size = 25)))

#4-correlation between Mean burst duration/active cell (sec) VS Number of bursts/active cells
# Enhanced Scatterplot of `AMean burst duration/active cell (sec)` on the y-axis vs. `Number of bursts/active cells` on the x-axis
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

basic <- ggplot(data4, aes(No_bursts, mean_burst_dur, colour = Genotype, shape = Genotype)) +
  geom_point(size=1.5)+
  geom_smooth(method=lm, se=FALSE)+
  stat_cor(method = "pearson", size = 5, label.x = 5, show.legend = FALSE)

print(basic+ labs(x="Number of bursts per active cell",y="Mean burst duration per active cell (sec)",color="Genotype")+
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
    mean_burst_dur ~ No_bursts + Genotype + No_bursts*Genotype
  )

# Normality of Residuals
# Fit the model, the covariate goes first
model <- lm(mean_burst_dur ~ No_bursts + Genotype + No_bursts*Genotype, data = data4)
# Inspect the model diagnostic metrics
model.metrics <- augment(model) %>%
  select(-.hat, -.sigma, -.fitted, -.se.fit) # Remove details
head(model.metrics, 3)

# when using ANCOVA, need to check the variance of the residuals is equal for all groups by Levene’s test
levene_test(.resid ~ Genotype, data = model.metrics)
# p=0.00114, p<0.05, The Levene’s test was significant (p < 0.05), so we can NOT assume homogeneity of the residual variances for all groups.

# ANCOVA2
# Or use 
Anova(lm(mean_burst_dur ~ No_bursts * Genotype, data=data4, type=2))
# Anova Table (Type II tests)
# Response: mean_burst_dur
#                     Sum Sq  Df F value            Pr(>F)  
# No_bursts            22.3   1  2.6780           0.10254  
# Genotype             23.2   1  2.7880           0.09576 .
# No_bursts:Genotype    3.7   1  0.4409           0.50708  
# Residuals          3294.2 396                       
#---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

