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
par(mar=c(3,5,4,5)+.1)
notboxplot(data3,ylab="Number of bursts per active cell",group.names= c("WT","KO"),cex.lab = 2.2, cex.axis=2)





data<-read.csv(file.choose(),header=T)
head(data)

#2-t-test
t.test(dv~grp,data=data,alternative="two.sided")
#p-value =  0.0216, significant difference
# the smaller the p-value, the stronger the evidence is that the two populations have different means
wilcox.test(dv~grp,data=data,alternative="two.sided")
#p-value = 0.006677, significant difference
# the p-value is the probability that observed difference in mean happens by chance
# Both the Mann-Whitney and the Kolmogorov-Smirnov tests are nonparametric tests to compare two unpaired groups of data. Both compute P values that test the null hypothesis that the two groups have the same distribution. But they work very differently:
# •The Mann-Whitney test first ranks all the values from low to high, and then computes a P value that depends on the discrepancy between the mean ranks of the two groups.
# •The Kolmogorov-Smirnov test compares the cumulative distribution of the two data sets, and computes a P value that depends on the largest discrepancy between distributions.
# Here are some guidelines for choosing between the two tests:
# •The KS test is sensitive to any differences in the two distributions. Substantial differences in shape, spread or median will result in a small P value. In contrast, the MW test is mostly sensitive to changes in the median.
# •The MW test is used more often and is recognized by more people, so choose it if you have no idea which to choose.
# •The MW test has been extended to handle tied values. The KS test does not handle ties so well. If your data are categorical, so has many ties, don't choose the KS test.
# •Some fields of science tend to prefer the KS test over the MW test. It makes sense to follow the traditions of your field.


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
#p-value =  0.1401, not significant difference
# You reject the null hypothesis that the two samples were drawn from the same distribution if the p-value is less than your significance level.




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
  geom_histogram(color=hue_pal()(4)[1], fill= hue_pal()(4)[1], alpha=0.5, breaks = seq(0,30,5))
print(a + labs(y="Number of active cells", x = "Number of bursts/active cells",fill = "Genotype", color = "Genotype") +
        theme(
          axis.title.x = element_text(size = 35),
          axis.title.y = element_text(size = 35),
          axis.text.x = element_text(size = 25),
          axis.text.y = element_text(size = 25)))
#histogram for C=KO
c<-ggplot(C,aes(x=dv)) +
  geom_histogram(color=hue_pal()(4)[3], fill= hue_pal()(4)[3], alpha=0.5, breaks = seq(0,30,5))
print(c + labs(y="Number of active cells", x = "Number of bursts/active cells",fill = "Genotype", color = "Genotype") +
        theme(
          axis.title.x = element_text(size = 35),
          axis.title.y = element_text(size = 35),
          axis.text.x = element_text(size = 25),
          axis.text.y = element_text(size = 25)))



#4-correlation between Number of bursts/active cell VS Number of active cells
# Enhanced Scatterplot of `Number of bursts/active cell` on the y-axis vs. `Number of active cells` on the x-axis
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

basic <- ggplot(data4, aes(Active_cells, No_bursts, colour = Genotype, shape = Genotype)) +
  geom_point(size=1.3)+
  geom_smooth(method=lm, se=FALSE,size=1.5)+
  stat_cor(method = "pearson", size = 5, label.x = 5, show.legend = FALSE)

print(basic+ labs(x="Number of active cells per recording",y="Number of bursts per active cell",color="Genotype")+
        theme(
          axis.title.x = element_text(size = 25),
          axis.title.y = element_text(size = 25),
          legend.title = element_text(size=20),
          legend.text = element_text(size = 20), 
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20)))
# ANCOVA1
data4 %>%
  anova_test(
    No_bursts ~ Active_cells + Genotype + Active_cells*Genotype
  )

# Normality of Residuals
# Fit the model, the covariate goes first
model <- lm(No_bursts ~ Active_cells + Genotype + Active_cells*Genotype, data = data4)
# Inspect the model diagnostic metrics
model.metrics <- augment(model) %>%
  select(-.hat, -.sigma, -.fitted, -.se.fit) # Remove details
head(model.metrics, 3)

# when using ANCOVA, need to check the variance of the residuals is equal for all groups by Levene’s test
levene_test(.resid ~ Genotype, data = model.metrics)
# p=0.668, p>0.05, The Levene’s test was not significant (p > 0.05), so we can assume homogeneity of the residual variances for all groups.

# ANCOVA2
# Or use 
Anova(lm(No_bursts ~ Active_cells * Genotype, data=data4, type=2))
# Anova Table (Type II tests)
# Response: No_bursts
#                        Sum Sq  Df F value       Pr(>F)    
# Active_cells           1317.9   1  34.039    1.124e-08 ***
# Genotype                604.1   1  15.602    9.253e-05 ***
# Active_cells:Genotype   571.0   1  14.747    0.0001431 ***
# Residuals             15332.3 396                      
#---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



#5-correlation between Number of bursts/active cell VS Percentage of active cells
# Enhanced Scatterplot of `Number of bursts/active cell` on the y-axis vs. `Percentage of active cells` on the x-axis
# by Genotype
data5<-read.csv(file.choose(),header=T)
head(data5)
library(car)
library(reshape2)
library(ggpubr)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(broom)

basic <- ggplot(data5, aes(percentage_Active_cells, No_bursts, colour = Genotype, shape = Genotype)) +
  geom_point(size=1.3)+
  geom_smooth(method=lm, se=FALSE,size=1.5)+
  stat_cor(method = "pearson", size = 5, label.x = 5, show.legend = FALSE)

print(basic+ labs(x="Percent active cells per recording (%)",y="Number of bursts per active cell",color="Genotype")+
        theme(
          axis.title.x = element_text(size = 25),
          axis.title.y = element_text(size = 25),
          legend.title = element_text(size=20),
          legend.text = element_text(size = 20), 
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20)))

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
data5 %>%
  anova_test(
    No_bursts ~ percentage_Active_cells + Genotype + percentage_Active_cells*Genotype
  )

# Normality of Residuals
# Fit the model, the covariate goes first
model <- lm(No_bursts ~ percentage_Active_cells + Genotype + percentage_Active_cells*Genotype, data = data5)
# Inspect the model diagnostic metrics
model.metrics <- augment(model) %>%
  select(-.hat, -.sigma, -.fitted, -.se.fit) # Remove details
head(model.metrics, 3)

# when using ANCOVA, need to check the variance of the residuals is equal for all groups by Levene’s test
levene_test(.resid ~ Genotype, data = model.metrics)
# p=0.475, p>0.05, The Levene’s test was not significant (p > 0.05), so we can assume homogeneity of the residual variances for all groups.

# ANCOVA2
# Or use 
Anova(lm(No_bursts ~ percentage_Active_cells * Genotype, data=data5, type=2))
# Anova Table (Type II tests)
# Response: No_bursts
#                                   Sum Sq  Df F value              Pr(>F)    
# percentage_Active_cells           2591.9   1 70.1657           9.566e-16 ***
# Genotype                           562.3   1 15.2208           0.0001124 ***
# percentage_Active_cells:Genotype     1.1   1  0.0288           0.8652595    
#Residuals                        14628.1 396                     
#---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1





