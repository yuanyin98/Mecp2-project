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
  stripchart(x,method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg='gray70',cex=1.5,lwd=2)
}
#notBoxplot
par(mar=c(3,5,3,3)+.1)
notboxplot(data3, ylim=c(0,250), ylab="Mean IEI per active cell per recording (sec)",group.names= c("WT","KO"),cex.lab = 2.1, cex.axis =2)






data<-read.csv(file.choose(),header=T)
head(data)
library(ggplot2)

#1-boxplot
#从下往上: “minimum”, first quartile (Q1=25%), median, third quartile (Q3=75%), and “maximum”
IEI_pacpr<- ggplot(data,aes(x= factor(grp, level = c('WT','KO')),y=dv))+geom_boxplot(aes(color=grp),size=2,position=position_dodge(width=0.5), outlier.shape = NA)+
  geom_jitter(colour = "grey60",size=3) #to prevent points overlay
print(IEI_pacpr + labs(y="Mean IEI per active cell per recording (sec)", x = NULL,colour = "Genotype") +
        theme(
          axis.title.x = element_text(size = 35),
          axis.title.y = element_text(size = 25),
          legend.title = element_text(size=23),
          legend.text = element_text(size = 23),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20)))



#2-t-test
t.test(dv~grp,data=data,alternative="two.sided")
#p-value = 0.6808, not significant difference
# the smaller the p-value, the stronger the evidence is that the two populations have different means
wilcox.test(dv~grp,data=data,alternative="two.sided")
#p-value = 0.7045, not significant difference
# the p-value is the probability that observed difference in mean happens by chance


#change C & A name
length(data$dv)
C<-subset(data,data$grp=="KO")
length(C$dv)
sapply(C, na.rm=TRUE)
summary(C)
#dv         grp    
#Min.   : 27.84   KO:11  
#1st Qu.: 46.94   WT: 0  
#Median : 69.30          
#Mean   : 85.20          
#3rd Qu.:108.24          
#Max.   :215.86  
A<-subset(data,data$grp=="WT")
length(A$dv)
sapply(A, na.rm=TRUE)
summary(A)
#dv         grp    
#Min.   : 26.26   KO: 0  
#1st Qu.: 44.64   WT:10  
#Median : 58.81          
#Mean   : 75.89          
#3rd Qu.:112.03          
#Max.   :153.38   


#3-normal distribution
#Normal Q-Q plot for A
qqnorm(A$dv,main = 'Normal Q-Q Plot for WT',cex.lab = 1.4)
qqline(A$dv)
#Normal Q-Q plot for C
qqnorm(C$dv,main = 'Normal Q-Q Plot for KO',cex.lab = 1.4)
qqline(C$dv)
#histogram for A=WT
par(mar=c(5,3,4,3)+.1)
a<-ggplot(A,aes(x=dv)) +
  geom_histogram(color=hue_pal()(4)[1], fill= hue_pal()(4)[1], alpha=0.5, breaks = seq(20,160,10))
print(a + labs(y="Number of recordings", x = "Average IEI/active cell/recording (sec)",fill = "Genotype", color = "Genotype") +
        theme(
          axis.title.x = element_text(size = 30),
          axis.title.y = element_text(size = 35),
          axis.text.x = element_text(size = 25),
          axis.text.y = element_text(size = 25)))
#histogram for C=KO
c<-ggplot(C,aes(x=dv)) +
  geom_histogram(color=hue_pal()(4)[3], fill= hue_pal()(4)[3], alpha=0.5, breaks = seq(0,250,10))
print(c + labs(y="Number of recordings", x = "Average IEI/active cell/recording (sec)",fill = "Genotype", color = "Genotype") +
        theme(
          axis.title.x = element_text(size = 35),
          axis.title.y = element_text(size = 35),
          axis.text.x = element_text(size = 25),
          axis.text.y = element_text(size = 25)))

#4-correlation between Average IEI/active cell/recording (sec) VS Average number of bursts/recording
# Enhanced Scatterplot of `Average IEI/active cell/recording (sec)` on the y-axis vs. `Average number of bursts/recording` on the x-axis
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

basic <- ggplot(data4, aes(ave_No_bursts, ave_IEI, colour = Genotype, shape = Genotype)) +
  geom_point(size=5)+
  geom_smooth(method=lm, se=FALSE)+
  stat_cor(method = "pearson", size = 5, label.x = 5, show.legend = FALSE)

print(basic+ labs(x="Mean number of bursts per active cell per recording",y="Mean IEI per active cell per recording (sec)",color="Genotype")+
        theme(
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          legend.title = element_text(size=15),
          legend.text = element_text(size = 15),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15)))

# ANCOVA1
data4 %>%
  anova_test(
    ave_IEI ~ ave_No_bursts + Genotype + ave_No_bursts*Genotype
  )

# Normality of Residuals
# Fit the model, the covariate goes first
model <- lm(ave_IEI ~ ave_No_bursts + Genotype + ave_No_bursts*Genotype, data = data4)
# Inspect the model diagnostic metrics
model.metrics <- augment(model) %>%
  select(-.hat, -.sigma, -.fitted, -.se.fit) # Remove details
head(model.metrics, 3)

# when using ANCOVA, need to check the variance of the residuals is equal for all groups by Levene’s test
levene_test(.resid ~ Genotype, data = model.metrics)
# p=0.269, p>0.05, The Levene’s test was not significant (p > 0.05), so we can assume homogeneity of the residual variances for all groups.

# ANCOVA2
# Or use 
Anova(lm(ave_IEI ~ ave_No_bursts * Genotype, data=data4, type=2))
# Anova Table (Type II tests)
# Response: ave_IEI
#                         Sum Sq Df F value          Pr(>F)    
# ave_No_bursts           34647  1 37.2487       1.171e-05 ***
# Genotype                  610  1  0.6556          0.4293    
# ave_No_bursts:Genotype     31  1  0.0334          0.8572    
# Residuals               15813 17
#---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



#5-correlation between Average IEI/active cell/recording (sec) VS Number of active cells/recording
# Enhanced Scatterplot of `Average IEI/active cell/recording (sec)` on the y-axis vs. `Number of active cells/recording` on the x-axis
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

basic <- ggplot(data5, aes(No_active_cells, ave_IEI, colour = Genotype, shape = Genotype)) +
  geom_point(size=5)+
  geom_smooth(method=lm, se=FALSE)+
  stat_cor(method = "pearson", size = 5, label.x = 15, show.legend = FALSE)

print(basic+ labs(x="Number of active cells/recording",y="Average IEI/active cell/recording (sec)",color="Genotype")+
        theme(
          axis.title.x = element_text(size = 22),
          axis.title.y = element_text(size = 22),
          legend.title = element_text(size=20),
          legend.text = element_text(size = 20),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15)))

#-1 indicates a strong negative correlation : this means that every time x increases, y decreases (left panel figure)
#0 means that there is no association between the two variables (x and y) (middle panel figure)
#1 indicates a strong positive correlation : this means that y increases with x (right panel figure)


# ANCOVA1
data5 %>%
  anova_test(
    ave_IEI ~ No_active_cells + Genotype + No_active_cells*Genotype
  )

# Normality of Residuals
# Fit the model, the covariate goes first
model <- lm(ave_IEI ~ No_active_cells + Genotype + No_active_cells*Genotype, data = data5)
# Inspect the model diagnostic metrics
model.metrics <- augment(model) %>%
  select(-.hat, -.sigma, -.fitted, -.se.fit) # Remove details
head(model.metrics, 3)

# when using ANCOVA, need to check the variance of the residuals is equal for all groups by Levene’s test
levene_test(.resid ~ Genotype, data = model.metrics)
# p=0.938, p>0.05, The Levene’s test was not significant (p > 0.05), so we can assume homogeneity of the residual variances for all groups.

# ANCOVA2
# Or use 
Anova(lm(ave_IEI ~ No_active_cells * Genotype, data=data5, type=2))
# Anova Table (Type II tests)
# Response: ave_IEI
#                           Sum Sq Df F value        Pr(>F)  
# No_active_cells           10030  1  4.5987       0.04674 *
# Genotype                    199  1  0.0910       0.76653  
# No_active_cells:Genotype   3381  1  1.5499       0.23003  
# Residuals                 37080 17 
#---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

