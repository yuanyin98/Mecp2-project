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
notboxplot(data3, ylim=c(0,12), ylab="Mean burst duration per active cell per recording (sec)",group.names= c("WT","KO"),cex.lab = 2.1, cex.axis =2)



data<-read.csv(file.choose(),header=T)
head(data)
library(ggplot2)

#1-boxplot
#从下往上: “minimum”, first quartile (Q1=25%), median, third quartile (Q3=75%), and “maximum”
dur_pacpr<- ggplot(data,aes(x= factor(grp, level = c('WT','KO')),y=dv))+geom_boxplot(aes(color=grp),size=2,position=position_dodge(width=0.5), outlier.shape = NA)+
  geom_jitter(colour = "grey60",size=3) #to prevent points overlay
print(dur_pacpr + labs(y="Mean burst duration per active cell per recording (sec)", x = NULL,colour = "Genotype") +
        theme(
          axis.title.x = element_text(size = 35),
          axis.title.y = element_text(size = 25),
          legend.title = element_text(size=23),
          legend.text = element_text(size = 23),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20)))


data<-read.csv(file.choose(),header=T)
head(data)

#2-t-test
t.test(dv~grp,data=data,alternative="two.sided")
#p-value =  0.3346, not significant difference
# the smaller the p-value, the stronger the evidence is that the two populations have different means
wilcox.test(dv~grp,data=data,alternative="two.sided")
#p-value = 0.3144, not significant difference
# the p-value is the probability that observed difference in mean happens by chance


#change C & A name
length(data$dv)
C<-subset(data,data$grp=="KO")
length(C$dv)
sapply(C, na.rm=TRUE)
summary(C)
#dv         grp    
#Min.   : 6.693   KO:11  
#1st Qu.: 7.441   WT: 0  
#Median : 8.248          
#Mean   : 8.449          
#3rd Qu.: 9.612          
#Max.   :10.864 
A<-subset(data,data$grp=="WT")
length(A$dv)
sapply(A, na.rm=TRUE)
summary(A)
#Min.   : 6.197   KO: 0  
#1st Qu.: 6.709   WT:10  
#Median : 7.651          
#Mean   : 7.819          
#3rd Qu.: 8.357          
#Max.   :11.462 


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
  geom_histogram(color=hue_pal()(4)[1], fill= hue_pal()(4)[1], alpha=0.5, breaks = seq(6,12,0.5))
print(a + labs(y="Number of recordings", x = "Average burst duration/active cell/recording (sec)",fill = "Genotype", color = "Genotype") +
        theme(
          axis.title.x = element_text(size = 30),
          axis.title.y = element_text(size = 35),
          axis.text.x = element_text(size = 25),
          axis.text.y = element_text(size = 25)))
#histogram for C=KO
c<-ggplot(C,aes(x=dv)) +
  geom_histogram(color=hue_pal()(4)[3], fill= hue_pal()(4)[3], alpha=0.5, breaks = seq(6,11,0.5))
print(c + labs(y="Number of recordings", x = "Average burst duration/active cell/recording (sec)",fill = "Genotype", color = "Genotype") +
        theme(
          axis.title.x = element_text(size = 35),
          axis.title.y = element_text(size = 35),
          axis.text.x = element_text(size = 25),
          axis.text.y = element_text(size = 25)))



#4-correlation between Average burst duration per active cell per recording (sec) VS Average number of bursts/recording
# Enhanced Scatterplot of `Average burst duration per active cell per recording (sec)` on the y-axis vs. `Average number of bursts/recording` on the x-axis
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

basic <- ggplot(data4, aes(ave_No_bursts, ave_burst_dur, colour = Genotype, shape = Genotype)) +
  geom_point(size=5)+
  geom_smooth(method=lm, se=FALSE)+
  stat_cor(method = "pearson", size = 5, label.x = 5, show.legend = FALSE)

print(basic+ labs(x="Mean number of bursts per active cell per recording",y="Mean burst duration per active cell per recording (sec)",color="Genotype")+
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
    ave_burst_dur ~ ave_No_bursts + Genotype + ave_No_bursts*Genotype
  )

# Normality of Residuals
# Fit the model, the covariate goes first
model <- lm(ave_burst_dur ~ ave_No_bursts + Genotype + ave_No_bursts*Genotype, data = data4)
# Inspect the model diagnostic metrics
model.metrics <- augment(model) %>%
  select(-.hat, -.sigma, -.fitted, -.se.fit) # Remove details
head(model.metrics, 3)

# when using ANCOVA, need to check the variance of the residuals is equal for all groups by Levene’s test
levene_test(.resid ~ Genotype, data = model.metrics)
# p= 0.796, p>0.05, The Levene’s test was not significant (p > 0.05), so we can assume homogeneity of the residual variances for all groups.

# ANCOVA2
# Or use 
Anova(lm(ave_burst_dur ~ ave_No_bursts * Genotype, data=data4, type=2))
# Anova Table (Type II tests)
# Response: ave_burst_dur
#                       Sum Sq Df F value         Pr(>F)
#ave_No_bursts           6.471  1  3.2959         0.08713 .
#Genotype                2.215  1  1.1281         0.30304  
#ave_No_bursts:Genotype  0.024  1  0.0122         0.91341  
#Residuals              33.375 17                       
#---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1




#5-correlation between Average burst duration per active cell per recording (sec) VS Number of active cells/recording
# Enhanced Scatterplot of `Average burst duration per active cell per recording (sec)` on the y-axis vs. `Number of active cells/recording` on the x-axis
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

basic <- ggplot(data5, aes(No_active_cells, ave_burst_dur, colour = Genotype, shape = Genotype)) +
  geom_point(size=5)+
  geom_smooth(method=lm, se=FALSE)+
  stat_cor(method = "pearson", size = 5, label.x = 15, show.legend = FALSE)

print(basic+ labs(x="Number of active cells/recording",y="Average burst duration/active cell/recording (sec)",color="Genotype")+
        theme(
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 17),
          legend.title = element_text(size=15),
          legend.text = element_text(size = 15)))

#-1 indicates a strong negative correlation : this means that every time x increases, y decreases (left panel figure)
#0 means that there is no association between the two variables (x and y) (middle panel figure)
#1 indicates a strong positive correlation : this means that y increases with x (right panel figure)



