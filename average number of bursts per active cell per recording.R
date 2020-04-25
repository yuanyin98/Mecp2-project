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
  stripchart(x,method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg='gray70',cex=2,lwd=2)
}
#notBoxplot
par(mar=c(2,5,1,3)+.1)
notboxplot(data3, ylim=c(0,12), ylab="Mean number of bursts per active cell per recording",group.names= c("WT","KO"),cex.lab = 2.1, cex.axis =2)





data<-read.csv(file.choose(),header=T)
head(data)

#2-t-test
t.test(dv~grp,data=data,alternative="two.sided")
#p-value =  0.937, not significant difference
# the smaller the p-value, the stronger the evidence is that the two populations have different means
wilcox.test(dv~grp,data=data,alternative="two.sided")
#p-value = 1, not significant difference
# the p-value is the probability that observed difference in mean happens by chance


#change C & A name
length(data$dv)
C<-subset(data,data$grp=="KO")
length(C$dv)
sapply(C, na.rm=TRUE)
summary(C)
#dv         grp    
#Min.   : 1.600   KO:11  
#1st Qu.: 2.310   WT: 0  
#Median : 4.882          
#Mean   : 5.341          
#3rd Qu.: 7.375          
#Max.   :11.714   
A<-subset(data,data$grp=="WT")
length(A$dv)
sapply(A, na.rm=TRUE)
summary(A)
#dv         grp    
#Min.   : 1.349   KO: 0  
#1st Qu.: 3.219   WT:10  
#Median : 5.195          
#Mean   : 5.229          
#3rd Qu.: 7.103          
#Max.   :10.267   


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
  geom_histogram(color=hue_pal()(4)[1], fill= hue_pal()(4)[1], alpha=0.5, breaks = seq(0,12,1))
print(a + labs(y="Number of recordings", x = "Average number of bursts/active cell/recording",fill = "Genotype", color = "Genotype") +
        theme(
          axis.title.x = element_text(size = 30),
          axis.title.y = element_text(size = 35),
          axis.text.x = element_text(size = 25),
          axis.text.y = element_text(size = 25)))
#histogram for C=KO
c<-ggplot(C,aes(x=dv)) +
  geom_histogram(color=hue_pal()(4)[3], fill= hue_pal()(4)[3], alpha=0.5, breaks = seq(0,12,2))
print(c + labs(y="Number of recordings", x = "Average number of bursts/active cell/recording",fill = "Genotype", color = "Genotype") +
        theme(
          axis.title.x = element_text(size = 35),
          axis.title.y = element_text(size = 35),
          axis.text.x = element_text(size = 25),
          axis.text.y = element_text(size = 25)))




#4-correlation between Average number of bursts/active cell/recording VS Number of cells 
# Enhanced Scatterplot of `Average number of bursts/active cell/recording ` on the y-axis vs. `Number of cells/recording` on the x-axis
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

basic <- ggplot(data4, aes(Cells, ave_No_bursts, colour = Genotype, shape = Genotype)) +
  geom_point(size=5)+
  geom_smooth(method=lm, se=FALSE)+
  stat_cor(method = "pearson", size = 5, show.legend = FALSE, label.x = 10)

print(basic+ labs(x="Number of cells/recording",y="Average number of bursts/active cell/recording",color="Genotype")+
        theme(
          axis.title.x = element_text(size = 28),
          axis.title.y = element_text(size = 30),
          legend.title = element_text(size=20),
          legend.text = element_text(size = 20),
          axis.text.x = element_text(size = 25),
          axis.text.y = element_text(size = 25)))
# ANCOVA1
data4 %>%
  anova_test(
    ave_No_bursts ~ Cells + Genotype + Cells*Genotype
  )

# Normality of Residuals
# Fit the model, the covariate goes first
model <- lm(ave_No_bursts ~ Cells + Genotype + Cells*Genotype, data = data4)
# Inspect the model diagnostic metrics
model.metrics <- augment(model) %>%
  select(-.hat, -.sigma, -.fitted, -.se.fit) # Remove details
head(model.metrics, 3)

# when using ANCOVA, need to check the variance of the residuals is equal for all groups by Levene’s test
levene_test(.resid ~ Genotype, data = model.metrics)
# p=0.971, p>0.05, The Levene’s test was not significant (p > 0.05), so we can assume homogeneity of the residual variances for all groups.

# ANCOVA2
# Or use 
Anova(lm(ave_No_bursts ~ Cells * Genotype, data=data4, type=2))
# Anova Table (Type II tests)
# Response: ave_No_bursts
#                 Sum Sq Df F value         Pr(>F)
# Cells            1.064  1  0.0924         0.7648
# Genotype         0.090  1  0.0078         0.9307
# Cells:Genotype   0.764  1  0.0663         0.7999
# Residuals      195.838 17                        
#---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



#5-correlation between Average number of bursts/active cell/recording VS number of active cells
# Enhanced Scatterplot of `Average number of bursts/active cell/recording` on the y-axis vs. `number of active cells/recording` on the x-axis
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

basic <- ggplot(data5, aes(Active_cells,ave_No_bursts, colour = Genotype, shape = Genotype)) +
  geom_point(size=5)+
  geom_smooth(method=lm, se=FALSE)+
  stat_cor(method = "pearson", size = 5, show.legend = FALSE)

print(basic+ labs(x="Number of active cells per recording",y="Mean number of bursts per active cell per recording",color="Genotype")+
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


# ANCOVA1
data5 %>%
  anova_test(
    ave_No_bursts ~ Active_cells + Genotype + Active_cells*Genotype
  )

# Normality of Residuals
# Fit the model, the covariate goes first
model <- lm(ave_No_bursts ~ Active_cells + Genotype + Active_cells*Genotype, data = data5)
# Inspect the model diagnostic metrics
model.metrics <- augment(model) %>%
  select(-.hat, -.sigma, -.fitted, -.se.fit) # Remove details
head(model.metrics, 3)

# when using ANCOVA, need to check the variance of the residuals is equal for all groups by Levene’s test
levene_test(.resid ~ Genotype, data = model.metrics)
# p=0.196, p>0.05, The Levene’s test was not significant (p > 0.05), so we can assume homogeneity of the residual variances for all groups.

# ANCOVA2
# Or use 
Anova(lm(ave_No_bursts ~ Active_cells * Genotype, data=data5, type=2))
# Anova Table (Type II tests)
# Response: ave_No_bursts
#                       Sum Sq Df F value              Pr(>F)    
# Active_cells          90.942  1 17.5195           0.0006206 ***
# Genotype              12.180  1  2.3464           0.1439683    
# Active_cells:Genotype 18.478  1  3.5597           0.0763997 .  
# Residuals             88.246 17                          14628.1 396                     
#---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



#6-correlation between Average number of bursts/active cell/recording VS percentage of active cells
# Enhanced Scatterplot of `Average number of bursts/active cell/recording` on the y-axis vs. `percentage of active cells/recording` on the x-axis
# by Genotype
data6<-read.csv(file.choose(),header=T)
head(data6)
library(car)
library(reshape2)
library(ggpubr)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(broom)

basic <- ggplot(data6, aes(p_Active_cells,ave_No_bursts, colour = Genotype, shape = Genotype)) +
  geom_point(size=5)+
  geom_smooth(method=lm, se=FALSE)+
  stat_cor(method = "pearson", size = 5, show.legend = FALSE)

print(basic+ labs(x="Percent active cells per recording (%)",y="Mean number of bursts per active cell per recording",color="Genotype")+
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


# ANCOVA1
data6 %>%
  anova_test(
    ave_No_bursts ~ p_Active_cells + Genotype + p_Active_cells*Genotype
  )

# Normality of Residuals
# Fit the model, the covariate goes first
model <- lm(ave_No_bursts ~ p_Active_cells + Genotype + p_Active_cells*Genotype, data = data6)
# Inspect the model diagnostic metrics
model.metrics <- augment(model) %>%
  select(-.hat, -.sigma, -.fitted, -.se.fit) # Remove details
head(model.metrics, 3)

# when using ANCOVA, need to check the variance of the residuals is equal for all groups by Levene’s test
levene_test(.resid ~ Genotype, data = model.metrics)
# p=0.507, p>0.05, The Levene’s test was not significant (p > 0.05), so we can assume homogeneity of the residual variances for all groups.

# ANCOVA2
# Or use 
Anova(lm(ave_No_bursts ~ p_Active_cells * Genotype, data=data6, type=2))
# Anova Table (Type II tests)
# Response: ave_No_bursts
#                         Sum Sq Df F value           Pr(>F)    
# p_Active_cells          126.061  1 30.0921         4.02e-05 ***
# Genotype                 21.945  1  5.2386          0.03516 *  
# p_Active_cells:Genotype   0.389  1  0.0927          0.76442    
# Residuals                71.216 17                    
#---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
