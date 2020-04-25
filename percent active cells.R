library(ggplot2)
data<-read.csv(file.choose(),header=T)
head(data)

#1-boxplot
#从下往上: “minimum”, first quartile (Q1=25%), median, third quartile (Q3=75%), and “maximum”
pac<- ggplot(data,aes(x= factor(grp, level = c('WT','KO')),y=dv))+geom_boxplot(aes(color=grp),size=2,position=position_dodge(width=0.5), outlier.shape = NA)+
  geom_jitter(colour = "grey60",size=3) #to prevent points overlay
print(pac + labs(y="Percent active cells per recording (%)", x = NULL,colour = "Genotype") +
        theme(
          axis.title.x = element_text(size = 30),
          axis.title.y = element_text(size = 30),
          legend.title = element_text(size=25),
          legend.text = element_text(size = 25),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20)))




data2<-read.csv(file.choose(),header=T)
head(data2)
data3 = list(as.double(data2$WT[!is.na(data2$WT)]),as.double(data2$KO[!is.na(data2$KO)]))
data3
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
par(mar=c(3,5,4,5)+.1)
notboxplot(data3,ylab="Percent active cells per recording (%)",group.names= c("WT","KO"),cex.lab = 2.1, cex.axis = 2)





data<-read.csv(file.choose(),header=T)
head(data)

#2-t-test
t.test(dv~grp,data=data,alternative="two.sided")
#p-value =  0.07193, not significant difference
# the smaller the p-value, the stronger the evidence is that the two populations have different means
wilcox.test(dv~grp,data=data,alternative="two.sided")
#p-value = 0.05134, not significant difference
# the p-value is the probability that observed difference in mean happens by chance


#change C & A name
length(data$dv)
C<-subset(data,data$grp=="KO")
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
  geom_histogram(color=hue_pal()(4)[1], fill= hue_pal()(4)[1], alpha=0.5, breaks = seq(0,70,5))
print(a + labs(y="Count", x = "Percentage of active cells/recording (%)",fill = "Genotype", color = "Genotype") +
        theme(
          axis.title.x = element_text(size = 25),
          axis.title.y = element_text(size = 35),
          axis.text.x = element_text(size = 25),
          axis.text.y = element_text(size = 25)))
#histogram for C=KO
c<-ggplot(C,aes(x=dv)) +
  geom_histogram(color=hue_pal()(4)[3], fill= hue_pal()(4)[3], alpha=0.5, breaks = seq(0,70,5))
print(c + labs(y="Count", x = "Percentage of active cells/recording (%)",fill = "Genotype", color = "Genotype") +
        theme(
          axis.title.x = element_text(size = 25),
          axis.title.y = element_text(size = 35),
          axis.text.x = element_text(size = 25),
          axis.text.y = element_text(size = 25)))


#4-correlation between Percentage of active cells VS Percentage of cells 
# Enhanced Scatterplot of `Percentage of active cells` on the y-axis vs. `Percentage of cells` on the x-axis
# by Genotype
data4<-read.csv(file.choose(),header=T)
head(data4)
library(ggplot2)
library(car)
library(reshape2)
library(ggpubr)

basic <- ggplot(data4, aes(Cells, Active_cells, colour = Genotype, shape = Genotype)) +
  geom_point(size=5)+
  geom_smooth(method=lm, se=FALSE)+
  stat_cor(method = "pearson", size = 5, show.legend = FALSE)

print(basic+ labs(x="Number of cells/recording",y="Percentage of active cells/recording (%)",color="Genotype")+
        theme(
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          legend.title = element_text(size=15),
          legend.text = element_text(size = 15)))


#-1 indicates a strong negative correlation : this means that every time x increases, y decreases (left panel figure)
#0 means that there is no association between the two variables (x and y) (middle panel figure)
#1 indicates a strong positive correlation : this means that y increases with x (right panel figure)
