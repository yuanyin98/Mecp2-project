my_data <- read.csv(file.choose())

library("dplyr")
sample_n(my_data, 10)

# F-test
res.ftest <- var.test(dv ~ grp, data = my_data, alternative = "less")
res.ftest

#Interpretation of the result
#A=KO, B=WT, null: A>B
#alternative hypothesis: true ratio of variances (A:B) is less than 1 = B>A
#The p-value of F-test is p = 0.0001557: smaller than the significance level 0.05, reject null hypothesis
#In conclusion, the variance of A is significantly smaller than the variance of B