my_data <- read.csv(file.choose())

library("dplyr")
sample_n(my_data, 10)

# F-test
res.ftest <- var.test(dv ~ grp, data = my_data, alternative = "less")
res.ftest

#Interpretation of the result
# A=KO, B=WT
#alternative hypothesis: true ratio of variances is less than 1
#The p-value of F-test is p = 0.003624: smaller than the significance level 0.05, reject null hypothesis
#In conclusion, KO has significantly smaller variances than WT