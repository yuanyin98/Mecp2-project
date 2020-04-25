my_data <- read.csv(file.choose())

library("dplyr")
sample_n(my_data, 10)

# F-test
res.ftest <- var.test(dv ~ grp, data = my_data, alternative = "two.sided")
res.ftest

#Interpretation of the result
# A=KO, B=WT
#alternative hypothesis: true ratio of variances is not equal to 1
#The p-value of F-test is p = 0.1203: greater than the significance level 0.05, accept null hypothesis
#In conclusion, no significant difference in variances