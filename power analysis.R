# http://powerandsamplesize.com/Calculators/Compare-2-Means/2-Sample-1-Sided

muA=23
muB=15
kappa=10/11
#kappa=n(A)/n(B); A=WT, B=KO
sdA=13
sdB=12
alpha=0.05
beta=0.20

#sample size calculation=nA
(nA=(sdA^2+sdB^2/kappa)*((qnorm(1-alpha)+qnorm(1-beta))/(muA-muB))^2)
ceiling(nA) # 32
#nB=kappa*nA
nB <- kappa*nA
ceiling(nB) # 29

#power calculation=Power
z=(muA-muB)/sqrt(sdA^2+sdB^2/kappa)*sqrt(nA)
(Power=pnorm(z-qnorm(1-alpha)))
## Note: Rosner example on p.303 is for 2-sided test.
## These formulas give the numbers in that example
## after dividing alpha by 2, to get 2-sided alpha.