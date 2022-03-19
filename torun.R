rm(list = ls())

# Set your work directory
source("minBIC_03182022.R")
library(dplyr)
library(readxl)


#####################################################################################
##############################BINARY DATA
# South Africa Heart Disease Data 
library(bestglm)
data(SAheart)
X = SAheart[,1:9]
X$famhist = X$famhist=="Present"
X$famhist = as.numeric(X$famhist)
y = SAheart[,10]

n=dim(X)[1]
p=dim(X)[2]
X = t((t(X)-apply(X,2,mean))/apply(X,2,sd))  # standardize the design matrix


# Specify if to include an intercept
logic = T

# Specify a decreasing sequence for tau
tau_seq=10^-seq(1,15,1)


pen = log(n)   # penalty parameter
start = rep(1,p)
weight_est = coef(glm(y~X-1, family="binomial"))
weight = log(1+1/tau_seq[1])/log(1+abs(weight_est)/tau_seq[1])
# weight = 1
run = minBIC(X,y,beta=start,interVal=1, inter=logic,tau3=tau_seq[1],lambda=pen,preweights=weight,thresh=1e-6, meth="IC", fam="binomial")
finish = run$BETA[dim(run$BETA)[1],]
	
for(i in 2:length(tau_seq)){
	# i=i+1
	# print(i)
	start = finish[-1]
	weight_est = finish[-1]
	weight_est[weight_est==0] = 1e-5
	weight = log(1+1/tau_seq[i])/log(1+abs(weight_est)/tau_seq[i])

	run = minBIC(X,y,beta=start,interVal=finish[1], inter=logic,tau3=tau_seq[i],lambda=pen,preweights=weight,thresh=1e-6, meth="IC", fam="binomial")
	finish = run$BETA[dim(run$BETA)[1],]
}
finish


#####################################################################################
##############################Continuous DATA


# Auto-mpg Data 
autoMpg = read.table("auto-mpg.data", sep="", header=F, dec=".")
colnames(autoMpg) = c("mpg","cylinders","displacement","horsepower","weight","acceleration","model year", "origin", "car name")

autoMpg$horsepower = as.numeric(autoMpg$horsepower)
idsNA = which(is.na(autoMpg$horsepower))
autoMpg = autoMpg[-idsNA,]

X = as.matrix(autoMpg[,-c(1,8,9)])
y = autoMpg[,"mpg"]


X = t((t(X)-apply(X,2,mean))/apply(X,2,sd))


# Specify intercept setting
logic=T

p=dim(X)[2]
n=dim(X)[1]


##########################
# Analysis with our method
tau_seq=10^-seq(1,15,1)
	
	
pen_seq = log(n)
start = rep(1,p)
weight_est = coef(lm(y~X-1))
weight = log(1+1/tau_seq[1])/log(1+abs(weight_est)/tau_seq[1])
run = minBIC(X,y,beta=start,interVal=mean(y), inter=logic,tau3=tau_seq[1],lambda=pen_seq,preweights=weight,thresh=1e-6, meth="IC", fam="gaussian")
finish = run$BETA[dim(run$BETA)[1],]

	
for(i in 2:length(tau_seq)){
	
	start = finish[-1]	
	weight_est = finish[-1]
	weight_est[weight_est==0] = 1e-5
	weight = log(1+1/tau_seq[i])/log(1+abs(weight_est)/tau_seq[i])
	
	run = minBIC(X,y,beta=start,interVal=finish[1], inter=logic,tau3=tau_seq[i],lambda=pen_seq,preweights=weight,thresh=1e-6, meth="IC", fam="gaussian")
	finish = run$BETA[dim(run$BETA)[1],]

}
finish


