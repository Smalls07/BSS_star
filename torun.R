rm(list = ls())

# Set your work directory
setwd("/Users/yangyuan/Desktop/BSS_star-main")
source("minBIC_03182022.R")
source("fit0norm_leastSquares.R")
source("BIC_cont.R")
source("BIC_binary.R")




# # # # # # # # # # 
# # CONTINUOUS DATA


# Load data and construct data
autoMpg = read.table("auto-mpg.data", sep="", header=F, dec=".")
colnames(autoMpg) = c("mpg","cylinders","displacement","horsepower","weight","acceleration","model year", "origin", "car name")

autoMpg$horsepower = as.numeric(autoMpg$horsepower)
idsNA = which(is.na(autoMpg$horsepower))
autoMpg = autoMpg[-idsNA,]

X = as.matrix(autoMpg[,-c(1,8,9)])
y = autoMpg[,"mpg"]


# Specify function parameters
p = dim(X)[2]
start = rep(1, p)
tau_seq = 10^-seq(1,15,1)


# 1. 0-norm penalized least squares
exampleRun1 = fit0norm_leastSquares(X,y,start,fitInter=F,lambda=6,tau_seq)


# 2. minimize BIC
exampleRun2 = BIC_cont(X,y,start,fitInter=T,tau_seq)




# # # # # # # # # # 
# # BINARY DATA

# Load data and construct data
library(bestglm)
data(SAheart)
X = SAheart[,1:9]
X$famhist = X$famhist=="Present"
X$famhist = as.numeric(X$famhist)
y = SAheart[,10]

# Specify function parameters
p=dim(X)[2]
tau_seq = 10^-seq(1,15,1)

exampleRun3 = BIC_binary(X,y,start=rep(1,p),fitInter=T,tau_seq)

