rm(list = ls())
library(readxl)
library(dplyr)

# Set your work directory
# setwd("/Users/yangyuan/Desktop/BSS_star-main")
source("minBIC_03182022.R")
source("l0norm_linear.R")
source("BIC_logistic.R")
source("BIC_linear.R")
tau_seq = 10^-seq(1,15,1)



# # # # # # # # # # 
# # CONTINUOUS DATA

# 1. Prostate Data  *****  (small p)
library(bestglm)
data(zprostate)

X=as.matrix(zprostate[,1:8])
y=zprostate[,9]
p = dim(X)[2]

fit1_inter = BIC_linear(X,y,rep(1,p),Intercept=T,tau_seq)
fit1_noInter = BIC_linear(X,y,rep(1,p),Intercept=F,tau_seq)


# 2. Auto-mpg Data  *****  (small p)
autoMpg = read.table("auto-mpg.data", sep="", header=F, dec=".")
colnames(autoMpg) = c("mpg","cylinders","displacement","horsepower","weight","acceleration","model year", "origin", "car name")

autoMpg$horsepower = as.numeric(autoMpg$horsepower)
idsNA = which(is.na(autoMpg$horsepower))
autoMpg = autoMpg[-idsNA,]

X = as.matrix(autoMpg[,-c(1,8,9)])
y = autoMpg[,"mpg"]
p = dim(X)[2]

fit2_inter = BIC_linear(X,y,rep(1,p),Intercept=T,tau_seq)
fit2_noInter = BIC_linear(X,y,rep(1,p),Intercept=F,tau_seq)


# 3. real estate data  *****  (small p)
realEstate = read_excel("RealEstateValuationDataSet.xlsx")
realEstate = as.matrix(realEstate)
X = realEstate[,3:7]
y = realEstate[,8]
colnames(X) = c("houseAge","distMRT","#convStores","lat","long")
p = dim(X)[2]

fit3_inter = BIC_linear(X,y,rep(1,p),Intercept=T,tau_seq)
fit3_noInter = BIC_linear(X,y,rep(1,p),Intercept=F,tau_seq)




# # # # # # # # # # 
# # BINARY DATA

# Load data and construct data

# 1. SA heart disease data
library(bestglm)
data(SAheart)
X = SAheart[,1:9]
X$famhist = X$famhist=="Present"
X$famhist = as.numeric(X$famhist)
y = SAheart[,10]

p=dim(X)[2]

fit4_inter = BIC_logistic(X,y,beta0=rep(1,p), Intercept=T,tau_seq)
fit4_noInter = BIC_logistic(X,y,beta0=rep(1,p), Intercept=F,tau_seq)



# 2. Indian Liver Patient Dataset *****
IL = read.csv("Indian Liver Patient Dataset (ILPD).csv", header = F)
myifelse = function(x) ifelse(x=="Female",1,0)
IL = mutate_at(IL, "V2", myifelse)

X = IL[,1:10]
X[is.na(X[,10]),10] = mean(X[,10], na.rm=T)
y = IL[,11]
y[y==2] = 0

p=dim(X)[2]

fit5_inter = BIC_logistic(X,y,beta0=rep(1,p), Intercept=T,tau_seq)
fit5_noInter = BIC_logistic(X,y,beta0=rep(1,p), Intercept=F,tau_seq)



# 3. Hungarian heart disease data *****
hungHeart = read.table("reprocessed.hungarian.data", sep=" ")
hungHeart = hungHeart[-295,]

X = hungHeart[,1:13]
y = as.numeric(hungHeart$V14!=0)

p=dim(X)[2]

fit6_inter = BIC_logistic(X,y,beta0=rep(1,p), Intercept=T,tau_seq)
fit6_noInter = BIC_logistic(X,y,beta0=rep(1,p), Intercept=F,tau_seq)


