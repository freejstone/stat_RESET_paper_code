
## logistic model with 2 side-info variable
args <- commandArgs(trailingOnly = TRUE)
ParamsRowIndex <- as.integer(args[1])
if(is.na(ParamsRowIndex)==1){
  ParamsRowIndex = 1 
}

#for parallelization
RNGkind(kind = "L'Ecuyer-CMRG")

####################################
## Libraries and sources
####################################
library("adaptMT")
library("splines")
library("knockoff")
library("SNPknock")
library("dplyr")
library("corpcor")
library("glmnet")
library("MASS")
library("R.matlab")
library("gam")
library("SNPknock")
library("randomForest")
library("devtools")
library("RCurl")
library("ggplot2")
library("httr")
library("mgcv")
library("flare")
library("expm")
library("rdist")

file_vec <- c("R_code/TDC_flex_c.R", 
              "R_code/ensemble_functions.R",
              "R_code/filter_glm.R",
              "R_code/filter_gam.R",
              "R_code/filter_randomForest.R",
              "R_code/filter_EM.R"
)
getfile <- sapply(file_vec,source)

####################################
## Parameters
####################################
set.seed(ParamsRowIndex)
n = 1000
p = 1600 # 30-by-30

alphalist = c(0.3,0.2,0.1,0.05)
Sigma = diag(rep(1,p))
amp = 25
settingName = "./results/simulation2_vary_ensemble"

####################################
## determining the signals
####################################
nonzero = c()
cor_hypothesis = expand.grid(-20:19,-20:19)
for (i in 1:p){
  xh = (cor_hypothesis[i,1])/9
  yh =(cor_hypothesis[i,2])/9
  dist = (xh^2+yh^2-2)^4-xh^3*yh^5
  dist3 = ((xh+1)^2+(yh-10/8)^2-0.015)
  dist4 = ((xh-1)^2+(yh+10/8)^2-0.015)
  
  if(dist<0 ) nonzero = c(nonzero,i)
  if(dist3<0 ) nonzero = c(nonzero,i)
  if(dist4<0 ) nonzero = c(nonzero,i)
}

# covariance matrix
Sigma = exp(-3*rdist::pdist(cor_hypothesis, metric = "euclidean", p = 2)^2)
k = length(nonzero)
beta0 = amp * (1:p %in% nonzero)*sign(rnorm(p)) / sqrt(n)
y.sample = function(X) rbinom(n,1,exp(X %*% beta0)/(1+exp(X %*% beta0)))

####################################
## Generating data
####################################
set.seed(ParamsRowIndex)
X = matrix(rnorm(n*p),n) %*% chol(Sigma)
y = y.sample(X)

####################################
## Vanilla  knockoff
####################################
Xk = create.gaussian(X,rep(0,p),Sigma)
mdl = cv.glmnet(cbind(X,Xk),y,alpha=1,family = "binomial")
cvlambda = mdl$lambda.min
beta = mdl$glmnet.fit$beta[,mdl$lambda ==mdl$lambda.min]
T = abs(beta[1:p])
T_tilde = abs(beta[(p+1):(2*p)])
T_max = pmax(T,T_tilde)
W = (T-T_tilde)
z = as.matrix(cor_hypothesis)

####################################
## KF
####################################
start.time = Sys.time()
q_vals = TDC_flex_c(W[W!=0] < 0, W[W!=0] > 0, score = abs(W[W!=0]))
end.time = Sys.time()
fdp = c()
power = c()
for (i in 1:length(alphalist)){
  fdp = c(fdp, sum((q_vals <= alphalist[i]) & (W[W!=0] > 0) & (beta0[W!=0] == 0))/max(sum((q_vals <= alphalist[i]) & (W[W!=0] > 0)), 1))
  power = c(power,sum((q_vals <= alphalist[i]) & (W[W!=0] > 0) & (beta0[W!=0] != 0))/max(k,1))
}
print(power)

savedir = paste(settingName,'/kn',k,'_',as.character(ParamsRowIndex),'.csv',sep = "")
write.csv(data.frame(fdp = fdp,power = power), savedir)

####################################
## filter RF
####################################
print('Doing RF adakn')
start.time = Sys.time()
res = filter_randomForest(W,z,alpha = alphalist,offset=1,reveal_prop = 0.1,mute = FALSE)
end.time = Sys.time()
fdp = c()
power = c()
for (i in 1:length(alphalist )){
  if(res$nrejs[[i]]>0){
    rej  = res$rejs[[i]]
    fdp = c(fdp,sum(beta0[rej]==0)/max(length(rej),1))
    power = c(power,sum(beta0[rej]!=0)/max(k,1))
    
  }else{
    fdp = c(fdp,0)
    power = c(power,0)
  }
}
print(power)

savedir = paste(settingName,'/adakn_RF',as.character(ParamsRowIndex),'.csv',sep = "")
write.csv(data.frame(fdp = fdp,power = power), savedir)


#####################################
### filter glm
#####################################
print('Doing glm adakn')
start.time = Sys.time()
res = filter_glm(W,z,alpha = alphalist,offset=1,reveal_prop = 0.1,mute = FALSE)
end.time = Sys.time()
fdp = c()
power = c()
for (i in 1:length(alphalist )){
  if(res$nrejs[[i]]>0){
    rej  = res$rejs[[i]]
    fdp = c(fdp,sum(beta0[rej]==0)/max(length(rej),1))
    power = c(power,sum(beta0[rej]!=0)/max(k,1))
    
  }else{
    fdp = c(fdp,0)
    power = c(power,0)
  }
}
print(power)

savedir = paste(settingName,'/adakn_glm', as.character(ParamsRowIndex),'.csv',sep = "")
write.csv(data.frame(fdp = fdp,power = power), savedir)


#####################################
### filter gam
#####################################
print('Doing gam adakn')
start.time = Sys.time()
res = filter_gam(W,z,alpha = alphalist,df_list=5,offset=1)
end.time = Sys.time()
fdp = c()
power = c()
for (i in 1:length(alphalist )){
  if(res$nrejs[[i]]>0){
    rej  = res$rejs[[i]]
    fdp = c(fdp,sum(beta0[rej]==0)/max(length(rej),1))
    power = c(power,sum(beta0[rej]!=0)/max(k,1))
    
  }else{
    fdp = c(fdp,0)
    power = c(power,0)
  }
}
print(power)

savedir = paste(settingName,'/adakn_gam', as.character(ParamsRowIndex),'.csv',sep = "")
write.csv(data.frame(fdp = fdp,power = power), savedir)



#####################################
### filter EM
#####################################
print('Doing EM Adakn')
start.time = Sys.time()
z = as.matrix(z)
s0 = quantile(abs(W)[W!=0], 0.1)
if (nrow(z) == 1 | ncol(z) == 1) {
  res = filter_EM(W,z,alpha =alphalist,offset=1,df = 2,s0 =s0,mute=FALSE)
} else {
  res = filter_EM(W,z,alpha =alphalist,offset=1,s0=s0,cutoff = 0,mute=FALSE)
}

end.time = Sys.time()
fdp = c()
power = c()
for (i in 1:length(alphalist )){
  if(res$nrejs[[i]]>0){
    rej  = res$rejs[[i]]
    fdp = c(fdp,sum(beta0[rej]==0)/max(length(rej),1))
    power = c(power,sum(beta0[rej]!=0)/max(k,1))
    
  }else{
    fdp = c(fdp,0)
    power = c(power,0)
  }
}
print(power)

savedir = paste(settingName,'/adakn_em', as.character(ParamsRowIndex),'.csv',sep = "")
write.csv(data.frame(fdp = fdp,power = power), savedir)

