
## linear model with 1 side-info variable
args <- commandArgs(trailingOnly = TRUE)
ParamsRowIndex <- as.integer(args[1])
k <- as.integer(args[2])
if(is.na(ParamsRowIndex)==1){
  ParamsRowIndex = 1 
}
if(is.na(k)==1){
  k = 150
}

#for parallelization
RNGkind(kind = "L'Ecuyer-CMRG")

print(paste('Using seed:', ParamsRowIndex))

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
## generating model
set.seed(ParamsRowIndex)
n = 1000
p = 900 


alphalist = c(0.3,0.2,0.1,0.05)
rho = 0.5

Sigma = toeplitz(rho^(0:(p-1)))
amp = 3.5
sigprob = rep(0,p)
sigprob[1:(2*k)] = 1/(1:(2*k))^2/(sum(1/(1:(2*k))^2))
nonzero = sample(1:p,k,prob = sigprob)
beta0 = amp * (1:p %in% nonzero)*sign(rnorm(p)) / sqrt(n)
y.sample = function(X) X%*%beta0+rnorm(n,0,1)
settingName = "./results/simulation1_vary_ensemble"

####################################
## HMM parameters
####################################
# Number of possible states for each variable
K=5;M=3;  
# Marginal distribution for the first variable
pInit = rep(1/K,K)
# Create p-1 transition matrices
Q = array(stats::runif((p-1)*K*K),c(p-1,K,K))
for(j in 1:(p-1)) { Q[j,,] = Q[j,,] / rowSums(Q[j,,]) }
pEmit = array(stats::runif(p*M*K),c(p,M,K))
for(j in 1:p) { pEmit[j,,] = pEmit[j,,] / rowSums(pEmit[j,,]) }

####################################
## Generating data
####################################
X = sampleHMM(pInit, Q, pEmit, n=n)
y = y.sample(X)

z <-1:p
Xk = knockoffHMM(X, pInit, Q,pEmit,seed = ParamsRowIndex+24601)
mdl = cv.glmnet(cbind(X,Xk),y,alpha=1)
cvlambda = mdl$lambda.min
beta = mdl$glmnet.fit$beta[,mdl$lambda ==mdl$lambda.min]
T = abs(beta[1:p])
T_tilde = abs(beta[(p+1):(2*p)])
T_max = pmax(T,T_tilde)
W = T-T_tilde

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

savedir = paste(settingName,'/adakn_em_test', k,'_', as.character(ParamsRowIndex),'.csv',sep = "")
write.csv(data.frame(fdp = fdp,power = power), savedir)


#####################################
### filter RF
#####################################
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

savedir = paste(settingName,'/adakn_RF_test',k,'_',as.character(ParamsRowIndex),'.csv',sep = "")
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

savedir = paste(settingName,'/adakn_glm_test', k,'_', as.character(ParamsRowIndex),'.csv',sep = "")
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

savedir = paste(settingName,'/adakn_gam_test', k,'_', as.character(ParamsRowIndex),'.csv',sep = "")
write.csv(data.frame(fdp = fdp,power = power), savedir)

