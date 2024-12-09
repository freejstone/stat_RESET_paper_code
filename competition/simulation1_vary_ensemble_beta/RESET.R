
## linear model with 2 side-info variable
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
####################################
set.seed(ParamsRowIndex)
n = 2*1000
p = 2*900


alphalist = c(0.3,0.2,0.1,0.05)
rho = 0.5

Sigma = toeplitz(rho^(0:(p-1)))
amp = 3.5
sigprob = rep(0,p)
sigprob[1:(2*k)] = 1/(1:(2*k))^2/(sum(1/(1:(2*k))^2))
nonzero = sample(1:p,k,prob = sigprob)
beta0 = amp * (1:p %in% nonzero)*sign(rnorm(p)) / sqrt(n)
y.sample = function(X) X%*%beta0+rnorm(n,0,1)
settingName = "./results/simulation1_vary_ensemble_beta"

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

z = 1:p
Xk = knockoffHMM(X, pInit, Q,pEmit,seed = ParamsRowIndex+14112024)
mdl = cv.glmnet(cbind(X,Xk),y,alpha=1)
cvlambda = mdl$lambda.min
beta = mdl$glmnet.fit$beta[,mdl$lambda ==mdl$lambda.min]
T = abs(beta[1:p])
T_tilde = abs(beta[(p+1):(2*p)])
T_max = pmax(T,T_tilde)
W = T-T_tilde

####################################
## Generating side-information
####################################

get_z2 = function(beta0, sig = 2) {
  z = mvrnorm(p, mu = c(0, 0), Sigma = matrix(c(1, 0, 0, 1), nrow = 2))
  z[beta0 != 0,] = z[beta0!=0,] + matrix(2, nrow = sum(beta0!=0), ncol = 2)*sample(c(-1, 1), size = sum(beta0!=0), replace = TRUE)
  return(z)
}

z = get_z2(beta0)

####################################
## ensemble RESET 
####################################

print('Doing RF')
fdp = c()
power = c()
for (i in 1:length(alphalist)){
  res = filter_ensemble_RESET(W, z, methods = 'rf', verbose = TRUE, test_alpha = alphalist[i], 
                              seed = ParamsRowIndex, mult = 1, reps = 10, n_nodes = c(2, 5, 10), num_cores = 10, get_nn = 20, train_alpha = 0.5)
  fdp = c(fdp, sum((res$q_vals <= alphalist[i]) & (res$pseudo_Labels == 1) & (res$Labels == 1) & (beta0[W !=0] == 0))/max(sum((res$q_vals <= alphalist[i]) & (res$Labels == 1)), 1))
  power = c(power,sum((res$q_vals <= alphalist[i]) & (res$pseudo_Labels == 1) & (res$Labels == 1) & (beta0[W !=0] != 0))/max(k,1))
}
print(power)


savedir = paste(settingName,'/RESET_RF_beta',k,'_',as.character(ParamsRowIndex),'.csv',sep = "")
write.csv(data.frame(fdp = fdp,power = power), savedir)


####################################
## ensemble RESET 
####################################

print('Doing NN')
fdp = c()
power = c()
for (i in 1:length(alphalist)){
  res = filter_ensemble_RESET(W, z, methods = 'nn', verbose = TRUE, test_alpha = alphalist[i], 
                              seed = ParamsRowIndex, mult = 1, reps = 10, n_nodes = c(2, 5, 10), num_cores = 10, get_nn = 20, train_alpha = 0.5)
  fdp = c(fdp, sum((res$q_vals <= alphalist[i]) & (res$pseudo_Labels == 1) & (res$Labels == 1) & (beta0[W !=0] == 0))/max(sum((res$q_vals <= alphalist[i]) & (res$Labels == 1)), 1))
  power = c(power,sum((res$q_vals <= alphalist[i]) & (res$pseudo_Labels == 1) & (res$Labels == 1) & (beta0[W !=0] != 0))/max(k,1))
}
print(power)


savedir = paste(settingName,'/RESET_NN_beta',k,'_',as.character(ParamsRowIndex),'.csv',sep = "")
write.csv(data.frame(fdp = fdp,power = power), savedir)


####################################
## ensemble RESET 
####################################

print('Doing GAM')
fdp = c()
power = c()
for (i in 1:length(alphalist)){
  res = filter_ensemble_RESET(W, z, methods = 'gam', verbose = TRUE, test_alpha = alphalist[i], 
                              seed = ParamsRowIndex, mult = 1, reps = 10, n_nodes = c(2, 5, 10), num_cores = 10, get_nn = 20, train_alpha = 0.5)
  fdp = c(fdp, sum((res$q_vals <= alphalist[i]) & (res$pseudo_Labels == 1) & (res$Labels == 1) & (beta0[W !=0] == 0))/max(sum((res$q_vals <= alphalist[i]) & (res$Labels == 1)), 1))
  power = c(power,sum((res$q_vals <= alphalist[i]) & (res$pseudo_Labels == 1) & (res$Labels == 1) & (beta0[W !=0] != 0))/max(k,1))
}
print(power)


savedir = paste(settingName,'/RESET_GAM_beta',k,'_',as.character(ParamsRowIndex),'.csv',sep = "")
write.csv(data.frame(fdp = fdp,power = power), savedir)


####################################
## ensemble RESET 
####################################

print('Doing ensemble')
fdp = c()
power = c()
for (i in 1:length(alphalist)){
  res = filter_ensemble_RESET(W, z, methods = c('rf', 'nn', 'gam'), verbose = TRUE, test_alpha = alphalist[i], 
                              seed = ParamsRowIndex, mult = 1, reps = 10, n_nodes = c(2, 5, 10), num_cores = 10, get_nn = 20, train_alpha = 0.5)
  fdp = c(fdp, sum((res$q_vals <= alphalist[i]) & (res$pseudo_Labels == 1) & (res$Labels == 1) & (beta0[W !=0] == 0))/max(sum((res$q_vals <= alphalist[i]) & (res$Labels == 1)), 1))
  power = c(power,sum((res$q_vals <= alphalist[i]) & (res$pseudo_Labels == 1) & (res$Labels == 1) & (beta0[W !=0] != 0))/max(k,1))
}
print(power)


savedir = paste(settingName,'/RESET_ensemble_beta',k,'_',as.character(ParamsRowIndex),'.csv',sep = "")
write.csv(data.frame(fdp = fdp,power = power), savedir)
