
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
## ensemble RESET 
####################################

print('Doing RF')
fdp = c()
power = c()
for (i in 1:length(alphalist)){
  res = filter_ensemble_RESET(W, z, methods = 'rf', verbose = TRUE, test_alpha = alphalist[i], 
                              seed = ParamsRowIndex, mult = 1, reps = 10, n_nodes = c(2, 5, 10), decays = c(0, 0.1, 1), num_cores = 10, get_nn = 20, train_alpha = 0.5)
  fdp = c(fdp, sum((res$q_vals <= alphalist[i]) & (res$pseudo_Labels == 1) & (res$Labels == 1) & (beta0[W!=0] == 0))/max(sum((res$q_vals <= alphalist[i]) & (res$Labels == 1)), 1))
  power = c(power,sum((res$q_vals <= alphalist[i]) & (res$pseudo_Labels == 1) & (res$Labels == 1) & (beta0[W!=0] != 0))/max(k,1))
}
print(power)


savedir = paste(settingName,'/RESET_RF_zero',k,'_',as.character(ParamsRowIndex),'.csv',sep = "")
write.csv(data.frame(fdp = fdp,power = power), savedir)


####################################
## ensemble RESET 
####################################

print('Doing NN')
fdp = c()
power = c()
for (i in 1:length(alphalist)){
  res = filter_ensemble_RESET(W, z, methods = 'nn', verbose = TRUE, test_alpha = alphalist[i], 
                              seed = ParamsRowIndex, mult = 1, reps = 10, n_nodes = c(2, 5, 10), decays = c(0, 0.1, 1), num_cores = 10, get_nn = 20, train_alpha = 0.5)
  fdp = c(fdp, sum((res$q_vals <= alphalist[i]) & (res$pseudo_Labels == 1) & (res$Labels == 1) & (beta0[W!=0] == 0))/max(sum((res$q_vals <= alphalist[i]) & (res$Labels == 1)), 1))
  power = c(power,sum((res$q_vals <= alphalist[i]) & (res$pseudo_Labels == 1) & (res$Labels == 1) & (beta0[W!=0] != 0))/max(k,1))
}
print(power)


savedir = paste(settingName,'/RESET_NN_zero',k,'_',as.character(ParamsRowIndex),'.csv',sep = "")
write.csv(data.frame(fdp = fdp,power = power), savedir)


####################################
## ensemble RESET 
####################################

print('Doing GAM')
fdp = c()
power = c()
for (i in 1:length(alphalist)){
  res = filter_ensemble_RESET(W, z, methods = c('gam'), verbose = TRUE, test_alpha = alphalist[i], 
                              seed = ParamsRowIndex, mult = 1, reps = 10, n_nodes = c(2, 5, 10), decays = c(0, 0.1, 1), num_cores = 10, get_nn = 20, train_alpha = 0.5)
  fdp = c(fdp, sum((res$q_vals <= alphalist[i]) & (res$pseudo_Labels == 1) & (res$Labels == 1) & (beta0[W!=0] == 0))/max(sum((res$q_vals <= alphalist[i]) & (res$Labels == 1)), 1))
  power = c(power,sum((res$q_vals <= alphalist[i]) & (res$pseudo_Labels == 1) & (res$Labels == 1) & (beta0[W!=0] != 0))/max(k,1))
}
print(power)


savedir = paste(settingName,'/RESET_GAM_zero',k,'_',as.character(ParamsRowIndex),'.csv',sep = "")
write.csv(data.frame(fdp = fdp,power = power), savedir)


####################################
## ensemble RESET 
####################################

print('Doing ensemble')
fdp = c()
power = c()
for (i in 1:length(alphalist)){
  res = filter_ensemble_RESET(W, z, methods = c('nn', 'rf', 'gam'), verbose = TRUE, test_alpha = alphalist[i], 
                              seed = ParamsRowIndex, mult = 1, reps = 10, n_nodes = c(2, 5, 10), decays = c(0, 0.1, 1), num_cores = 10, get_nn = 20, train_alpha = 0.5)
  fdp = c(fdp, sum((res$q_vals <= alphalist[i]) & (res$pseudo_Labels == 1) & (res$Labels == 1) & (beta0[W!=0] == 0))/max(sum((res$q_vals <= alphalist[i]) & (res$Labels == 1)), 1))
  power = c(power,sum((res$q_vals <= alphalist[i]) & (res$pseudo_Labels == 1) & (res$Labels == 1) & (beta0[W!=0] != 0))/max(k,1))
}
print(power)


savedir = paste(settingName,'/RESET_ensemble_zero',k,'_',as.character(ParamsRowIndex),'.csv',sep = "")
write.csv(data.frame(fdp = fdp,power = power), savedir)

