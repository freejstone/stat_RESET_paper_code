library(stepdownfdp)
source('../../competition/R_code/ensemble_functions.R')
source('../functions/guo.R')

files <- c(
  "airway",
  "pasilla",
  "bottomly",
  "microbiome_enigma_ph",
  "microbiome_enigma_al",
  "proteomics",
  "fmri_auditory","fmri_imagination"
#  "adipose_subcutaneous",
 # "adipose_visceral_omentum"
)
data_path <- "./data_files/"

preprocess <- function(file){
  if(grepl("adipose",file)){
    x <- read.csv(paste0(data_path,file,"_x.csv"),header=F)
    pvals <- read.csv(paste0(data_path,file,"_p.csv"),header=F)$V1
  }else{
    
    
    data <- read.csv(paste0(data_path,file),header=T)
    pvals <- data$p_val
    # If the columns are not labeled
    if(is.null(pvals)){
      pvals <- data[,1]
      if(max(pvals)>1 | min(pvals)<0){
        stop("Invalid dataset found.")
      } 
      # Rest of columns are covariates
      x <- data[,2:ncol(data)]
    }else{
      x <- data[,!names(data) =="p_val"]
    }
    if(is.null(ncol(x))){
      x <- data.frame(x=x)
    }
    x <- Filter(function(y) sd(y) != 0, x)
  }
  
  colnames(x) <- paste0("x",seq(1:ncol(x)))
  return(list(x=x,pvals=pvals))
}

alphas <- c(0.01,0.05,0.1,0.2)

set.seed(22052024)

for (file in files) {
  print(paste("File:",file))
  out <- preprocess(file)
  x <- out$x
  pvals <- out$pvals
  pvals1 <- pvals
  
  x = as.matrix(x)[pvals <= 0.9,]
  pvals = pvals[pvals <= 0.9]
  W = pvals
  W[pvals > 0.3] = 0.3 - (pvals[pvals > 0.3] - 0.3)/2
  W = abs(qnorm(W))
  W = W*((pvals <= 0.3) - (pvals > 0.3))
  W[W == Inf] = max(W[W!=Inf])
  W[W == -Inf] = min(W[W!=-Inf])
  x = data.frame(x)
  mult = 2
  reps = 10
  if (file %in% c('airway', 'bottomly', 'pasilla')) {
    W=W[x>=0]
    pvals = pvals[x>=0]
    x=x[x>=0]
    x = as.matrix(x)
  } else if (file %in% c('proteomics')) {
    x = as.matrix(x)
  } else if (file %in% c('microbiome_enigma_ph','microbiome_enigma_al', 'adipose_subcutaneous', 'adipose_visceral_omentum')) {
    x = as.matrix(x)
  }
  
  out = data.frame(alpha = numeric(), 
                   type = character(), 
                   power = numeric(), 
                   time = numeric(), 
                   conf = numeric(),
                   seed = numeric())
  
  for (seed in 22072024:(22072024 + 9)) {
    for (alpha in alphas) {
      start.time = Sys.time()
      res = filter_ensemble_RESET(W, x, verbose = TRUE, test_alpha = alpha, seed = seed, mult = mult, reps = reps, get_nn = NULL, num_cores = 10)
      scores = rank(res$score[res$pseudo_Labels == 1], ties.method = 'random')
      labels = res$Labels[res$pseudo_Labels == 1]
      end.time = Sys.time()
      
      for (conf in c(0.1, 0.2, 0.5)) {
        res = fdp_sd(cbind(scores, labels), alpha = alpha, conf = conf, c = 1/2, lambda = 1/2, procedure = 'coinflip')
        total.time = difftime(end.time, start.time, units='mins')
        power = length(res$discoveries)
        out = rbind(out, list(alpha = alpha, type = 'reset_ensemble', power = power, time = total.time, conf = conf, seed = seed))
      }
      
      for (conf in c(0.1, 0.2, 0.5)) {
        start.time = Sys.time()
        res = fdp_sd(cbind(abs(W)[W!=0], sign(W)[W!=0]), alpha = alpha, conf = conf, c = 1/3, lambda = 1/3, procedure = 'coinflip')
        end.time = Sys.time()
        total.time = difftime(end.time, start.time, units='mins')
        power = length(res$discoveries_ind)
        out = rbind(out, list(alpha = alpha, type = 'FDP-SD', power = power, time = total.time, conf = conf, seed = seed))
      }
      
      for (conf in c(0.1, 0.2, 0.5)) {
        start.time = Sys.time()
        res = fdp_sdp_guo(pvals1, alpha, conf)
        end.time = Sys.time()
        total.time = difftime(end.time, start.time, units='mins')
        power = sum(res$which_less)
        out = rbind(out, list(alpha = alpha, type = 'GR-SD', power = power, time = total.time, conf = conf, seed = seed))
      }
      
    }
  }
  
  write.csv(out, paste('results/', file,'_fdp_reset_all.csv', sep = ''))
}
