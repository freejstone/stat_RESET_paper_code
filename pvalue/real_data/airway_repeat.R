source('../../R_code/ensemble_functions.R')
source('../../R_code/helper_functions.R')

files <- c(
  "airway"
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

alphas <- c(0.1)

set.seed(22052024)

for (file in files) {
  print(paste("File:",file))
  out <- preprocess(file)
  x <- out$x
  pvals <- out$pvals
  W = abs(qnorm(pvals))*((pvals < 0.5) - (pvals > 0.5))
  W[W == Inf] = max(W[W!=Inf])
  W[W == -Inf] = min(W[W!=-Inf])
  mult = 1
  decoy_int = c(1/2, 1)
  target_upper = 1/2
  if (grepl('fmri', file)) {
    x = x[pvals <= 0.9,]
    pvals = pvals[pvals <= 0.9]
    W = pvals
    W[pvals > 0.3] = 0.3 - (pvals[pvals > 0.3] - 0.3)/2
    W = abs(qnorm(W))
    W = W*((pvals <= 0.3) - (pvals > 0.3))
    W[W == Inf] = max(W[W!=Inf])
    W[W == -Inf] = min(W[W!=-Inf])
    x = data.frame(x)
    mult = 2
    decoy_int = c(0.3, 0.9)
    target_upper = 0.3
  }
  reps = 10
  if (file %in% c('airway', 'bottomly', 'pasilla')) {
    W=W[x>=0]
    pvals = pvals[x>=0]
    x=x[x>=0]
    x = as.matrix(x)
  } else if (file %in% c('proteomics')) {
    x = x$x1
    x = as.matrix(x)
  } else if (file %in% c('microbiome_enigma_ph','microbiome_enigma_al', 'adipose_subcutaneous', 'adipose_visceral_omentum')) {
    x = as.matrix(x)
  }
  
  out = data.frame(alpha = numeric(), 
                   type = character(), 
                   power = numeric(), 
                   time = numeric(),
                   flag = logical(),
                   seed = numeric())
  
  if (file %in% c("adipose_subcutaneous", "adipose_visceral_omentum")) {
    seeds = 22072024
  } else {
    seeds = 22072024:(22072024 + 99)
  }
  
  for (seed in seeds) {
    for (alpha in alphas) {
      start = Sys.time()
      res = filter_ensemble_RESET(W, x, decoy_int = decoy_int, target_upper = target_upper,
                                  verbose = TRUE, test_alpha = alpha, seed = seed,
                                  mult = mult, reps = reps, num_cores = 20)
      end = Sys.time()
      out = rbind(out, list(alpha = alpha, type = 'reset_ensemble', power = sum(res$q_vals <= alpha & res$Labels == 1), 
                            time = difftime(end, start, units='mins'), flag = res$flag, seed = seed))
    }
  }
  write.csv(out, paste('results/', file,'_rep_reset_all.csv', sep = ''))
}
