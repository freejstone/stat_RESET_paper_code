library(zap)
library(splines)
RNGkind("L'Ecuyer-CMRG")
files <- c(
  "airway",
  "pasilla",
  "bottomly",
  "microbiome_enigma_ph",
  "microbiome_enigma_al",
  "proteomics",
  "fmri_auditory","fmri_imagination",
  "adipose_subcutaneous",
  "adipose_visceral_omentum"
)
data_path <- "./data_files/"


#extract the zap thresholding routine at the end of their function (so I don't need to rerun it for each alpha separately)
zap_thresholding = function(blfdr_vec, blfdr_flip_vec, alpha, ep = 1) {
  bottom <- sapply(X = blfdr_vec, FUN = function(cutpt, vec) {
    max(1, sum(vec <= cutpt))
  }, vec = blfdr_vec)
  top_BC <- ep + sapply(X = blfdr_vec, FUN = function(cutpt, 
                                                      vec) {
    sum(vec <= cutpt)
  }, vec = blfdr_flip_vec)
  FDR_est <- (top_BC/bottom)
  cutpt_max_BC <- max(c(blfdr_vec[FDR_est <= alpha], -Inf))
  rej_index_blfdr_BC <- which(blfdr_vec <= cutpt_max_BC)
  return(rej_index_blfdr_BC)
}

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
for(file in files){
  print(paste("File:",file))
  out <- preprocess(file)
  x <- out$x
  pvals <- out$pvals
  if (file %in% 'bottomly') {
    pvals = pvals[x > 0]
    x =as.matrix(x[x > 0])
  }
  
  if (file %in% c('airway', 'pasilla')) {
    pvals = pvals[x>=0]
    x=as.matrix(x[x>=0])
  }
  if (ncol(x) == 1) {
    colnames(x) = 'x'
    fs = 'ns(x, df=6)'
  }
  if (ncol(x) == 2) {
    fs = 'ns(x1, df=6) + ns(x2, df=6)'
  }
  if (ncol(x) == 4) {
    fs = 'ns(x1, df=6) + ns(x2, df=6) + ns(x3, df=6) + x4' #last one doesn't allow for spline transformation
  }
  if (grepl('fmri', file)) {
    fs = "ns(x1,df=6) + ns(x2, df=6) + ns(x3, df = 6) + ns(x4,df=6)"
  }
  
  df_final = data.frame(alpha = numeric(), 
                        type = character(), 
                        power = numeric(), 
                        time = numeric(), 
                        seed = numeric())
  set.seed(24062024)
  z = qnorm(pvals/2)*sign(rnorm(length(pvals)))
  z[z == -Inf] = min(z[z!= -Inf])
  z[z == Inf] = max(z[z!= Inf])
  start = Sys.time()
  x = model.matrix(formula(paste("~", fs,"-1")), data = data.frame(x))
  res_zap_asymp = zap_asymp(z, as.matrix(x), alpha = Inf) #Get mirror statistics
  end = Sys.time()
  mirror_zap_statistics = res_zap_asymp$mirror_stat
  zap_statistics = res_zap_asymp$stat
  res_zap_asymp = lapply(alphas, zap_thresholding, blfdr_vec = zap_statistics, blfdr_flip_vec = mirror_zap_statistics)
  for (i in 1:length(alphas)) {
    alpha = alphas[i]
    df_final = rbind(df_final, list(alpha = alpha, type = 'zap', power = length(res_zap_asymp[i][[1]]), time = difftime(end, start, units='mins'), seed = 24062024))
  }
  write.csv(df_final, paste('results/', file,'_zap.csv', sep = ''))
}