library(adaptMT)
library(splines)
library(mgcv)
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
  
  if (file %in% c('airway', 'pasilla','bottomly')) {
    pvals = pvals[x>=0]
    x=as.matrix(x[x>=0])
  }
  if (ncol(x) == 1) {
    colnames(x) = 'x'
    fs = 'ns(x, df=5)'
  }
  if (ncol(x) == 2) {
    fs = 'ns(x1, df=5) + ns(x2, df=5)'
  }
  if (ncol(x) == 4) {
    fs = 'ns(x1, df=5) + ns(x2, df=5) + ns(x3, df=5) + x4' #last one doesn't allow for spline transformation
  }
  if (grepl('fmri', file)) {
    fs = "ns(x1, 5) + ns(x2, df=5) + ns(x3, df = 5) + ns(x4,df=5)"
  }
  
  df_final = data.frame(alpha = numeric(), 
                        type = character(), 
                        power = numeric(), 
                        time = numeric(), 
                        seed = numeric())
  set.seed(24062024)
  start = Sys.time()
  out <- adapt_gam(pvals=pvals, x= data.frame(x),
                   pi_formulas = fs,
                   mu_formulas = fs,
                   dist = beta_family(),
                   alphas = alphas)
  end = Sys.time()
  for (i in 1:length(alphas)) {
    alpha = alphas[i]
    df_final = rbind(df_final, list(alpha = alpha, type = 'AdaPT', power = out$nrejs[i], time = difftime(end, start, units='mins'), seed = 24062024))
  }
  write.csv(df_final, paste('results/', file,'_adapt.csv', sep = ''))
  
}