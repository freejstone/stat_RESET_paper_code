library(AdaPTGMM)
library(RadaFDR)
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
  
  # Reduce model selection for large files
  if(grepl("adipose",file)){
    nclasses <- 4
    spline_df <- 3
  }else{
    nclasses <- c(2,4,6)
    spline_df <- c(3,5)
  }
  
  
  # These p-value distributions violate assumptions, we resample the p-values from the interval [0.3,0.9]
  # randomize <- (grepl("airway",file) | grepl("pasilla",file))
  # seems more natural to delete them.
  if (file %in% c('airway', 'pasilla', 'bottomly')) {
    pvals = pvals[x>=0]
    x=x[x>=0]
    x=data.frame(x1 = x)
  }
  
  if(ncol(x)==1){
    spline_formulas <- paste0("splines::ns(x1, df = ",seq(2,8,2) ,")")
  }else{
    spline_formulas <- c(0)
    for(df in spline_df){
      curr_formula <- paste0("splines::ns(x",1,", df = ",df ,")")
      for(ind in 2:ncol(x)){
        if(grepl("adipose",file) & ind == 4){
          df <- 1
        }
        curr_formula <- paste(curr_formula,"+",paste0("splines::ns(x",ind,", df = ",df ,")"))
      }
      spline_formulas <- c(spline_formulas,curr_formula)
    }
    spline_formulas <- spline_formulas[2:length(spline_formulas)]
    
  }
  
  target_alpha_level <- 0.1
  if(grepl("microbiome", file, fixed = TRUE)){
    target_alpha_level <- 0.2
  }
  df_final = data.frame(alpha = numeric(), 
                        type = character(), 
                        power = numeric(), 
                        time = numeric(), 
                        seed = numeric())
  start = Sys.time()
  set.seed(24062024)
  out <- adapt_gmm(pvals=pvals, x= x,
                   beta_formulas = spline_formulas,
                   nclasses=nclasses, model_type="neural",
                   alphas = alphas,target_alpha_level = target_alpha_level,intercept_model = F)
  end = Sys.time()
  for (i in 1:length(alphas)) {
    alpha = alphas[i]
    df_final = rbind(df_final, list(alpha = alpha, type = 'AdaPT_GMM', power = out$nrejs[i], time = difftime(end, start, units='mins'), seed = 24062024))
  }
  write.csv(df_final, paste('results/', file,'_adapt_gmm.csv', sep = ''))
  
}


# [1] "File: airway"
# [1] "splines::ns(x1, df = 2)" "splines::ns(x1, df = 4)"
# [3] "splines::ns(x1, df = 6)" "splines::ns(x1, df = 8)"
# Model selection starting. Shrink the set of candidate models if it is too time-consuming.
# |==================================================| 100%
# alpha = 0.2: FDPhat 0.2, Number of Rej. 7468
# alpha = 0.1: FDPhat 0.1, Number of Rej. 5731
# alpha = 0.05: FDPhat 0.05, Number of Rej. 4683
# alpha = 0.01: FDPhat 0.0099, Number of Rej. 3283
# Complete.
# $testing
# [1] "one_sided"
# 
# $model_type
# [1] "neural"
# 
# $rendpoint
# NULL
# 
# $lendpoint
# NULL
# 
# $alpha_m
# [1] 0.3
# 
# $zeta
# [1] 2
# 
# $lambda
# [1] 0.3
# 
# $niter_fit
# [1] 5
# 
# $niter_ms
# [1] 10
# 
# $nfits
# [1] 20
# 
# $masking_shape
# [1] "tent"
# 
# $p_to_z
# function (p)
#   qnorm(p, lower.tail = FALSE)
# <bytecode: 0x555e89f2e580>
#   <environment: 0x555e89f30740>
#   
#   $z_to_p
# function (z)
#   pnorm(z, lower.tail = FALSE)
# <bytecode: 0x555e89f2e740>
#   <environment: 0x555e89f30740>
#   
#   $beta_formula
# class ~ splines::ns(x1, df = 6)
# <environment: 0x555e89824220>
#   
#   $nclasses
# [1] 6
# 
# $all_a
# [1] "s" "b"
# 
# $n
# [1] 33468
# 
# $symmetric_modeling
# [1] FALSE
# 
# $base_prob
# function (z, se)
# {
#   return(dnorm(z, 0, se))
# }
# <bytecode: 0x555e89f33ee0>
#   <environment: namespace:AdaPTGMM>
#   
#   attr(,"class")
# [1] "File: adipose_subcutaneous"
# [1] "splines::ns(x1, df = 3) + splines::ns(x2, df = 3) + splines::ns(x3, df = 3) + splines::ns(x4, df = 1)"
# p
# 4.521527e-19
# Model selection starting. Shrink the set of candidate models if it is too time-consuming.
# |==================================================| 100%
# alpha = 0.2: FDPhat 0.2, Number of Rej. 3318
# alpha = 0.1: FDPhat 0.0999, Number of Rej. 2358
# alpha = 0.05: FDPhat 0.0498, Number of Rej. 1728
# alpha = 0.01: FDPhat 0.0098, Number of Rej. 1279
# Complete.
# $testing
# [1] "one_sided"
# 
# $model_type
# [1] "neural"
# 
# $rendpoint
# NULL
# 
# $lendpoint
# NULL
# 
# $alpha_m
# [1] 0.3
# 
# $zeta
# [1] 2
# 
# $lambda
# [1] 0.3
# 
# $niter_fit
# [1] 5
# 
# $niter_ms
# [1] 10
# 
# $nfits
# [1] 20
# 
# $masking_shape
# [1] "tent"
# 
# $p_to_z
# function (p)
#   qnorm(p, lower.tail = FALSE)
# <bytecode: 0x555e89f2e580>
#   <environment: 0x555e9103b518>
#   
#   $z_to_p
# function (z)
#   pnorm(z, lower.tail = FALSE)
# <bytecode: 0x555e89f2e740>
#   <environment: 0x555e9103b518>
#   
#   $beta_formula
# class ~ splines::ns(x1, df = 3) + splines::ns(x2, df = 3) + splines::ns(x3,
#                                                                         df = 3) + splines::ns(x4, df = 1)
# <environment: 0x555e8dafc8e0>
#   
#   $nclasses
# [1] 4
# 
# $all_a
# [1] "s" "b"
# 
# $n
# [1] 300000
# 
# $symmetric_modeling
# [1] FALSE
# 
# $base_prob
# function (z, se)
# {
#   return(dnorm(z, 0, se))
# }
# <bytecode: 0x555e89f33ee0>
#   <environment: namespace:AdaPTGMM>
#   
#   attr(,"class")
# [1] "args"
# [1] "File: adipose_visceral_omentum"
# [1] "splines::ns(x1, df = 3) + splines::ns(x2, df = 3) + splines::ns(x3, df = 3) + splines::ns(x4, df = 1)"
# p
# 3.574731e-18
# Model selection starting. Shrink the set of candidate models if it is too time-consuming.
# |==================================================| 100%
# alpha = 0.2: FDPhat 0.2, Number of Rej. 2673
# alpha = 0.1: FDPhat 0.0999, Number of Rej. 1676
# alpha = 0.05: FDPhat 0.0499, Number of Rej. 1213
# alpha = 0.01: FDPhat 0.0099, Number of Rej. 707
# Complete.
# $testing
# [1] "one_sided"
# 
# $model_type
# [1] "neural"
# 
# $rendpoint
# NULL
# 
# $lendpoint
# NULL
# 
# $alpha_m
# [1] 0.3
# 
# $zeta
# [1] 2
# 
# $lambda
# [1] 0.3
# 
# $niter_fit
# [1] 5
# 
# $niter_ms
# [1] 10
# 
# $nfits
# [1] 20
# 
# $masking_shape
# [1] "tent"
# 
# $p_to_z
# function (p)
#   qnorm(p, lower.tail = FALSE)
# <bytecode: 0x555e89f2e580>
#   <environment: 0x555ecbade780>
#   
#   $z_to_p
# function (z)
#   pnorm(z, lower.tail = FALSE)
# <bytecode: 0x555e89f2e740>
#   <environment: 0x555ecbade780>
#   
#   $beta_formula
# class ~ splines::ns(x1, df = 3) + splines::ns(x2, df = 3) + splines::ns(x3,
#                                                                         df = 3) + splines::ns(x4, df = 1)
# <environment: 0x555e8fe4fe50>
#   
#   $nclasses
# [1] 4
# 
# $all_a
# [1] "s" "b"
# 
# $n
# [1] 300000
# 
# $symmetric_modeling
# [1] FALSE
# 
# $base_prob
# function (z, se)
# {
#   return(dnorm(z, 0, se))
# }
# <bytecode: 0x555e89f33ee0>
#   <environment: namespace:AdaPTGMM>
# [1] "File: bottomly"
# Model selection starting. Shrink the set of candidate models if it is too time-consuming.
# |==================================================| 100%
# alpha = 0.2: FDPhat 0.2, Number of Rej. 3293
# alpha = 0.1: FDPhat 0.0999, Number of Rej. 2142
# alpha = 0.05: FDPhat 0.0499, Number of Rej. 1523
# Complete.
# $testing
# [1] "one_sided"
# 
# $model_type
# [1] "neural"
# 
# $rendpoint
# NULL
# 
# $lendpoint
# NULL
# 
# $alpha_m
# [1] 0.3
# 
# $zeta
# [1] 2
# 
# $lambda
# [1] 0.3
# 
# $niter_fit
# [1] 5
# 
# $niter_ms
# [1] 10
# 
# $nfits
# [1] 20
# 
# $masking_shape
# [1] "tent"
# 
# $p_to_z
# function (p)
#   qnorm(p, lower.tail = FALSE)
# <bytecode: 0x56369a6eaed8>
#   <environment: 0x56369acc4650>
#   
#   $z_to_p
# function (z)
#   pnorm(z, lower.tail = FALSE)
# <bytecode: 0x56369a6eab90>
#   <environment: 0x56369acc4650>
#   
#   $beta_formula
# class ~ splines::ns(x1, df = 4)
# <environment: 0x56369a862c98>
#   
#   $nclasses
# [1] 6
# 
# $all_a
# [1] "s" "b"
# 
# $n
# [1] 13931
# 
# $symmetric_modeling
# [1] FALSE
# 
# $base_prob
# function (z, se)
# {
#   return(dnorm(z, 0, se))
# }
# <bytecode: 0x56369a6e9bf8>
#   <environment: namespace:AdaPTGMM>
#   
#   attr(,"class")
# [1] "args"
# [1] "File: pasilla"
# [1] "splines::ns(x1, df = 2)" "splines::ns(x1, df = 4)"
# [3] "splines::ns(x1, df = 6)" "splines::ns(x1, df = 8)"
# Model selection starting. Shrink the set of candidate models if it is too time-consuming.
# |==================================================| 100%
# alpha = 0.2: FDPhat 0.2, Number of Rej. 1205
# alpha = 0.1: FDPhat 0.0995, Number of Rej. 879
# alpha = 0.05: FDPhat 0.05, Number of Rej. 660
# alpha = 0.01: FDPhat 0.0095, Number of Rej. 420
# Complete.
# $testing
# [1] "one_sided"
# 
# $model_type
# [1] "neural"
# 
# $rendpoint
# NULL
# 
# $lendpoint
# NULL
# 
# $alpha_m
# [1] 0.3
# 
# $zeta
# [1] 2
# 
# $lambda
# [1] 0.3
# 
# $niter_fit
# [1] 5
# 
# $niter_ms
# [1] 10
# 
# $nfits
# [1] 20
# 
# $masking_shape
# [1] "tent"
# 
# $p_to_z
# function (p)
#   qnorm(p, lower.tail = FALSE)
# <bytecode: 0x5598975a8740>
#   <environment: 0x5598975aa900>
#   
#   $z_to_p
# function (z)
#   pnorm(z, lower.tail = FALSE)
# <bytecode: 0x5598975a8900>
#   <environment: 0x5598975aa900>
#   
#   $beta_formula
# class ~ splines::ns(x1, df = 6)
# <environment: 0x559896e9e530>
#   
#   $nclasses
# [1] 6
# 
# $all_a
# [1] "s" "b"
# 
# $n
# [1] 11831
# 
# $symmetric_modeling
# [1] FALSE
# 
# $base_prob
# function (z, se)
# {
#   return(dnorm(z, 0, se))
# }
# <bytecode: 0x5598975ae0a0>
#   <environment: namespace:AdaPTGMM>
#   
#   attr(,"class")
# [1] "args"
# [1] "File: microbiome_enigma_ph"
# Model selection starting. Shrink the set of candidate models if it is too time-consuming.
# |==================================================| 100%
# alpha = 0.2: FDPhat 0.1982, Number of Rej. 169
# alpha = 0.1: FDPhat 0.0969, Number of Rej. 98
# alpha = 0.05: FDPhat 0.0479, Number of Rej. 73
# Complete.

#me
#[1] "File: microbiome_enigma_ph"
#Model selection starting. Shrink the set of candidate models or reduce niter_fit/niter_ms/nfits if it is too time-consuming.
#|==================================================| 100%
#alpha = 0.2: FDPhat 0.2, Number of Rej. 140
#alpha = 0.1: FDPhat 0.1, Number of Rej. 95
#alpha = 0.05: FDPhat 0.0424, Number of Rej. 59
#alpha = 0.01: FDPhat 0.02, Number of Rej. 0

# $testing
# [1] "one_sided"
# 
# $model_type
# [1] "neural"
# 
# $rendpoint
# NULL
# 
# $lendpoint
# NULL
# 
# $alpha_m
# [1] 0.3
# 
# $zeta
# [1] 2
# 
# $lambda
# [1] 0.3
# 
# $niter_fit
# [1] 5
# 
# $niter_ms
# [1] 10
# 
# $nfits
# [1] 20
# 
# $masking_shape
# [1] "tent"
# 
# $p_to_z
# function (p)
#   qnorm(p, lower.tail = FALSE)
# <bytecode: 0x56369a6eaed8>
#   <environment: 0x56369e6c57d8>
#   
#   $z_to_p
# function (z)
#   pnorm(z, lower.tail = FALSE)
# <bytecode: 0x56369a6eab90>
#   <environment: 0x56369e6c57d8>
#   
#   $beta_formula
# class ~ splines::ns(x1, df = 3) + splines::ns(x2, df = 3)
# <environment: 0x56369e1d0f48>
#   
#   $nclasses
# [1] 6
# 
# $all_a
# [1] "s" "b"
# 
# $n
# [1] 4007
# 
# $symmetric_modeling
# [1] FALSE
# 
# $base_prob
# function (z, se)
# {
#   return(dnorm(z, 0, se))
# }
# <bytecode: 0x56369a6e9bf8>
#   <environment: namespace:AdaPTGMM>
#   
#   attr(,"class")
# [1] "args"
# [1] "File: microbiome_enigma_al"
# Model selection starting. Shrink the set of candidate models if it is too time-consuming.
# |==================================================| 100%
# alpha = 0.2: FDPhat 0.1998, Number of Rej. 478
# alpha = 0.1: FDPhat 0.0991, Number of Rej. 217
# alpha = 0.05: FDPhat 0.0491, Number of Rej. 112
# Complete.

#[1] "File: microbiome_enigma_al"
#Model selection starting. Shrink the set of candidate models or reduce niter_fit/niter_ms/nfits if it is too time-consuming.
#|==================================================| 100%
#alpha = 0.2: FDPhat 0.1998, Number of Rej. 488
#alpha = 0.1: FDPhat 0.1, Number of Rej. 220
#alpha = 0.05: FDPhat 0.0474, Number of Rej. 116
#alpha = 0.01: FDPhat 0.0111, Number of Rej. 0

# $testing
# [1] "one_sided"
# 
# $model_type
# [1] "neural"
# 
# $rendpoint
# NULL
# 
# $lendpoint
# NULL
# 
# $alpha_m
# [1] 0.3
# 
# $zeta
# [1] 2
# 
# $lambda
# [1] 0.3
# 
# $niter_fit
# [1] 5
# 
# $niter_ms
# [1] 10
# 
# $nfits
# [1] 20
# 
# $masking_shape
# [1] "tent"
# 
# $p_to_z
# function (p)
#   qnorm(p, lower.tail = FALSE)
# <bytecode: 0x56369a6eaed8>
#   <environment: 0x56369f33b4d8>
#   
#   $z_to_p
# function (z)
#   pnorm(z, lower.tail = FALSE)
# <bytecode: 0x56369a6eab90>
#   <environment: 0x56369f33b4d8>
#   
#   $beta_formula
# class ~ splines::ns(x1, df = 3) + splines::ns(x2, df = 3)
# <environment: 0x56369f4961f0>
#   
#   $nclasses
# [1] 4
# 
# $all_a
# [1] "s" "b"
# 
# $n
# [1] 3993
# 
# $symmetric_modeling
# [1] FALSE
# 
# $base_prob
# function (z, se)
# {
#   return(dnorm(z, 0, se))
# }
# <bytecode: 0x56369a6e9bf8>
#   <environment: namespace:AdaPTGMM>
#   
#   attr(,"class")
# [1] "args"
# [1] "File: proteomics"
# Model selection starting. Shrink the set of candidate models if it is too time-consuming.
# |==================================================| 100%
# alpha = 0.2: FDPhat 0.2, Number of Rej. 615
# alpha = 0.1: FDPhat 0.0995, Number of Rej. 387
# alpha = 0.05: FDPhat 0.0486, Number of Rej. 257
# Complete.
# $testing
# [1] "one_sided"
# 
# $model_type
# [1] "neural"
# 
# $rendpoint
# NULL
# 
# $lendpoint
# NULL
# 
# $alpha_m
# [1] 0.3
# 
# $zeta
# [1] 2
# 
# $lambda
# [1] 0.3
# 
# $niter_fit
# [1] 5
# 
# $niter_ms
# [1] 10
# 
# $nfits
# [1] 20
# 
# $masking_shape
# [1] "tent"
# 
# $p_to_z
# function (p)
#   qnorm(p, lower.tail = FALSE)
# <bytecode: 0x56369a6eaed8>
#   <environment: 0x56369fe7ee58>
#   
#   $z_to_p
# function (z)
#   pnorm(z, lower.tail = FALSE)
# <bytecode: 0x56369a6eab90>
#   <environment: 0x56369fe7ee58>
#   
#   $beta_formula
# class ~ splines::ns(x1, df = 2)
# <environment: 0x56369ea783e8>
#   
#   $nclasses
# [1] 4
# 
# $all_a
# [1] "s" "b"
# 
# $n
# [1] 2666
# 
# $symmetric_modeling
# [1] FALSE
# 
# $base_prob
# function (z, se)
# {
#   return(dnorm(z, 0, se))
# }
# <bytecode: 0x56369a6e9bf8>
#   <environment: namespace:AdaPTGMM>
#   
#   attr(,"class")
# [1] "args"
# [1] "File: fmri_auditory"
# Model selection starting. Shrink the set of candidate models if it is too time-consuming.
# |==================================================| 100%
# alpha = 0.2: FDPhat 0.1998, Number of Rej. 1652
# alpha = 0.1: FDPhat 0.1, Number of Rej. 1125
# alpha = 0.05: FDPhat 0.05, Number of Rej. 800
# Complete.
# $testing
# [1] "one_sided"
# 
# $model_type
# [1] "neural"
# 
# $rendpoint
# NULL
# 
# $lendpoint
# NULL
# 
# $alpha_m
# [1] 0.3
# 
# $zeta
# [1] 2
# 
# $lambda
# [1] 0.3
# 
# $niter_fit
# [1] 5
# 
# $niter_ms
# [1] 10
# 
# $nfits
# [1] 20
# 
# $masking_shape
# [1] "tent"
# 
# $p_to_z
# function (p)
#   qnorm(p, lower.tail = FALSE)
# <bytecode: 0x56369a6eaed8>
#   <environment: 0x56369cffb5d0>
#   
#   $z_to_p
# function (z)
#   pnorm(z, lower.tail = FALSE)
# <bytecode: 0x56369a6eab90>
#   <environment: 0x56369cffb5d0>
#   
#   $beta_formula
# class ~ splines::ns(x1, df = 3) + splines::ns(x2, df = 3) + splines::ns(x3,
#                                                                         df = 3) + splines::ns(x4, df = 3)
# <environment: 0x56369daa2dc8>
#   
#   $nclasses
# [1] 6
# 
# $all_a
# [1] "s" "b"
# 
# $n
# [1] 56682
# 
# $symmetric_modeling
# [1] FALSE
# 
# $base_prob
# function (z, se)
# {
#   return(dnorm(z, 0, se))
# }
# <bytecode: 0x56369a6e9bf8>
#   <environment: namespace:AdaPTGMM>
#   
#   attr(,"class")
# [1] "args"
# [1] "File: fmri_imagination"
# Model selection starting. Shrink the set of candidate models if it is too time-consuming.
# |==================================================| 100%
# alpha = 0.2: FDPhat 0.2, Number of Rej. 4055
# alpha = 0.1: FDPhat 0.0999, Number of Rej. 2982
# alpha = 0.05: FDPhat 0.0499, Number of Rej. 2324
# Complete.
# $testing
# [1] "one_sided"
# 
# $model_type
# [1] "neural"
# 
# $rendpoint
# NULL
# 
# $lendpoint
# NULL
# 
# $alpha_m
# [1] 0.3
# 
# $zeta
# [1] 2
# 
# $lambda
# [1] 0.3
# 
# $niter_fit
# [1] 5
# 
# $niter_ms
# [1] 10
# 
# $nfits
# [1] 20
# 
# $masking_shape
# [1] "tent"
# 
# $p_to_z
# function (p)
#   qnorm(p, lower.tail = FALSE)
# <bytecode: 0x56369a6eaed8>
#   <environment: 0x56369d0fbba8>
#   
#   $z_to_p
# function (z)
#   pnorm(z, lower.tail = FALSE)
# <bytecode: 0x56369a6eab90>
#   <environment: 0x56369d0fbba8>
#   
#   $beta_formula
# class ~ splines::ns(x1, df = 5) + splines::ns(x2, df = 5) + splines::ns(x3,
#                                                                         df = 5) + splines::ns(x4, df = 5)
# <environment: 0x56369d0bdad0>
#   
#   $nclasses
# [1] 4
# 
# $all_a
# [1] "s" "b"
# 
# $n
# [1] 15648
# 
# $symmetric_modeling
# [1] FALSE
# 
# $base_prob
# function (z, se)
# {
#   return(dnorm(z, 0, se))
# }
# <bytecode: 0x56369a6e9bf8>
#   <environment: namespace:AdaPTGMM>
#   
#   attr(,"class")
# [1] "args"
