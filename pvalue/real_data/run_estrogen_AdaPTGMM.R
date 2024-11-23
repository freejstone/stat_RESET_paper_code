library(AdaPTGMM)
library(splines)

estrogen = read.csv('data_files/estrogen.csv')

set.seed(14072024)
alphas = seq(0.01, 0.3, 0.01)

df_final = data.frame(alpha = numeric(), 
                      type = character(), 
                      power = numeric(), 
                      time = numeric(), 
                      seed = numeric(),
                      order = character())
for (i in 1:2){
  if (i == 1){
    dataset <- "estrogen_high"
    x <- estrogen$ord_high
  } else {
    dataset <- "estrogen_mod"
    x <- estrogen$ord_mod
  }
  pvals <- estrogen$pvals
  
  reorder <- order(x)
  x <- x[reorder]
  pvals <- pvals[reorder]
  x <- data.frame(x = x)
  nclasses <- c(2,4,6)
  target_alpha_level <- 0.1
  
  start = Sys.time()
  out <- adapt_gmm(pvals=pvals, x= x,
                   beta_formulas = paste0("splines::ns(x, df = ",seq(2,8,2) ,")"),
                   nclasses=nclasses, model_type="neural",
                   alphas = alphas,target_alpha_level = target_alpha_level,intercept_model = F)
  end = Sys.time()
  for (j in 1:length(alphas)) {
    alpha = alphas[j]
    df_final = rbind(df_final, list(alpha = alpha, type = 'AdaPT_GMM', power = out$nrejs[j], time = difftime(end, start, units='mins'), seed = 24062024, order = i))
  }
  
}

write.csv(df_final, paste('results/', 'estrogen_adaptgmm.csv', sep = ''))