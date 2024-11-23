library(zap)
library(splines)
library(adaptMT)
data(estrogen)

set.seed(14072024)
alphas = seq(0.01, 0.3, 0.01)

df_final = data.frame(alpha = numeric(), 
                      type = character(), 
                      power = numeric(), 
                      time = numeric(), 
                      seed = numeric(),
                      order = character())

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

  fs = 'ns(x, df=6)'

  z = qnorm(pvals/2)*sign(rnorm(length(pvals)))
  start = Sys.time()
  x = model.matrix(formula(paste("~", fs,"-1")), data = data.frame(x))
  res_zap_asymp = zap_asymp(z, as.matrix(x), alpha = 0.1) #Get mirror statistics
  end = Sys.time()
  
  mirror_zap_statistics = res_zap_asymp$mirror_stat
  zap_statistics = res_zap_asymp$stat
  res_zap_asymp = lapply(alphas, zap_thresholding, blfdr_vec = zap_statistics, blfdr_flip_vec = mirror_zap_statistics)
  for (j in 1:length(alphas)) {
    alpha = alphas[j]
    df_final = rbind(df_final, list(alpha = alpha, type = 'zap', power = length(res_zap_asymp[j][[1]]), time = difftime(end, start, units='mins'), seed = 24062024, order = i))
  }
  
}

write.csv(df_final, paste('results/', 'estrogen_zap.csv', sep = ''))