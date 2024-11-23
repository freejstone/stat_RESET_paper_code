library(adaptgMT)
library(splines)
data(estrogen)

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
  
  pi.formulas <- paste0("ns(x, df = ", 6:10, ")")
  mu.formulas <- paste0("ns(x, df = ", 6:10, ")")
  formulas <- expand.grid(pi.formulas, mu.formulas)
  pi_formulas <- as.character(formulas[, 1])
  mu_formulas <- as.character(formulas[, 1])
  
  start = Sys.time()
  out <- adapt_glm(x, pvals,
                   dist = beta_family(),
                   pi_formulas = pi_formulas,
                   mu_formulas = mu_formulas,
                   nfits = 50)
  end = Sys.time()
  for (j in 1:length(alphas)) {
    alpha = alphas[j]
    df_final = rbind(df_final, list(alpha = alpha, type = 'AdaPTg', power = out$nrejs[j], time = difftime(end, start, units='mins'), seed = 24062024, order = i))
  }
}

write.csv(df_final, paste('results/', 'estrogen_adaptg.csv', sep = ''))