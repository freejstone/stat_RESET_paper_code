library(adaptMT)
library(stepdownfdp)
source('../functions/guo.R')
source('../../competition/R_code/ensemble_functions.R')
data(estrogen)

set.seed(14072024)
alphas = c(0.01, 0.05, 0.1, 0.2)

out = data.frame(alpha = numeric(), 
                 type = character(), 
                 power = numeric(), 
                 time = numeric(),
                 conf = numeric(),
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
        out = rbind(out, list(alpha = alpha, type = 'reset_ensemble', power = power, time = total.time, conf = conf, seed = seed, order = i))
      }
      
      for (conf in c(0.1, 0.2, 0.5)) {
        start.time = Sys.time()
        res = fdp_sd(cbind(abs(W)[W!=0], sign(W)[W!=0]), alpha = alpha, conf = conf, c = 1/3, lambda = 1/3, procedure = 'coinflip')
        end.time = Sys.time()
        total.time = difftime(end.time, start.time, units='mins')
        power = length(res$discoveries_ind)
        out = rbind(out, list(alpha = alpha, type = 'FDP-SD', power = power, time = total.time, conf = conf, seed = seed, order = i))
      }
      
      for (conf in c(0.1, 0.2, 0.5)) {
        start.time = Sys.time()
        res = fdp_sdp_guo(pvals1, alpha, conf)
        end.time = Sys.time()
        total.time = difftime(end.time, start.time, units='mins')
        power = sum(res$which_less)
        out = rbind(out, list(alpha = alpha, type = 'GR-SD', power = power, time = total.time, conf = conf, seed = seed, order = i))
      }
      
    }
  }
  
  write.csv(out, paste('results/', 'estrogen_fdp_reset_ensemble.csv', sep = ''))
}