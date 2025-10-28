library(adaptMT)
source('../../R_code/ensemble_functions.R')
source('../../R_code/helper_functions.R')
data(estrogen)

set.seed(14072024)
alphas = c(0.01, 0.05, 0.1, 0.2)

out = data.frame(alpha = numeric(), 
                 type = character(), 
                 power = numeric(), 
                 time = numeric(), 
                 flag = logical(),
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
  
  W = abs(qnorm(pvals))*((pvals < 0.5) - (pvals > 0.5))
  W[W == Inf] = max(W[W!=Inf])
  W[W == -Inf] = min(W[W!=-Inf])
  
  mult = 1
  reps = 10
  
  for (seed in 22072024:(22072024 + 9)) {
    for (alpha in alphas) {
      start = Sys.time()
      res = filter_ensemble_RESET(W, x,
                                  verbose = TRUE, test_alpha = alpha, seed = seed,
                                  mult = mult, reps = reps, num_cores = 20)
      end = Sys.time()
      out = rbind(out, list(alpha = alpha, type = 'reset_ensemble', 
                            power = sum(res$q_vals <= alpha & res$Labels == 1), 
                            time = difftime(end, start, units='mins'), 
                            flag = res$flag, seed = seed, order = i))
    }
  }
  
  for (seed in 22072024:(22072024 + 9)) {
    for (alpha in alphas) {
      start = Sys.time()
      res = filter_ensemble_RESET(W, x, methods = 'gam',
                                  verbose = TRUE, test_alpha = alpha, seed = seed,
                                  mult = mult, reps = reps, num_cores = 20)
      end = Sys.time()
      out = rbind(out, list(alpha = alpha, type = 'reset_gam', 
                            power = sum(res$q_vals <= alpha & res$Labels == 1), 
                            time = difftime(end, start, units='mins'), 
                            flag = res$flag, seed = seed, order = i))
    }
  }
  
  count = 1
  for (seed in 22072024:(22072024 + 9)) {
    for (alpha in alphas) {
      start = Sys.time()
      res = filter_ensemble_RESET(W, x, methods = 'rf', 
                                  verbose = TRUE, test_alpha = alpha, seed = seed,
                                  mult = mult, reps = reps, num_cores = 20)
      end = Sys.time()
      out = rbind(out, list(alpha = alpha, type = 'reset_rf', 
                            power = sum(res$q_vals <= alpha & res$Labels == 1), 
                            time = difftime(end, start, units='mins'), 
                            flag = res$flag, seed = seed, order = i))
    }
  }
  
  count = 1
  for (seed in 22072024:(22072024 + 9)) {
    for (alpha in alphas) {
      start = Sys.time()
      res = filter_ensemble_RESET(W, x, methods = 'nn',
                                  verbose = TRUE, test_alpha = alpha, seed = seed,
                                  mult = mult, reps = reps, num_cores = 20)
      end = Sys.time()
      out = rbind(out, list(alpha = alpha, type = 'reset_nn', 
                            power = sum(res$q_vals <= alpha & res$Labels == 1), 
                            time = difftime(end, start, units='mins'), 
                            flag = res$flag, seed = seed, order = i))
    }
  }
}
write.csv(out, paste('results/', 'estrogen_dep_reset_ensemble.csv', sep = ''))