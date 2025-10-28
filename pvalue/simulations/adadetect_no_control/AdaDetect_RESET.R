source('../../../R_code/ensemble_functions.R')
source('../../../R_code/helper_functions.R')


RNGkind(kind = "L'Ecuyer-CMRG")

args = commandArgs(trailingOnly = TRUE)
path = args[1]
file = args[2]
out = args[3]
seed = as.numeric(args[4])
set.seed(seed + 22052024)
alphas = c(0.01, 0.05, 0.1, 0.2)

storey_estimator = function(pvalues, threshold) {
  return((1 + sum(pvalues >= threshold))/ (length(pvalues)*(1-threshold)))
}
get_pvals = function(Labels, score = NULL) {
  if (!is.null(score)) {
    indxs = order(score, decreasing = TRUE)
    Labels = Labels[indxs]
  }
  nDD = cumsum(Labels == -1)
  pvals = (1 + nDD)/(1 + sum(Labels == -1))
  if (!is.null(score)) {
    pvals[indxs] = pvals
  }
  return(pvals)
}
BH = function(pvalues, level) {
  indxs = order(pvalues, decreasing = FALSE)
  pvalues = pvalues[indxs]
  threshs = (1:length(pvalues))*level/length(pvalues)
  discs = rep(FALSE, length(pvalues))
  last_indx = dplyr::last(which(pvalues <= threshs))
  if (!is.na(last_indx)) {
    discs[1:last_indx] = TRUE
    discs[indxs] = discs
  }
  return(discs)
}

df = read.csv(file = paste(path, file, sep = '/'))
x = df[, grepl('x', colnames(df))]
pvals = df$pvals
H0 = df$H0
W = pmax(pvals, 1 - pvals)*((pvals < 0.5) - (pvals > 0.5))
W = abs(qnorm(pvals))*((pvals < 0.5) - (pvals > 0.5))
W[W == Inf] = max(W[W!=Inf])
W[W == -Inf] = min(W[W!=-Inf])

print('running RESET ensemble')

summary_FDP = rep(0, length(alphas))
summary_power = rep(0, length(alphas))

for (j in 1:length(alphas)) {
  alpha = alphas[j]
  res = filter_ensemble_RESET(W, x, methods = c('rf', 'nn', 'gam'), verbose = TRUE, test_alpha = alpha, 
                              seed = seed + 10052025, mult = 1, reps = 10, num_cores = 20)
  considered = (res$pseudo_Labels == 1)
  Labels = res$Labels[considered]
  score = res$score[considered]
  empirical_pvals = get_pvals(Labels, score)
  l = sum(Labels == -1)
  threshold = floor(l/2)/(l + 1)
  test_set = (Labels == 1)
  empirical_pvals = empirical_pvals[test_set]
  pi0 = storey_estimator(empirical_pvals, threshold)
  print(paste('estimated pi0:', pi0))
  discs = BH(empirical_pvals, alpha/pi0)
  false_disc <- sum(discs & (!H0[considered][test_set] == 1))
  true_disc <- sum(discs & (H0[considered][test_set] == 1))
  print(paste('real pi0:', 1 - sum(Labels == 1 & H0[considered] == 1)/length(Labels)))
  print(true_disc)
  disc <- false_disc + true_disc
  summary_FDP[j] <- false_disc/max(disc, 1)
  summary_power[j] <- true_disc/sum(H0)
}

df_out1 = data.frame(FDP = summary_FDP, power = summary_power, type = rep('reset_adadetect', length(alphas)), alpha = alphas)

write.csv(df_out1, out)
