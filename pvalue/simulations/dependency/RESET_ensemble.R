source('../../../R_code/ensemble_functions.R')
source('../../../R_code/helper_functions.R')

RNGkind(kind = "L'Ecuyer-CMRG")

args = commandArgs(trailingOnly = TRUE)
path = args[1]
file = args[2]
out = args[3]
seed = as.numeric(args[4])
set.seed(seed + 10052025)
alphas = c(0.01, 0.05, 0.1, 0.2)

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
flag = rep(0, length(alphas))
for (j in 1:length(alphas)) {
  alpha = alphas[j]
  res = filter_ensemble_RESET(W, x, methods = c('rf', 'nn', 'gam'), verbose = TRUE, test_alpha = alpha, 
                              seed = seed + 22052024, reps = 10, num_cores = 20)
  false_disc <- sum(res$q_vals <= alpha & res$Labels == 1 & !H0)
  true_disc <- sum(res$q_vals <= alpha & res$Labels == 1 & H0)
  print(true_disc)
  disc <- false_disc + true_disc
  summary_FDP[j] <- false_disc/max(disc, 1)
  summary_power[j] <- true_disc/sum(H0)
  flag[j] <- res$flag
}
df_out1 = data.frame(FDP = summary_FDP, power = summary_power, 
                     type = rep('reset_ensemble', length(alphas)), flag = flag, alpha = alphas)

##########################################################################################################################################################
print('running RESET nn')

summary_FDP = rep(0, length(alphas))
summary_power = rep(0, length(alphas))
flag = rep(0, length(alphas))
for (j in 1:length(alphas)) {
  alpha = alphas[j]
  res = filter_ensemble_RESET(W, x, methods = c('nn'), verbose = TRUE, test_alpha = alpha, 
                              seed = seed + 22052024, reps = 10, num_cores = 20)
  false_disc <- sum(res$q_vals <= alpha & res$Labels == 1 & !H0)
  true_disc <- sum(res$q_vals <= alpha & res$Labels == 1 & H0)
  print(true_disc)
  disc <- false_disc + true_disc
  summary_FDP[j] <- false_disc/max(disc, 1)
  summary_power[j] <- true_disc/sum(H0)
  flag[j] <- res$flag
}

df_out2 = data.frame(FDP = summary_FDP, power = summary_power, 
                     type = rep('reset_nn', length(alphas)), flag = flag, alpha = alphas)

############################################################################################################################################################
print('running RESET rf')

summary_FDP = rep(0, length(alphas))
summary_power = rep(0, length(alphas))
flag = rep(0, length(alphas))
for (j in 1:length(alphas)) {
  alpha = alphas[j]
  res = filter_ensemble_RESET(W, x, methods = c('rf'), verbose = TRUE, test_alpha = alpha, 
                              seed = seed + 22052024, reps = 10, num_cores = 20)
  false_disc <- sum(res$q_vals <= alpha & res$Labels == 1 & !H0)
  true_disc <- sum(res$q_vals <= alpha & res$Labels == 1 & H0)
  print(true_disc)
  disc <- false_disc + true_disc
  summary_FDP[j] <- false_disc/max(disc, 1)
  summary_power[j] <- true_disc/sum(H0)
  flag[j] <- res$flag
}

df_out3 = data.frame(FDP = summary_FDP, power = summary_power, 
                     type = rep('reset_rf', length(alphas)), flag = flag, alpha = alphas)

############################################################################################################################################################

print('running RESET gam')

summary_FDP = rep(0, length(alphas))
summary_power = rep(0, length(alphas))
flag = rep(0, length(alphas))
for (j in 1:length(alphas)) {
  alpha = alphas[j]
  res = filter_ensemble_RESET(W, x, methods = c('gam'), verbose = TRUE, test_alpha = alpha, 
                              seed = seed + 22052024, reps = 10, num_cores = 20)
  false_disc <- sum(res$q_vals <= alpha & res$Labels == 1 & !H0)
  true_disc <- sum(res$q_vals <= alpha & res$Labels == 1 & H0)
  print(true_disc)
  disc <- false_disc + true_disc
  summary_FDP[j] <- false_disc/max(disc, 1)
  summary_power[j] <- true_disc/sum(H0)
  flag[j] <- res$flag
}

df_out4 = data.frame(FDP = summary_FDP, power = summary_power, 
                     type = rep('reset_gam', length(alphas)), flag = flag, alpha = alphas)

write.csv(rbind(df_out1, df_out2, df_out3, df_out4), out)
