library("adaptMT")
library("mgcv")
source('../../functions/summarize_methods.R')


args = commandArgs(trailingOnly = TRUE)
path = args[1]
file = args[2]
type = args[3]
out = args[4]
seed = as.numeric(args[5])
set.seed(seed + 10052025)
alphas = seq(0.01, 0.3, 0.01)

df = read.csv(file = paste(path, file, sep = '/'))

x = as.data.frame(df[, grepl('x', colnames(df))])
names(x) = 'x'
pvals = df$pvals
H0 = (df$H0 == 1)

if (type == 'gam') {
  pi_formula = args[6]
  mu_formula = args[7]
  print(class(x))
  res = try(adapt_gam(x, pvals,
                      pi_formula = pi_formula,
                      mu_formula = mu_formula,
                      alphas = alphas))
} else if (type == 'glmnet') {
  res = try(adapt_glmnet(as.matrix(x), pvals, alphas = alphas))
}

if (class(res) != "try-error"){
  adapt_result = summary_adapt(res, H0, pvals)
  df_out = data.frame(FDP = adapt_result[, 2], power =  adapt_result[, 3], type = rep('adapt', length(alphas)), alpha = alphas)
  write.csv(df_out, out)
}
