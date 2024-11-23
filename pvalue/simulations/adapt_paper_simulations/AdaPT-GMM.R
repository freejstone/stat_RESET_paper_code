library("AdaPTGMM")


args = commandArgs(trailingOnly = TRUE)
path = args[1]
file = args[2]
type = args[3]
out = args[4]
seed = as.numeric(args[5])
set.seed(seed + 22052024)
alphas = seq(0.01, 0.3, 0.01)

df = read.csv(file = paste(path, file, sep = '/'))

x = df[, grepl('x', colnames(df))]
pvals = df$pvals
z = df$z
H0 = df$H0
if (type == 'gam') {
  formulas = c("mgcv::s(x1, x2)") #Best to use s() to capture interaction and to compare to AdaPT
  res = try(adapt_gmm(x = x, z = z, beta_formulas = formulas, model_type = "mgcv", alphas = alphas))
} else if (type == 'glmnet') {
  formulas = paste(paste("x", seq(100), sep = ""), collapse = '+')
  res = try(adapt_gmm(x = x, pvals = pvals, beta_formulas = formulas, model_type = "glmnet", alphas = alphas))
}

if (class(res) != "try-error"){
  power = sapply(seq(length(alphas)), function(x)sum(H0[res$rejs[[x]]]))/max(sum(H0), 1)
  FDP = sapply(seq(length(alphas)), function(x)sum(!H0[res$rejs[[x]]]))/pmax(res$nrejs, 1)
  df_out = data.frame(FDP = FDP, power =  power, type = rep('adapt-gmm', length(alphas)), alpha = alphas)
  write.csv(df_out, out)
}
