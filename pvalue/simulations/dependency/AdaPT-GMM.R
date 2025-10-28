library("AdaPTGMM")


args = commandArgs(trailingOnly = TRUE)
path = args[1]
file = args[2]
out = args[3]
seed = as.numeric(args[4])
set.seed(seed + 10052025)
alphas = seq(0.01, 0.3, 0.01)

df = read.csv(file = paste(path, file, sep = '/'))

x = df[, grepl('x', colnames(df))]
x = data.frame(x = x)
pvals = df$pvals
z = df$z
H0 = df$H0
formulas = c("mgcv::s(x)") 
res = try(adapt_gmm(x = x, z = z, beta_formulas = formulas, model_type = "mgcv", alphas = alphas))

if (class(res) != "try-error"){
  power = sapply(seq(length(alphas)), function(x)sum(H0[res$rejs[[x]]]))/max(sum(H0), 1)
  FDP = sapply(seq(length(alphas)), function(x)sum(!H0[res$rejs[[x]]]))/pmax(res$nrejs, 1)
  df_out = data.frame(FDP = FDP, power =  power, type = rep('adapt-gmm', length(alphas)), alpha = alphas)
  write.csv(df_out, out)
}
