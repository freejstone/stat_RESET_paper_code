library(zap)
RNGkind("L'Ecuyer-CMRG")
files = list.files('data', pattern = 'sim1')
alphas = seq(0.01, 0.3, 0.01)
set.seed(10052025)

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

for (file in files) {
  print(file)
  df = read.csv(paste('data', file, sep = '/'))
  x = df[, grepl('x', colnames(df))]
  x = cbind(splines::ns(x, df = 6)) #copying the spline basis expansion from the real data set applications
  z = df$z
  H0 = df$H0
  summary_FDP = rep(0, length(alphas))
  summary_power = rep(0, length(alphas))
  res_zap_asymp = zap_asymp(z, as.matrix(x), alpha = Inf) #Get mirror statistics
  mirror_zap_statistics = res_zap_asymp$mirror_stat
  zap_statistics = res_zap_asymp$stat
  res_zap_asymp = lapply(alphas, zap_thresholding, blfdr_vec = zap_statistics, blfdr_flip_vec = mirror_zap_statistics)
  for (i in 1:length(alphas)) {
    disc_indices = res_zap_asymp[i][[1]]
    summary_power[i] = sum(H0[disc_indices])/max(sum(H0), 1)
    summary_FDP[i] = sum(!H0[disc_indices])/max(length(disc_indices), 1)
  }
  df_out = data.frame(FDP = summary_FDP, power = summary_power, type = rep('zap', length(alphas)), alpha = alphas)
  write.csv(df_out, paste('result_data/zap_', file, sep = ''))
}

files = list.files('data', pattern = 'sim2')
alphas = seq(0.01, 0.3, 0.01)
set.seed(10052025)
for (file in files) {
  print(file)
  df = read.csv(paste('data', file, sep = '/'))
  x = df[, grepl('x', colnames(df))]
  x = cbind(splines::ns(x, df = 6)) #copying the spline basis expansion from the real data set applications
  z = df$z
  H0 = df$H0
  summary_FDP = rep(0, length(alphas))
  summary_power = rep(0, length(alphas))
  res_zap_asymp = zap_asymp(z, as.matrix(x), alpha = Inf) #Get mirror statistics
  mirror_zap_statistics = res_zap_asymp$mirror_stat
  zap_statistics = res_zap_asymp$stat
  res_zap_asymp = lapply(alphas, zap_thresholding, blfdr_vec = zap_statistics, blfdr_flip_vec = mirror_zap_statistics)
  for (i in 1:length(alphas)) {
    disc_indices = res_zap_asymp[i][[1]]
    summary_power[i] = sum(H0[disc_indices])/max(sum(H0), 1)
    summary_FDP[i] = sum(!H0[disc_indices])/max(length(disc_indices), 1)
  }
  df_out = data.frame(FDP = summary_FDP, power = summary_power, type = rep('zap', length(alphas)), alpha = alphas)
  write.csv(df_out, paste('result_data/zap_', file, sep = ''))
}
