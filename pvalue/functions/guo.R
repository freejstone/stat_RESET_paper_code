#Guo and Romano: A Generalized Sidak-Holm Procedure
#kth order statistic of s uniform RVs ~ beta(k, s + 1 - k)

qH = function(x, k, s) {
  qbeta(x, k, s + 1 - k)
}

fdp_sdp_guo = function(pvals, alpha, gamma) {
  s = length(pvals)
  ks = floor(alpha*(1:s)) + 1
  threshs = sapply(1:s, function(i) qH(gamma, ks[i], s - i + ks[i]))
  indxs = order(pvals)
  pvals = pvals[indxs]
  which_less = as.logical(cumprod(pvals <= threshs))
  which_less[indxs] = which_less
  return(list(which_less = which_less, threshs = threshs))
}
