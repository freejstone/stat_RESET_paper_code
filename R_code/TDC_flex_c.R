TDC_flex_c = function(KF, obs, score = NULL, BC1 = 1, c = 1/2, lambda = 1/2) {
  if (!is.null(score)) {
    indxs = order(score, decreasing = TRUE)
    KF = KF[indxs]
    obs = obs[indxs]
  }
  nTD = cumsum(obs)
  nDD = cumsum(KF)
  fdps = (pmin(1, ((BC1 + nDD)/ nTD) * (c / (1-lambda)))) #nDD, nTD have same length
  qvals = rev(cummin(rev(fdps)))
  if (!is.null(score)) {
    qvals[indxs] = qvals
  }
  return(qvals)
}