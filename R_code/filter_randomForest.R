filter_randomForest <- function (W, z, alpha = 0.1, offset = 1, reveal_prop = 0.5, mute = TRUE) 
{
  if (is.numeric(W)) {
    W = as.vector(W)
  }else {
    stop("W is not a numeric vector")
  }
  if (is.numeric(z) == 1) {
    z = as.matrix(z)
  }else {
    stop("z is not numeric")
  }
  if (is.numeric(reveal_prop) == 0) 
    stop("reveal_prop should be a numeric.")
  if (reveal_prop > 1) 
    stop("reveal_prop should be a numeric between 0 and 1")
  if (reveal_prop < 0) 
    stop("reveal_prop should be a numeric between 0 and 1")
  p = length(W)
  if (dim(z)[1] != p) {
    if (dim(z)[2] == p) {
      z = t(z)
    }
    else {
      stop("Please check the dimensionality of the side information!")
    }
  }
  pz = dim(z)[2]
  tau = rep(0, p)
  rejs = list()
  nrejs = rep(0, length(alpha))
  ordered_alpha = sort(alpha, decreasing = TRUE)
  rej.path = c()
  W_abs = abs(W)
  W_sign = as.numeric(W > 0)
  revealed_sign = rep(1, p)
  all_id = 1:p
  tau.sel = c()
  acc = c()
  revealed_id = which(W_abs <= quantile(abs(W[W!=0]),reveal_prop))
  revealed_sign[revealed_id] = W_sign[revealed_id]
  unrevealed_id = all_id[-revealed_id]
  tau[revealed_id] = W_abs[revealed_id] + 1
  rej.path = c(rej.path, revealed_id)
  for (talpha in 1:length(alpha)) {
    fdr = ordered_alpha[talpha]
    n_unrevealed = length(unrevealed_id)
    for (i in 1:n_unrevealed) {
      fdphat = adaptiveKnockoff:::calculate.fdphat(W, tau, offset = offset)
      if (mute == "FALSE") {
        print(fdphat)
      }
      if (fdphat <= fdr) {
        break
      }
      mdl = randomForest(y = as.factor(revealed_sign), 
                         x = cbind(W_abs, z), norm.votes = TRUE, ntree = 1000)
      fitted.pval = mdl$votes[, ncol(mdl$votes)]
      fitted.pval = fitted.pval[unrevealed_id]
      predicted.sign = fitted.pval > 0.5
      acc = c(acc, sum(predicted.sign == W_sign[unrevealed_id])/length(unrevealed_id))
      ind.min = which(fitted.pval == min(fitted.pval))
      if (length(ind.min) == 1) {
        ind.reveal = ind.min
      } else {
        ind.reveal = ind.min[which.min(W_abs[ind.min])]
      }
      ind.reveal = unrevealed_id[ind.reveal]
      revealed_id = c(revealed_id, ind.reveal)
      rej.path = c(rej.path, ind.reveal)
      unrevealed_id = all_id[-revealed_id]
      revealed_sign[ind.reveal] = W_sign[ind.reveal]
      tau[ind.reveal] = W_abs[ind.reveal] + 1
      
    }
    rej = which(W >= tau)
    rejs[[talpha]] = rej
    nrejs[talpha] = length(rej)
    tau.sel = c(tau.sel, length(revealed_id))
  }
  result = list(rejs = rejs, fdphat = fdphat, nrejs = nrejs, 
                rej.path = c(rej.path, unrevealed_id), unrevealed_id = unrevealed_id, 
                tau = tau.sel, acc = acc, W = W, alpha = alpha)
  return(result)
}

