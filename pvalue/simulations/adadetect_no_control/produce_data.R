set.seed(10052025)

simul1 <- function(n = 2e4, repeat_times = 100, out = 'sim1'){
  w = c(0, 0, 1)
  a = 0.5; mu1 = 0.25; mu2 = 0.75; s1 = s2 = 0.05
  fslope = function(x) {
    return(exp(a*x)*a/(exp(a) - 1))
  } 
  fbump = function(x, mu, s) {
    normalisation = pnorm(1, mean = mu, sd = s) - pnorm(0, mean = mu, sd = s)
    return(dnorm(x, mean = mu, sd = s)/normalisation)
  }
  ph = function(x) {
    return(sum(w*c(fslope(x), fbump(x, mu1, s1), fbump(x, mu2, s2)))*0.1)
  }
  for (i in 1:repeat_times){
    x = runif(n)
    phs = sapply(x, ph)
    H0 = as.integer(runif(n) <= phs)
    zs = rep(0, n)
    zs[H0 == 1] = rnorm(n = sum(H0 == 1), mean = 3)
    means = sample(c(0, -0.5), size = sum(H0 == 0 & x < 0.25), replace = TRUE)
    zs[H0 == 0 & x < 0.25] = rnorm(n = sum(H0 == 0 & x < 0.25), mean = means)
    zs[H0 == 0 & x >= 0.25] = rnorm(n = sum(H0 == 0 & x >= 0.25))
    pvals = 1 - pnorm(zs)
    df = data.frame(pvals, x, zs, H0)
    write.csv(df, paste('data1/', out, '_', i, '.csv', sep = ''))
    print(paste0(i, "-th step finishes!"))
  }
  return(df)
}

####### Generate x
simul1()
