# Note that when I originally ran the data in R, RNGkind was set to L'Euyer-CMRG.
# Hence the following line is needed to reproduce the 100 runs.
RNGkind(kind = "L'Ecuyer-CMRG")
set.seed(15052024)
repeat_times <- 100

simul2 <- function(x, mu, pi,
                   repeat_times, out){
  n <- length(mu)
  for (i in 1:repeat_times){
    H0 <- as.logical(ifelse(runif(n) < pi, 1, 0))
    z <- ifelse(H0, rexp(n, 1/mu), rexp(n, 1))
    pvals <- exp(-z)
    
    df <- data.frame(pvals, z, x, mu, H0)
    write.csv(df, paste('data2/', out, '_', i, '.csv', sep = ''))
    print(paste0(i, "-th step finishes!"))
  }
}

#### Simulation 2
m <- 100
n <- 2000

x <- matrix(runif(n * m), n, m)
pi1 <- 0.3

inv_logit <- function(x) {exp(x) / (1 + exp(x))}
beta_pi <- c(3, 3, rep(0, m-2))
beta0_pi <- uniroot(function(b){
  mean(inv_logit(x %*% beta_pi + b)) - pi1
}, c(-100, 100))$root
pi <- inv_logit(x %*% beta_pi + beta0_pi)
beta_mu <- c(2, 2, rep(0, m-2))
beta0_mu <- 0
mu <- pmax(1, x %*% beta_mu + beta0_mu)
colnames(x) <- paste('x', 1:100, sep = '')

result <- simul2(x, mu, pi, repeat_times, 'sim')
