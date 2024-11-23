set.seed(14052024)
repeat_times <- 100

simul1 <- function(x, mu, H0, repeat_times, out){
  n <- length(mu)
  for (i in 1:repeat_times){
    z <- rnorm(n) + mu
    pvals <- 1 - pnorm(z)
    df <- data.frame(pvals, x, z, mu, H0)
    write.csv(df, paste('data1/', out, '_', i, '.csv', sep = ''))
    print(paste0(i, "-th step finishes!"))
  }
  return(df)
}

####### Generate x
n <- 2500
x1 <- seq(-100, 100, length.out = 50)
x2 <- seq(-100, 100, length.out = 50)
x <- expand.grid(x1, x2)
colnames(x) <- c("x1", "x2")


## Case 1: a circle in the center
H0 <- apply(x, 1, function(coord){sum(coord^2) < 900})
mu <- ifelse(H0, 2, 0)

simul1(x, mu, H0, repeat_times,
                  'sim1')

## Case 2: a circle in the corner
H0 <- apply(x, 1, function(coord){sum((coord - 65)^2) < 900})
mu <- ifelse(H0, 2, 0)

simul1(x, mu, H0, repeat_times,
       'sim2')

## Case 3: a thin ellipsoid
shape_fun <- function(coord){
  transform_coord <- c(coord[1] + coord[2], coord[2] - coord[1])/sqrt(2)
  transform_coord[1]^2 / 100^2 + transform_coord[2]^2 / 15^2 < 1
}
H0 <- apply(x, 1, shape_fun)
mu <- ifelse(H0, 2, 0)

simul1(x, mu, H0, repeat_times,
       'sim3')
