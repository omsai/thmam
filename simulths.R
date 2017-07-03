source("/Users/ChaoranHu/Desktop/thmam/fitMovResHun.R")
source("/Users/ChaoranHu/Desktop/thmam/rMovResHun.R")

lam0 <- 4
lam1 <- .5
lam2 <- .1
p <- .8
sigma <- 1

grid <- seq(0, 8, length.out = 10)

data <- rMovResHun(grid, lam0, lam1, lam2, sigma, p, "m")

fitMovResHun2(data, c(lam0, lam1, lam2, sigma, p))
