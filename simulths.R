source("/Users/ChaoranHu/Desktop/thmam/fitMovResHun.R")
source("/Users/ChaoranHu/Desktop/thmam/rMovResHun.R")

lam0 <- 4
lam1 <- .5
lam2 <- .1
p <- .8
sigma <- 25

grid <- seq(0, 200, length.out = 11)

data <- rMovResHun(grid, lam0, lam1, lam2, sigma, p, "m")

fitMovResHun3(data, c(lam0, lam1, lam2, sigma, p), lower = c(3,.4,.05,20,.5), upper = c(6,.6,.2,27,.9), optim.control = list(iprint = 4))