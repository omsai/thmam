options(echo = TRUE)
args <- commandArgs(trailingOnly = TRUE)


simul1 <-  function(){
    source("fitMovResHun.R")
    source("rMovResHun.R")

    lam0 <- 4
    lam1 <- .5
    lam2 <- .1
    p <- .8
    sigma <- 25

    grid <- seq(0, 20, length.out = 2)

    data <- rMovResHun(grid, lam0, lam1, lam2, sigma, p, "m")

    fit <- fitMovResHun5(data,c(lam0, lam1, lam2, sigma, p),
                         lower = c(0.001, 0.001, 0.001, 1, 0.001),
                         upper = c(10, 10, 50, ,50, 0.999))
    fit$solution
}

if (args == "parallel") {
    library(parallel)
    library(snow)
    cl <- makeCluster(Sys.getenv()["SLURM_NTASKS"], type = "MPI")
    clusterSetupRNG(cl, seed = c(1,2,3,4,5,6))
    result <- clusterCall(cl, simul1)
    save(result, file =  "hpc.Rdata")
    stopCluster(cl)
}
