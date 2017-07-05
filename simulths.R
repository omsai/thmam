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

    grid <- seq(0, 200, length.out = 11)

    data <- rMovResHun(grid, lam0, lam1, lam2, sigma, p, "m")

    fitMovResHun3(data,
                  c(lam0, lam1, lam2, sigma, p),
                  lower = c(3,.4,.05,20,.5),
                  upper = c(6,.6,.2,27,.9),
                  optim.control = list(iprint = 4))
}

if (args == "parallel") {
    library(parallel)
    cl <- makeCluster(Sys.getenv()["SLURM_NTASKS"], type = "MPI")
    result <- clusterCall(cl, simul1)
    save(result,file =  "hpc.Rdata")
    stopCluster(cl)
}
