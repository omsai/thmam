options(echo = TRUE)
library(methods)                   # Silence argparser loading message
library(argparser)
parser <- arg_parser("Test animal movement model using simulated data")
parser <- add_argument(parser, "--replicates",
                       help = "Number of replicates to generate to fit the model",
                       default = 100)
parser <- add_argument(parser, "--parallel",
                       help = "Run using MPI",
                       flag = TRUE)
args <- parse_args(parser)
args

simul1 <-  function(replicate_number = NULL) {
    source("fitMovResHun.R")
    source("rMovResHun.R")

    lam0 <- 4
    lam1 <- .5
    lam2 <- .1
    p <- .8
    sigma <- 25

    grid <- seq(0, 2000, length.out = 101)

    data <- rMovResHun(grid, lam0, lam1, lam2, sigma, p, "m")

    fit <- fitMovResHun5(data,c(lam0, lam1, lam2, sigma, p),
                         lower = c(0.001, 0.001, 0.001, 1, 0.001),
                         upper = c(10, 10, 50, ,50, 0.999))
    fit$solution
}

if (args$parallel) {
    library(parallel)
    library(snow)

    nslaves <- Sys.getenv("SLURM_NTASKS",
                          unset = parallel::detectCores())
    nslaves
    cl <- makeCluster(nslaves, type = "MPI")
    cl
    clusterSetupRNG(cl, seed = c(1,2,3,4,5,6))
    result <- clusterApply(cl, 1:args$replicates, simul1)
    save(result, file =  "hpc.Rdata")
    stopCluster(cl)
} else {
    result <- replicate(args$replicates, simul1(), simplify = FALSE)
    save(result, file =  "hpc.Rdata")
}
