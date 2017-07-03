Rcpp::sourceCpp("/Users/ChaoranHu/Desktop/thmam/thmam_hyper.cpp")

fitMovResHun1 <- function(data, start,
                          optim.control = list(),
                          integrControl = integr.control()) {
    if (!is.matrix(data)) data <- as.matrix(data)
    dinc <- apply(data, 2, diff)
    integrControl <- unlist(integrControl)

    foo <- function (p) {nllk_fwd_ths(c(4,.5,.1,1,p), dinc, integrControl)}
    
    fit <- optim(start, foo, method = "L-BFGS-B",
                 lower = 0.5, upper = 0.9,
                 control = optim.control)
    
    list(estimate    = fit$par,
         loglik      = -fit$value,
         convergence = fit$convergence,
         counts = fit$counts
         )
}

fitMovResHun2 <- function(data, start,
                          optim.control = list(),
                          integrControl = integr.control()) {
    if (!is.matrix(data)) data <- as.matrix(data)
    dinc <- apply(data, 2, diff)
    integrControl <- unlist(integrControl)

    fit <- BB::spg(par = start, fn = nllk_fwd_ths, data = dinc,
                   integrControl = integrControl,
                   lower = rep(0.0001, 5),
                   upper = c(Inf, Inf, Inf, Inf, 0.9999),
                   control = optim.control)

    fit
}











integr.control <- function(rel.tol = .Machine$double.eps^.25,
                           abs.tol = rel.tol, subdivisions = 100L) {
    if (!is.numeric(rel.tol) || rel.tol <= 0) 
        stop("value of 'rel.tol' must be > 0")
    if (!is.numeric(abs.tol) || abs.tol <= 0) 
        stop("value of 'abs.tol' must be > 0")
    if (!is.numeric(subdivisions) || subdivisions <= 0) 
        stop("maximum number of subintervals must be > 0")
    list(rel.tol = rel.tol, abs.tol = abs.tol, subdivisions = subdivisions)
}

