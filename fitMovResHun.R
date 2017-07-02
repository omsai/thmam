Rcpp::sourceCpp("/Users/ChaoranHu/Desktop/thmam/thmam_hyper.cpp")

fitMovResHun <- function(data, start,
                         method = "Nelder-Mead",
                         optim.control = list(),
                         integrControl = integr.control()) {
    if (!is.matrix(data)) data <- as.matrix(data)
    dinc <- apply(data, 2, diff)
    ## objfun <- switch(likelihood,
    ##                  composite = ncllk_m1_inc,
    ##                  full = nllk_inc,
    ##                  stop("Not valid likelihood type.")
    ##                  )
    integrControl <- unlist(integrControl)
    fit <- optim(start, nllk_fwd_ths, data = dinc, method = method,
                 control = optim.control, 
                 integrControl = integrControl)
    ## ## get variance estimate
    ## varest <- matrix(NA_real_, 3, 3)
    ## estimate <- fit$par
    ## ## if (likelihood == "full") {
    ##     ## always do this without log transformation
    ##     hess <- tryCatch(numDeriv::hessian(
    ##         nllk_fwd_ths, estimate, data = dinc,
    ##         integrControl = integrControl, logtr = FALSE),
    ##                      error = function(e) e)
    ##     if (!is(hess, "error")) varest <- solve(hess)
    ##     ## not -hess because objfun is negative llk already
    ## ## }
    
    list(estimate    = fit$par,
         #varest      = varest,
         loglik      = -fit$value,
         convergence = fit$convergence,
         #likelihood  = likelihood
         )
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
