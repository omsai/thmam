Rcpp::sourceCpp("/Users/ChaoranHu/Desktop/thmam/thmam_hyper.cpp")

## use optim L-BFGS-B
fitMovResHun1 <- function(data, start,
                          optim.control = list(),
                          integrControl = integr.control()) {
    if (!is.matrix(data)) data <- as.matrix(data)
    dinc <- apply(data, 2, diff)
    integrControl <- unlist(integrControl)
    
    fit <- optim(par = start, fn = nllk_fwd_ths, data = dinc,
                 integrControl = integrControl,
                 method = "L-BFGS-B",
                 lower = rep(0.0001, 5),
                 upper = c(Inf, Inf, Inf, Inf, 0.9999),
                 control = optim.control)

    fit
}

## use BB::spg
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

## use minqa::bobyqa
fitMovResHun3 <- function(data, start, lower, upper,
                          optim.control = list(),
                          integrControl = integr.control()) {
    if (!is.matrix(data)) data <- as.matrix(data)
    dinc <- apply(data, 2, diff)
    integrControl <- unlist(integrControl)

    fit <- minqa::bobyqa(par = start, fn = nllk_fwd_ths, data = dinc,
                         integrControl = integrControl,
                         lower = lower,
                         upper = upper,
                         control = optim.control)

    fit
}

## use lbfgsb3::lbfgsb3
fitMovResHun4 <- function(data, start,
                          optim.control = list(),
                          integrControl = integr.control()) {
    if (!is.matrix(data)) data <- as.matrix(data)
    dinc <- apply(data, 2, diff)
    integrControl <- unlist(integrControl)

    fit <- lbfgsb3::lbfgsb3(prm = start, fn = nllk_fwd_ths, data = dinc,
                            integrControl = integrControl,
                            lower = rep(0.0001, 5),
                            upper = c(Inf, Inf, Inf, Inf, 0.9999),
                            control = optim.control)

    fit
}

## use nlopt::cobyla
fitMovResHun5 <- function(data, start, lower, upper,
                          integrControl = integr.control()) {
    if (!is.matrix(data)) data <- as.matrix(data)
    dinc <- apply(data, 2, diff)
    integrControl <- unlist(integrControl)

    fit <- nloptr::nloptr(x0 = start, eval_f = nllk_fwd_ths,
                          data = dinc,
                          integrControl = integrControl,
                          lb = lower,
                          ub = upper,
                          opts = list("algorithm"   = "NLOPT_LN_COBYLA",
                                      "print_level" = 3,
                                      "maxeval" = 0))

    fit
}

## use nlopt::bobyqa
fitMovResHun6 <- function(data, start, lower, upper,
                          integrControl = integr.control()) {
    if (!is.matrix(data)) data <- as.matrix(data)
    dinc <- apply(data, 2, diff)
    integrControl <- unlist(integrControl)

    fit <- nloptr::nloptr(x0 = start, eval_f = nllk_fwd_ths,
                          data = dinc,
                          integrControl = integrControl,
                          lb = lower,
                          ub = upper,
                          opts = list("algorithm"   = "NLOPT_LN_BOBYQA",
                                      "print_level" = 3,
                                      "maxeval" = 0))

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

