// [[Rcpp::depends(RcppGSL)]]

#include <Rcpp.h>
#include <RcppGSL.h>
#include <gsl/gsl_sf_hyperg.h>

using namespace Rcpp;

// [[Rcpp::export]]
double dcoga2dim(double x, double shape1, double shape2,
		 double rate1, double rate2) {
  // transfer rate to scale
  double beta1 = 1 / rate1;
  double beta2 = 1 / rate2;
  // handle one shape is 0
  if (shape1 == 0) return R::dgamma(x, shape2, beta2, 0);
  if (shape2 == 0) return R::dgamma(x, shape1, beta1, 0);
  // determine min beta
  if (beta1 > beta2) {
    double beta_cart = beta1;
    beta1 = beta2;
    beta2 = beta_cart;
    double shape_cart = shape1;
    shape1 = shape2;
    shape2 = shape_cart;
  }

  double lgam = shape1 + shape2;
  double cart = 0;
  double result = 0;
  int r = 0;
  
  while (TRUE) {
    cart = R::choose(shape2 + r - 1, r) / exp(R::lgammafn(lgam + r));
    cart *= pow(x * (1/beta1 - 1/beta2), r);
    if (cart == R_PosInf || R_IsNaN(cart)) {
      warning("Inf or NaN happened, not converge!");
      break;
    }
    result += cart;
    if (cart == 0) break;
    r++;
  }
  result *= pow(x, lgam - 1) * exp(-x / beta1);
  result /= pow(beta1, shape1) * pow(beta2, shape2);
  if (result == R_PosInf) warning("Overflow!");
  return result;
}

// [[Rcpp::export]]
double pcoga2dim(double x, double shape1, double shape2,
		 double rate1, double rate2) {
  // transfer rate to scale
  double beta1 = 1 / rate1;
  double beta2 = 1 / rate2;
  // handle one shape is 0
  if (shape1 == 0) return R::pgamma(x, shape2, beta2, 1, 0);
  if (shape2 == 0) return R::pgamma(x, shape1, beta1, 1, 0);
  // make convergence faster

  // determine min beta
  if (beta1 > beta2) {
    double beta_cart = beta1;
    beta1 = beta2;
    beta2 = beta_cart;
    double shape_cart = shape1;
    shape1 = shape2;
    shape2 = shape_cart;
  }

  double lgam = shape1 + shape2;
  double cart = 0;
  double result = 0;
  int r = 0;

  while (TRUE) {
    cart = R::choose(shape2 + r - 1, r) * pow(1 - beta1/beta2, r);
    cart *= R::pgamma(x/beta1, lgam + r, 1, 1, 0);
    if (cart == R_PosInf || R_IsNaN(cart)) {
      warning("Inf or NaN happened, not converge!");
      break;
    }
    result += cart;
    if (cart == 0) break;
    r++;
  }
  return result * pow(beta1/beta2, shape2);
}

// best choose
// [[Rcpp::export]]
double dcoga2dim_hyper(double x, double shape1, double shape2,
		       double rate1, double rate2) {
  // transfer rate to scale
  double beta1 = 1 / rate1;
  double beta2 = 1 / rate2;
  /*
  // handle one shape is 0
  if (shape1 == 0) return R::dgamma(x, shape2, beta2, 0);
  if (shape2 == 0) return R::dgamma(x, shape1, beta1, 0);
  
  // determine min beta
  if (beta1 > beta2) {
    double beta_cart = beta1;
    beta1 = beta2;
    beta2 = beta_cart;
    double shape_cart = shape1;
    shape1 = shape2;
    shape2 = shape_cart;
  }
  */
  double lgam = shape1 + shape2;
  double parx = (1/beta1 - 1/beta2) * x;
  double result = pow(x, lgam - 1) * exp(-x / beta1);
  result *= gsl_sf_hyperg_1F1(shape2, lgam, parx);
  result /= pow(beta1, shape1) * pow(beta2, shape2);
  result /= exp(R::lgammafn(lgam));
  return result;
}


// [[Rcpp::export]]
double dcoga2dim_recur(double x, double shape1, double shape2,
		 double rate1, double rate2) {
  // transfer rate to scale
  double beta1 = 1 / rate1;
  double beta2 = 1 / rate2;
  // handle one shape is 0
  if (shape1 == 0) return R::dgamma(x, shape2, beta2, 0);
  if (shape2 == 0) return R::dgamma(x, shape1, beta1, 0);
  // determine min beta
  if (beta1 > beta2) {
    double beta_cart = beta1;
    beta1 = beta2;
    beta2 = beta_cart;
    double shape_cart = shape1;
    shape1 = shape2;
    shape2 = shape_cart;
  }

  double lgam = shape1 + shape2;
  double cart = 1 / exp(R::lgammafn(lgam));
  double result = 0;
  int r = 0;
  double difbex = (1/beta1 - 1/beta2) * x;
  
  while (TRUE) {
    if (cart == R_PosInf || R_IsNaN(cart)) {
      warning("Inf or NaN happened, not converge!");
      break;
    }
    result += cart;
    if (cart == 0) break;
    cart *= (shape2 + r) * difbex / ((r + 1) * (lgam + r));
    r++;
  }
  result *= pow(x, lgam - 1) * exp(-x / beta1);
  result /= pow(beta1, shape1) * pow(beta2, shape2);
  return result;
}


// [[Rcpp::export]]
double pcoga2dim_recur(double x, double shape1, double shape2,
		 double rate1, double rate2) {
  // transfer rate to scale
  double beta1 = 1 / rate1;
  double beta2 = 1 / rate2;
  // handle one shape is 0
  if (shape1 == 0) return R::pgamma(x, shape2, beta2, 1, 0);
  if (shape2 == 0) return R::pgamma(x, shape1, beta1, 1, 0);
  // make convergence faster

  // determine min beta
  if (beta1 > beta2) {
    double beta_cart = beta1;
    beta1 = beta2;
    beta2 = beta_cart;
    double shape_cart = shape1;
    shape1 = shape2;
    shape2 = shape_cart;
  }

  double lgam = shape1 + shape2;
  double sun = 1 - beta1 / beta2;
  double moon = exp(-x/beta1) * pow(x/beta1, lgam);
  
  
  double cartB = 1.;
  double cartD = R::pgamma(x/beta1, lgam, 1, 1, 0);
  double cartE = 1 / exp(R::lgammafn(lgam + 1));
  double cart = cartD;
  double result = 0.;
  int r = 0;

  while (TRUE) {
    if (cart == R_PosInf || R_IsNaN(cart)) {
      //warning("Inf or NaN happened, not converge!");
      break;
    }
    result += cart;
    if (cart == 0) break;
    cartB *= sun * (shape2 + r) / (r + 1);
    cartD -= moon * cartE;
    cartE *= x / (beta1 * (r + lgam + 1));
    cart = cartB * cartD;
    r++;
  }
  return result * pow(beta1/beta2, shape2);
}


// [[Rcpp::export]]
NumericVector pcoga2dim_recur_v(NumericVector x,
				double shape1, double shape2,
				double rate1, double rate2) {
  double n = x.size();
  NumericVector result(n);
  for (int i = 0; i < n; ++i) {
    result[i] = pcoga2dim_recur(x[i], shape1, shape2, rate1, rate2);
  }
  return result;
}

// This is the best version until now
// [[Rcpp::export]]
double pcoga2dim_recur_nopgamma(double x, double shape1, double shape2,
		 double rate1, double rate2) {
  // transfer rate to scale
  double beta1 = 1 / rate1;
  double beta2 = 1 / rate2;
  // handle one shape is 0
  if (shape1 == 0) return R::pgamma(x, shape2, beta2, 1, 0);
  if (shape2 == 0) return R::pgamma(x, shape1, beta1, 1, 0);
  // make convergence faster
  /*
  // determine min beta
  if (beta1 > beta2) {
    double beta_cart = beta1;
    beta1 = beta2;
    beta2 = beta_cart;
    double shape_cart = shape1;
    shape1 = shape2;
    shape2 = shape_cart;
  }
  */
  double lgam = shape1 + shape2;
  double sun = 1 - beta1 / beta2;
  
  double cartB = 1.;
  double cartD = R::pgamma(x/beta1, lgam, 1, 1, 0);
  double cart = cartD;
  double result = 0.;
  int r = 0;

  while (TRUE) {
    if (cart == R_PosInf || R_IsNaN(cart)) {
      warning("Inf or NaN happened, not converge!");
      break;
    }
    result += cart;
    if (cart == 0) break;
    cartB *= sun * (shape2 + r) / (r + 1);
    r++;
    cartD = R::pgamma(x/beta1, lgam + r, 1, 1, 0);
    cart = cartB * cartD;
  }
  return result * pow(beta1/beta2, shape2);
}


// [[Rcpp::export]]
double pcoga2dim_diff_shape (double x,
			     double shape1, double shape2,
			     double rate1, double rate2) {
  double result = pow(rate1, shape1) * pow(rate2, shape2);
  double lgam = shape1 + shape2;
  result *= pow(x, lgam);
  result /= exp(R::lgammafn(lgam + 1));
  result *= exp(-x * rate1);
  result *= gsl_sf_hyperg_1F1(shape2, lgam + 1, x * (rate1 - rate2));
  return result;
}
