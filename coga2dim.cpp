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

// [[Rcpp::export]]
double dcoga2dim_hyper(double x, double shape1, double shape2,
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
  double parx = (1/beta1 - 1/beta2) * x;
  double result = pow(x, lgam - 1) * exp(-x / beta1);
  result /= pow(beta1, shape1) * pow(beta2, shape2) * exp(R::lgammafn(lgam));
  result *= gsl_sf_hyperg_1F1(shape2, lgam, parx);

  return result;
}
