#include <Rcpp.h>

using namespace Rcpp;

double prod_to_m(double x, int m) {
  if (m == 0) return 1;
  double result = 1;
  for (int i = 0; i < m; ++i) {
    result *= x + i;
  }
  return result;
}

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
    cart = prod_to_m(shape2, r) * pow((1 / beta1 - 1 / beta2) * x, r);
    cart /= exp(R::lgammafn(r + 1)) * prod_to_m(lgam, r);
    if (cart == R_PosInf) break;
    if (R_IsNaN(cart)) break;
    result += cart;
    if (cart == 0) break;
    r++;
  }
  result *= pow(x, lgam - 1) * exp(-x / beta1);
  result /= (pow(beta1, shape1) * pow(beta2, shape2) * R::gammafn(lgam));
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
    cart = prod_to_m(shape2, r) * pow(1 / beta1 - 1 / beta2, r);
    cart /= exp(R::lgammafn(r + 1));
    cart *= pow(beta1, r + lgam) * R::pgamma(x / beta1, r + lgam, 1, 1, 0);
    cart *= exp(R::lgammafn(lgam));
    if (cart == R_PosInf) break;
    if (R_IsNaN(cart)) break;
    result +=  cart;
    if (cart == 0) break;
    r++;
  }
  result /= pow(beta1, shape1) * pow(beta2, shape2) * R::gammafn(lgam);
  return result;
}


/**************************************************************************
 ****************************** calculate p  ******************************
 **************************************************************************/

/****************************** calculate p00 *****************************/
double sumT_p00(double s, double t,
	    double lambda1, double lambda2,
	    double p, int n) {
  double result = 0;
  for (int k = 0; k < n + 1; ++k) {
    result += dcoga2dim(t - s, k, n - k, lambda1, lambda2) * R::choose(n, k) * pow(p, k) * pow(1 - p, n - k);
  }
  return result;
}

double new_p00(double s, double t,
	       double lambda0, double lambda1, double lambda2,
	       double p) {
  int n = 1;
  double result = 0;
  double cart = 0;

  while (TRUE) {
    cart = R::pgamma(s, n, 1/lambda0, 1, 0) - R::pgamma(s, n + 1, 1/lambda0, 1, 0);
    cart *= sumT_p00(s, t, lambda1, lambda2, p, n);
    result += cart;
    if (cart == 0 && n > 50) break;
    n++;
  }

  return result;
}

// [[Rcpp::export]]
NumericVector new_vp00(NumericVector vs, double t,
		       double lambda0, double lambda1, double lambda2,
		       double p) {
  int n = vs.size();
  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    result[i] = new_p00(vs[i], t, lambda0, lambda1, lambda2, p);
  }
  return result;
}

/****************************** calculate p01 *****************************/
// [[Rcpp::export]]
double sumT_p01(double s, double t,
		double lambda1, double lambda2,
		double p, int n) {
  double result = 0;
  double cart = 0;
  for (int k = 0; k < n + 1; ++k) {
    cart = pcoga2dim(t - s, k, n - k, lambda1, lambda2) - pcoga2dim(t - s, k + 1, n - k, lambda1, lambda2);
    cart *= R::choose(n, k) * pow(p, k) * pow(1 - p, n - k);
    result += cart;
  }
  return result;
}

double new_p01(double s, double t,
	       double lambda0, double lambda1, double lambda2,
	       double p) {
  int n = 1;
  double result = 0;
  double cart = 0;

  while (TRUE) {
    cart = sumT_p01(s, t, lambda1, lambda2, p, n);
    cart *= p * R::dgamma(s, n + 1, 1/lambda0, 0);
    return cart;
    result += cart;
    if (cart == 0 && n > 50) break;
    n++;
  }

  return result + (p * R::dgamma(s, 1, 1/lambda0, 0) * (1 - R::pgamma(t - s, 1, 1/lambda1, 1, 0)));
}

// [[Rcpp::export]]
NumericVector new_vp01(NumericVector vs, double t,
		       double lambda0, double lambda1, double lambda2,
		       double p) {
  int n = vs.size();
  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    result[i] = new_p01(vs[i], t, lambda0, lambda1, lambda2, p);
  }
  return result;
}

// pcoga2dim still have some potential problem!
