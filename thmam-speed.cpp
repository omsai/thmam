#include <Rcpp.h>
#include <coga.h>

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
NumericVector old_vp00(NumericVector vs, double t,
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
// n starts from 1
NumericVector get_F2(double s, double t,
		     double lambda1, double lambda2,
		     int n) {
  NumericVector result(n + 1);
  double pco_x = t - s;
  for (int k = 0; k < n + 1; ++k) {
    result[k] = pcoga2dim(pco_x, k + 1, n - k, lambda1, lambda2);
  }
  return result;
}
// n starts from 2
NumericVector get_F1(double s, double t, double lambda2,
		     int n, NumericVector F2) {
  double firstF1 = R::pgamma(t - s, n, 1 / lambda2, 1, 0);
  F2.insert(0, firstF1);
  return F2;
}

double sumT_p01(double p, NumericVector F1, NumericVector F2) {
  double result = 0;
  double cart = 0;
  int n = F1.size() - 1;
  for (int k = 0; k < n + 1; ++k) {
    cart = F1[k] - F2[k];
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
  NumericVector F1(2);
  F1[0] = R::pgamma(t - s, 1, 1 / lambda2, 1, 0);
  F1[1] = R::pgamma(t - s, 1, 1 / lambda1, 1, 0);
  NumericVector F2 = get_F2(s, t, lambda1, lambda2, 1);

  while (TRUE) {
    cart = p * R::dgamma(s, n + 1, 1/lambda0, 0) * sumT_p01(p, F1, F2);
    result += cart;
    if (cart == 0 && n > 50) break;
    n++;
    F1 = get_F1(s, t, lambda2, n, F2);
    F2 = get_F2(s, t, lambda1, lambda2, n);
  }

  result += p * R::dgamma(s, 1, 1/lambda0, 0) * (1 - pcoga2dim(t - s, 1, 0, lambda1, lambda2));
  return result;
}

// [[Rcpp::export]]
NumericVector old_vp01(NumericVector vs, double t,
		       double lambda0, double lambda1, double lambda2,
		       double p) {
  int n = vs.size();
  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    result[i] = new_p01(vs[i], t, lambda0, lambda1, lambda2, p);
  }
  return result;
}

