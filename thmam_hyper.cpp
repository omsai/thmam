// [[Rcpp::depends(RcppGSL)]]

#include <Rcpp.h>
#include <RcppGSL.h>
#include <gsl/gsl_sf_hyperg.h>
#include <R_ext/Applic.h>

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

  gsl_set_error_handler_off();
  double lgam = shape1 + shape2;
  double parx = (1/beta1 - 1/beta2) * x;
  double result = pow(x, lgam - 1) * exp(-x / beta1);
  result *= gsl_sf_hyperg_1F1(shape2, lgam, parx);
  result /= pow(beta1, shape1) * pow(beta2, shape2);
  result /= exp(R::lgammafn(lgam));
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
  // handle special situations
  if (shape1 == 0 && shape2 == 0) {
    return 1 - pcoga2dim(x, 1, 0, rate1, rate2);
  }
  if (shape1 == 0 && shape2 != 0) {
    return pcoga2dim(x, 0, shape2, rate1, rate2) - pcoga2dim(x, 1, shape2, rate1, rate2);
  }
  if (shape1 != 0 && shape2 == 0) {
    return pcoga2dim(x, shape1, 0, rate1, rate2) - pcoga2dim(x, shape1 + 1, 0, rate1, rate2);
  }

  gsl_set_error_handler_off();
  double result = pow(rate1, shape1) * pow(rate2, shape2);
  double lgam = shape1 + shape2;
  result *= pow(x, lgam);
  result /= exp(R::lgammafn(lgam + 1));
  result *= exp(-x * rate1);
  result *= gsl_sf_hyperg_1F1(shape2, lgam + 1, x * (rate1 - rate2));
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

double ths_p00(double s, double t,
	       double lambda0, double lambda1, double lambda2,
	       double p) {
  int n = 1;
  double result = 0;
  double cart = 0;
  double cartlast = 0;

  while (TRUE) {
    cart = R::pgamma(s, n, 1/lambda0, 1, 0) - R::pgamma(s, n + 1, 1/lambda0, 1, 0);
    cart *= sumT_p00(s, t, lambda1, lambda2, p, n);
    if (cart == R_PosInf || R_IsNaN(cart)) {
      warning("Inf or NaN happened, not converge!");
      break;
    }
    result += cart;
    if (cart == 0 && cartlast > cart) break;
    cartlast = cart;
    n++;
  }

  return result;
}

// [[Rcpp::export]]
NumericVector ths_vp00(NumericVector vs, double t,
		       double lambda0, double lambda1, double lambda2,
		       double p) {
  int n = vs.size();
  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    result[i] = ths_p00(vs[i], t, lambda0, lambda1, lambda2, p);
  }
  return result;
}

/****************************** calculate p01 *****************************/

double sumT_p01(double s, double t,
		double lambda1, double lambda2,
		double p, int n) {
  double result = 0;
  double pdsx = t - s;
  double cart = pow(1 - p, n);
  for (int k = 0; k < n + 1; ++k) {
    result += pcoga2dim_diff_shape(pdsx, k, n - k, lambda1, lambda2) * cart;
    cart *= (n - k) * p / ((k + 1) * (1 - p));
  }
  return result;
}

double ths_p01(double s, double t,
	       double lambda0, double lambda1, double lambda2,
	       double p) {
  int n = 0;
  double result = 0;
  double cart = 0;
  double cartlast = 0;

  while (TRUE) {
    cart = p * R::dgamma(s, n + 1, 1/lambda0, 0);
    cart *= sumT_p01(s, t, lambda1, lambda2, p, n);
    if (cart == R_PosInf || R_IsNaN(cart)) {
      warning("Inf or NaN happened, not converge!");
      break;
    }
    result += cart;
    if (cart == 0 && cartlast > cart) break;
    cartlast = cart;
    n++;
  }

  return result;
}

// [[Rcpp::export]]
NumericVector ths_vp01(NumericVector vs, double t,
		       double lambda0, double lambda1, double lambda2,
		       double p) {
  int n = vs.size();
  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    result[i] = ths_p01(vs[i], t, lambda0, lambda1, lambda2, p);
  }
  return result;
}
