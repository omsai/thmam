#include <Rcpp.h>
#include <coga.h>

using namespace Rcpp;

/**************************************************************************
 ****************************** calculate p  ******************************
 **************************************************************************/

/****************************** calculate p00 *****************************/

double sumT_p00(double s, double t,
	    double lambda1, double lambda2,
	    double p, int n) {
  double result = 0;
  NumericVector dco_x = NumericVector::create(t - s);
  NumericVector rate  = NumericVector::create(lambda1, lambda2);
  NumericVector shape(2);
  for (int k = 0; k < n + 1; ++k) {
    shape[0] = k;
    shape[1] = n - k;
    result += coga::dcoga(dco_x, shape, rate)[0] * R::choose(n, k) * pow(p, k) * pow(1 - p, n - k);
  }
  return result;
}

// [[Rcpp::export]]
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
    if (cart == 0) break;
    n++;
  }

  return result;
}


/****************************** calculate p01 *****************************/

double sumT_p01(double s, double t,
		double lambda1, double lambda2,
		double p, int n) {
  double result = 0;
  NumericVector pco_x = NumericVector::create(t - s);
  NumericVector rate = NumericVector::create(lambda1, lambda2);
  NumericVector shape1(2);
  NumericVector shape2(2);
  double cart = 0;
  for (int k = 0; k < n + 1; ++k) {
    shape1[0] = k;
    shape1[1] = n - k;
    shape2[0] = k + 1;
    shape2[1] = n - k;
    cart = coga::pcoga(pco_x, shape1, rate)[0] - coga::pcoga(pco_x, shape2, rate)[0];
    cart *= R::choose(n, k) * pow(p, k) * pow(1 - p, n - k);
    result += cart;
  }
  return result;
}

// [[Rcpp::export]]
double new_p01(double s, double t,
	       double lambda0, double lambda1, double lambda2,
	       double p) {
  int n = 1;
  double result = 0;
  double cart = 0;

  while (TRUE) {
    cart = p * R::dgamma(s, n + 1, 1/lambda0, 0) * sumT_p01(s, t, lambda1, lambda2, p, n);
    result += cart;
    if (cart == 0) break;
    n++;
  }

  NumericVector shape = NumericVector::create(1, 0);
  NumericVector pco_x = NumericVector::create(t - s);
  NumericVector rate = NumericVector::create(lambda1, lambda2);
  return result + (p * R::dgamma(s, 1, 1/lambda0, 0) * (1 - coga::pcoga(pco_x, shape, rate)[0]));
}

/****************************** calculate p02 *****************************/

// [[Rcpp::export]]
double new_p02(double s, double t,
	       double lambda0, double lambda1, double lambda2,
	       double p) {
  p = 1 - p;
  double cart = lambda2;
  lambda2 = lambda1;
  lambda1 = cart;
  return new_p01(s, t, lambda0, lambda1, lambda2, p);
}

/****************************** calculate p10 *****************************/

double sumT_p10(double s, double t,
	    double lambda1, double lambda2,
	    double p, int n) {
  double result = 0;
  NumericVector dco_x = NumericVector::create(t - s);
  NumericVector rate  = NumericVector::create(lambda1, lambda2);
  NumericVector shape(2);
  for (int k = 0; k < n + 1; ++k) {
    shape[0] = k + 1;
    shape[1] = n - k;
    result += coga::dcoga(dco_x, shape, rate)[0] * R::choose(n, k) * pow(p, k) * pow(1 - p, n - k);
  }
  return result;
}

// [[Rcpp::export]]
double new_p10(double s, double t,
	       double lambda0, double lambda1, double lambda2,
	       double p) {
  int n = 0;
  double result = 0;
  double cart = 0;

  while (TRUE) {
    cart = R::pgamma(s, n, 1/lambda0, 1, 0) - R::pgamma(s, n + 1, 1/lambda0, 1, 0);
    cart *= sumT_p10(s, t, lambda1, lambda2, p, n);
    result += cart;
    if (cart == 0) break;
    n++;
  }

  return result;
}


/****************************** calculate p11 *****************************/

double sumT_p11(double s, double t,
		double lambda1, double lambda2,
		double p, int n) {
  double result = 0;
  NumericVector pco_x = NumericVector::create(t - s);
  NumericVector rate = NumericVector::create(lambda1, lambda2);
  NumericVector shape1(2);
  NumericVector shape2(2);
  double cart = 0;
  for (int k = 0; k < n; ++k) {
    shape1[0] = k + 1;
    shape1[1] = n - 1 - k;
    shape2[0] = k + 2;
    shape2[1] = n - 1 - k;
    cart = coga::pcoga(pco_x, shape1, rate)[0] - coga::pcoga(pco_x, shape2, rate)[0];
    cart *= R::choose(n - 1, k) * pow(p, k) * pow(1 - p, n - 1 - k);
    result += cart;
  }
  return result;
}

// [[Rcpp::export]]
double new_p11(double s, double t,
	       double lambda0, double lambda1, double lambda2,
	       double p) {
  int n = 1;
  double result = 0;
  double cart = 0;

  while (TRUE) {
    cart = p * R::dgamma(s, n, 1/lambda0, 0) * sumT_p11(s, t, lambda1, lambda2, p, n);
    result += cart;
    if (cart == 0) break;
    n++;
  }

  return result;
}

/****************************** calculate p12 *****************************/

double sumT_p12(double s, double t,
		double lambda1, double lambda2,
		double p, int n) {
  double result = 0;
  NumericVector pco_x = NumericVector::create(t - s);
  NumericVector rate = NumericVector::create(lambda1, lambda2);
  NumericVector shape1(2);
  NumericVector shape2(2);
  double cart = 0;
  for (int k = 0; k < n; ++k) {
    shape1[0] = k + 1;
    shape1[1] = n - 1 - k;
    shape2[0] = k + 1;
    shape2[1] = n - k;
    cart = coga::pcoga(pco_x, shape1, rate)[0] - coga::pcoga(pco_x, shape2, rate)[0];
    cart *= R::choose(n - 1, k) * pow(p, k) * pow(1 - p, n - 1 - k);
    result += cart;
  }
  return result;
}

// [[Rcpp::export]]
double new_p12(double s, double t,
	       double lambda0, double lambda1, double lambda2,
	       double p) {
  int n = 1;
  double result = 0;
  double cart = 0;

  while (TRUE) {
    cart = (1 - p) * R::dgamma(s, n, 1/lambda0, 0) * sumT_p12(s, t, lambda1, lambda2, p, n);
    result += cart;
    if (cart == 0) break;
    n++;
  }

  return result;
}

/****************************** calculate p20 *****************************/

// [[Rcpp::export]]
double new_p20(double s, double t,
	       double lambda0, double lambda1, double lambda2,
	       double p) {
  double cart = lambda1;
  lambda1 = lambda2;
  lambda2 = cart;
  p = 1 - p;

  return new_p10(s, t, lambda0, lambda1, lambda2, p);
}

/****************************** calculate p21 *****************************/

// [[Rcpp::export]]
double new_p21(double s, double t,
	       double lambda0, double lambda1, double lambda2,
	       double p) {
  double cart = lambda1;
  lambda1 = lambda2;
  lambda2 = cart;
  p = 1 - p;

  return new_p11(s, t, lambda0, lambda1, lambda2, p);
}

/****************************** calculate p22 *****************************/

// [[Rcpp::export]]
double new_p22(double s, double t,
	       double lambda0, double lambda1, double lambda2,
	       double p) {
  double cart = lambda1;
  lambda1 = lambda2;
  lambda2 = cart;
  p = 1 - p;

  return new_p12(s, t, lambda0, lambda1, lambda2, p);
}
