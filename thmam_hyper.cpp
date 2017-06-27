// [[Rcpp::depends(RcppGSL)]]

#include <Rcpp.h>
#include <RcppGSL.h>
#include <gsl/gsl_sf_hyperg.h>
#include <R_ext/Applic.h>

using namespace Rcpp;

/***************** basic formula for coga *****************/

// [[Rcpp::export]]
double dcoga2dim(double x, double shape1, double shape2,
		       double rate1, double rate2) {
  // transfer rate to scale
  double beta1 = 1 / rate1;
  double beta2 = 1 / rate2;
  
  /*
  // handle one shape is 0
  if (shape1 == 0) return R::dgamma(x, shape2, beta2, 0);
  if (shape2 == 0) return R::dgamma(x, shape1, beta1, 0);
  */
  
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
double pcoga2dim_diff_shape (double x,
			     double shape1, double shape2,
			     double rate1, double rate2) {
  /*
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
  */
  
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
  double cart = pow(1 - p, n);
  for (int k = 0; k < n + 1; ++k) {
    result += dcoga2dim(t - s, k, n - k, lambda1, lambda2) * cart;
    cart *= p * (n - k) / ((k + 1) * (1 - p));
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

/****************************** calculate p02 *****************************/

double ths_p02(double s, double t,
	       double lambda0, double lambda1, double lambda2,
	       double p) {
  return ths_p01(s, t, lambda0, lambda2, lambda1, 1 - p);
}

// [[Rcpp::export]]
NumericVector ths_vp02(NumericVector vs, double t,
		       double lambda0, double lambda1, double lambda2,
		       double p) {
  return ths_vp01(vs, t, lambda0, lambda2, lambda1, 1 - p);
}


/****************************** calculate p10 *****************************/
double sumT_p10(double s, double t,
	    double lambda1, double lambda2,
	    double p, int n) {
  double result = 0;
  double cart = pow(1 - p, n);
  for (int k = 0; k < n + 1; ++k) {
    result += dcoga2dim(t - s, k + 1, n - k, lambda1, lambda2) * cart;
    cart *= p * (n - k) / ((k + 1) * (1 - p));
  }
  return result;
}

double ths_p10(double s, double t,
	       double lambda0, double lambda1, double lambda2,
	       double p) {
  int n = 0;
  double result = 0;
  double cart = 0;
  double cartlast = 0;

  while (TRUE) {
    cart = R::pgamma(s, n, 1/lambda0, 1, 0) - R::pgamma(s, n + 1, 1/lambda0, 1, 0);
    cart *= sumT_p10(s, t, lambda1, lambda2, p, n);
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
NumericVector ths_vp10(NumericVector vs, double t,
		       double lambda0, double lambda1, double lambda2,
		       double p) {
  int n = vs.size();
  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    result[i] = ths_p10(vs[i], t, lambda0, lambda1, lambda2, p);
  }
  return result;
}


/****************************** calculate p11 *****************************/

double sumT_p11(double s, double t,
		double lambda1, double lambda2,
		double p, int n) {
  double result = 0;
  double pdsx = t - s;
  double cart = pow(1 - p, n - 1);
  for (int k = 0; k < n; ++k) {
    result += pcoga2dim_diff_shape(pdsx, k + 1, n - k - 1, lambda1, lambda2) * cart;
    cart *= (n - k - 1) * p / ((k + 1) * (1 - p));
  }
  return result;
}

double ths_p11(double s, double t,
	       double lambda0, double lambda1, double lambda2,
	       double p) {
  int n = 1;
  double result = 0;
  double cart = 0;
  double cartlast = 0;

  while (TRUE) {
    cart = p * R::dgamma(s, n, 1/lambda0, 0);
    cart *= sumT_p11(s, t, lambda1, lambda2, p, n);
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
NumericVector ths_vp11(NumericVector vs, double t,
		       double lambda0, double lambda1, double lambda2,
		       double p) {
  int n = vs.size();
  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    result[i] = ths_p11(vs[i], t, lambda0, lambda1, lambda2, p);
  }
  return result;
}


/****************************** calculate p12 *****************************/

double sumT_p12(double s, double t,
		double lambda1, double lambda2,
		double p, int n) {
  double result = 0;
  double pdsx = t - s;
  double cart = pow(1 - p, n - 1);
  for (int k = 0; k < n; ++k) {
    result += pcoga2dim_diff_shape(pdsx, n - k - 1, k + 1, lambda2, lambda1) * cart;
    cart *= (n - k - 1) * p / ((k + 1) * (1 - p));
  }
  return result;
}

double ths_p12(double s, double t,
	       double lambda0, double lambda1, double lambda2,
	       double p) {
  int n = 1;
  double result = 0;
  double cart = 0;
  double cartlast = 0;

  while (TRUE) {
    cart = (1 - p) * R::dgamma(s, n, 1/lambda0, 0);
    cart *= sumT_p12(s, t, lambda1, lambda2, p, n);
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
NumericVector ths_vp12(NumericVector vs, double t,
		       double lambda0, double lambda1, double lambda2,
		       double p) {
  int n = vs.size();
  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    result[i] = ths_p12(vs[i], t, lambda0, lambda1, lambda2, p);
  }
  return result;
}

/****************************** calculate p20 *****************************/

double ths_p20(double s, double t,
	       double lambda0, double lambda1, double lambda2,
	       double p) {
  return ths_p10(s, t, lambda0, lambda2, lambda1, 1 - p);
}

// [[Rcpp::export]]
NumericVector ths_vp20(NumericVector vs, double t,
			double lambda0, double lambda1, double lambda2,
			double p) {
  return ths_vp10(vs, t, lambda0, lambda2, lambda1, 1 - p);
}


/****************************** calculate p21 *****************************/

double ths_p21(double s, double t,
	       double lambda0, double lambda1, double lambda2,
	       double p) {
  return ths_p11(s, t, lambda0, lambda2, lambda1, 1 - p);
}

// [[Rcpp::export]]
NumericVector ths_vp21(NumericVector vs, double t,
			double lambda0, double lambda1, double lambda2,
			double p) {
  return ths_vp11(vs, t, lambda0, lambda2, lambda1, 1 - p);
}

/****************************** calculate p22 *****************************/

double ths_p22(double s, double t,
	       double lambda0, double lambda1, double lambda2,
	       double p) {
  return ths_p12(s, t, lambda0, lambda2, lambda1, 1 - p);
}

// [[Rcpp::export]]
NumericVector ths_vp22(NumericVector vs, double t,
			double lambda0, double lambda1, double lambda2,
			double p) {
  return ths_vp12(vs, t, lambda0, lambda2, lambda1, 1 - p);
}


/************************************************************************************
 ****************************** prepare for integrate  ******************************
 ************************************************************************************/

void ths_f00(double *s, int n, void *ex) {
  double *ptr = (double *) ex;
  double t       = ptr[0];
  double sigma   = ptr[1];
  double lambda0 = ptr[2];
  double lambda1 = ptr[3];
  double lambda2 = ptr[4];
  double p       = ptr[5];
  int    dim     = (int) ptr[6];
  double *x      = ptr + 7;
  for (int i = 0; i < n; i++) {
    double temp = ths_p00(s[i], t, lambda0, lambda1, lambda2, p);
    double sd = sigma * sqrt(s[i]);
    for (int j = 0; j < dim; j++) temp *= R::dnorm(x[j], 0.0, sd, 0);
    s[i] = temp;
  }
}

void ths_f01(double *s, int n, void *ex) {
  double *ptr = (double *) ex;
  double t       = ptr[0];
  double sigma   = ptr[1];
  double lambda0 = ptr[2];
  double lambda1 = ptr[3];
  double lambda2 = ptr[4];
  double p       = ptr[5];
  int    dim     = (int) ptr[6];
  double *x      = ptr + 7;
  for (int i = 0; i < n; i++) {
    double temp = ths_p01(s[i], t, lambda0, lambda1, lambda2, p);
    double sd = sigma * sqrt(s[i]);
    for (int j = 0; j < dim; j++) temp *= R::dnorm(x[j], 0.0, sd, 0);
    s[i] = temp;
  }
}

void ths_f02(double *s, int n, void *ex) {
  double *ptr = (double *) ex;
  double t       = ptr[0];
  double sigma   = ptr[1];
  double lambda0 = ptr[2];
  double lambda1 = ptr[3];
  double lambda2 = ptr[4];
  double p       = ptr[5];
  int    dim     = (int) ptr[6];
  double *x      = ptr + 7;
  for (int i = 0; i < n; i++) {
    double temp = ths_p02(s[i], t, lambda0, lambda1, lambda2, p);
    double sd = sigma * sqrt(s[i]);
    for (int j = 0; j < dim; j++) temp *= R::dnorm(x[j], 0.0, sd, 0);
    s[i] = temp;
  }
}

void ths_f10(double *s, int n, void *ex) {
  double *ptr = (double *) ex;
  double t       = ptr[0];
  double sigma   = ptr[1];
  double lambda0 = ptr[2];
  double lambda1 = ptr[3];
  double lambda2 = ptr[4];
  double p       = ptr[5];
  int    dim     = (int) ptr[6];
  double *x      = ptr + 7;
  for (int i = 0; i < n; i++) {
    double temp = ths_p10(s[i], t, lambda0, lambda1, lambda2, p);
    double sd = sigma * sqrt(s[i]);
    for (int j = 0; j < dim; j++) temp *= R::dnorm(x[j], 0.0, sd, 0);
    s[i] = temp;
  }
}

void ths_f11(double *s, int n, void *ex) {
  double *ptr = (double *) ex;
  double t       = ptr[0];
  double sigma   = ptr[1];
  double lambda0 = ptr[2];
  double lambda1 = ptr[3];
  double lambda2 = ptr[4];
  double p       = ptr[5];
  int    dim     = (int) ptr[6];
  double *x      = ptr + 7;
  for (int i = 0; i < n; i++) {
    double temp = ths_p11(s[i], t, lambda0, lambda1, lambda2, p);
    double sd = sigma * sqrt(s[i]);
    for (int j = 0; j < dim; j++) temp *= R::dnorm(x[j], 0.0, sd, 0);
    s[i] = temp;
  }
}

void ths_f12(double *s, int n, void *ex) {
  double *ptr = (double *) ex;
  double t       = ptr[0];
  double sigma   = ptr[1];
  double lambda0 = ptr[2];
  double lambda1 = ptr[3];
  double lambda2 = ptr[4];
  double p       = ptr[5];
  int    dim     = (int) ptr[6];
  double *x      = ptr + 7;
  for (int i = 0; i < n; i++) {
    double temp = ths_p12(s[i], t, lambda0, lambda1, lambda2, p);
    double sd = sigma * sqrt(s[i]);
    for (int j = 0; j < dim; j++) temp *= R::dnorm(x[j], 0.0, sd, 0);
    s[i] = temp;
  }
}

void ths_f20(double *s, int n, void *ex) {
  double *ptr = (double *) ex;
  double t       = ptr[0];
  double sigma   = ptr[1];
  double lambda0 = ptr[2];
  double lambda1 = ptr[3];
  double lambda2 = ptr[4];
  double p       = ptr[5];
  int    dim     = (int) ptr[6];
  double *x      = ptr + 7;
  for (int i = 0; i < n; i++) {
    double temp = ths_p20(s[i], t, lambda0, lambda1, lambda2, p);
    double sd = sigma * sqrt(s[i]);
    for (int j = 0; j < dim; j++) temp *= R::dnorm(x[j], 0.0, sd, 0);
    s[i] = temp;
  }
}

void ths_f21(double *s, int n, void *ex) {
  double *ptr = (double *) ex;
  double t       = ptr[0];
  double sigma   = ptr[1];
  double lambda0 = ptr[2];
  double lambda1 = ptr[3];
  double lambda2 = ptr[4];
  double p       = ptr[5];
  int    dim     = (int) ptr[6];
  double *x      = ptr + 7;
  for (int i = 0; i < n; i++) {
    double temp = ths_p21(s[i], t, lambda0, lambda1, lambda2, p);
    double sd = sigma * sqrt(s[i]);
    for (int j = 0; j < dim; j++) temp *= R::dnorm(x[j], 0.0, sd, 0);
    s[i] = temp;
  }
}

void ths_f22(double *s, int n, void *ex) {
  double *ptr = (double *) ex;
  double t       = ptr[0];
  double sigma   = ptr[1];
  double lambda0 = ptr[2];
  double lambda1 = ptr[3];
  double lambda2 = ptr[4];
  double p       = ptr[5];
  int    dim     = (int) ptr[6];
  double *x      = ptr + 7;
  for (int i = 0; i < n; i++) {
    double temp = ths_p22(s[i], t, lambda0, lambda1, lambda2, p);
    double sd = sigma * sqrt(s[i]);
    for (int j = 0; j < dim; j++) temp *= R::dnorm(x[j], 0.0, sd, 0);
    s[i] = temp;
  }
}


/*****************************************************************************
 ****************************** integrate for h ******************************
 *****************************************************************************/

// [[Rcpp::export]]
NumericVector ths_h00(NumericMatrix x, NumericVector t, NumericVector theta,
		  NumericVector integrControl) {
  int dim = x.ncol(), n = x.nrow();
  double lambda0 = theta[0], lambda1 = theta[1], lambda2 = theta[2];
  double sigma = theta[3], p = theta[4];

  /* set up for Rdqags */
  double *ex = Calloc(7 + dim, double);
  double a = 0., b; // = t;
  // input
  double epsabs = integrControl[0], epsrel = integrControl[1];
  int limit = (int) integrControl[2]; // subdivision in argument of R integrate
  // output
  double result, abserr; // integrate()$value and integrate()$abs.error
  int last, ier;  // integrate()$subdivision and integrate()$message
  int neval; // number of evaluation of integrand
  // working arrays
  int lenw = 4 * limit, *iwork = Calloc(limit, int);
  double *work = Calloc(lenw,  double);
  /* done setting for Rdqags */

  ex[1] = sigma; ex[2] = lambda0; ex[3] = lambda1; ex[4] = lambda2;
  ex[5] = p; ex[6] = (double) dim;
  NumericVector value(n);
  for (int i = 0; i < n; i++) {
    // for atom
    double sd = sigma * sqrt(t[i]);
    double prod = exp(-lambda0 * t[i]);
    for (int j = 0; j < dim; j++) {
      ex[7 + j] = x(i, j);
      prod *= R::dnorm(x(i, j), 0.0, sd, 0);
    }
    // do integrate
    b = t[i]; ex[0] = t[i];
    Rdqags(ths_f00, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	   &limit, &lenw, &last, iwork, work);
    value[i] = result + prod;
  }
  // free memory
  Free(ex); Free(iwork); Free(work);
  return(value);
}


// [[Rcpp::export]]
NumericVector ths_h01(NumericMatrix x, NumericVector t, NumericVector theta,
		  NumericVector integrControl) {
  int dim = x.ncol(), n = x.nrow();
  double lambda0 = theta[0], lambda1 = theta[1], lambda2 = theta[2];
  double sigma = theta[3], p = theta[4];

  /* set up for Rdqags */
  double *ex = Calloc(7 + dim, double);
  double a = 0., b; // = t;
  // input
  double epsabs = integrControl[0], epsrel = integrControl[1];
  int limit = (int) integrControl[2]; // subdivision in argument of R integrate
  // output
  double result, abserr; // integrate()$value and integrate()$abs.error
  int last, ier;  // integrate()$subdivision and integrate()$message
  int neval; // number of evaluation of integrand
  // working arrays
  int lenw = 4 * limit, *iwork = Calloc(limit, int);
  double *work = Calloc(lenw,  double);
  /* done setting for Rdqags */

  ex[1] = sigma; ex[2] = lambda0; ex[3] = lambda1; ex[4] = lambda2;
  ex[5] = p; ex[6] = (double) dim;
  NumericVector value(n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < dim; j++) ex[7 + j] = x(i, j);
    // do integrate
    b = t[i]; ex[0] = t[i];
    Rdqags(ths_f01, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	   &limit, &lenw, &last, iwork, work);
    value[i] = result;
  }
  // free memory
  Free(ex); Free(iwork); Free(work);
  return(value);
}

// [[Rcpp::export]]
NumericVector ths_h02(NumericMatrix x, NumericVector t, NumericVector theta,
		  NumericVector integrControl) {
  int dim = x.ncol(), n = x.nrow();
  double lambda0 = theta[0], lambda1 = theta[1], lambda2 = theta[2];
  double sigma = theta[3], p = theta[4];

  /* set up for Rdqags */
  double *ex = Calloc(7 + dim, double);
  double a = 0., b; // = t;
  // input
  double epsabs = integrControl[0], epsrel = integrControl[1];
  int limit = (int) integrControl[2]; // subdivision in argument of R integrate
  // output
  double result, abserr; // integrate()$value and integrate()$abs.error
  int last, ier;  // integrate()$subdivision and integrate()$message
  int neval; // number of evaluation of integrand
  // working arrays
  int lenw = 4 * limit, *iwork = Calloc(limit, int);
  double *work = Calloc(lenw,  double);
  /* done setting for Rdqags */

  ex[1] = sigma; ex[2] = lambda0; ex[3] = lambda1; ex[4] = lambda2;
  ex[5] = p; ex[6] = (double) dim;
  NumericVector value(n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < dim; j++) ex[7 + j] = x(i, j);
    // do integrate
    b = t[i]; ex[0] = t[i];
    Rdqags(ths_f02, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	   &limit, &lenw, &last, iwork, work);
    value[i] = result;
  }
  // free memory
  Free(ex); Free(iwork); Free(work);
  return(value);
}


// [[Rcpp::export]]
NumericVector ths_h10(NumericMatrix x, NumericVector t, NumericVector theta,
		  NumericVector integrControl) {
  int dim = x.ncol(), n = x.nrow();
  double lambda0 = theta[0], lambda1 = theta[1], lambda2 = theta[2];
  double sigma = theta[3], p = theta[4];

  /* set up for Rdqags */
  double *ex = Calloc(7 + dim, double);
  double a = 0., b; // = t;
  // input
  double epsabs = integrControl[0], epsrel = integrControl[1];
  int limit = (int) integrControl[2]; // subdivision in argument of R integrate
  // output
  double result, abserr; // integrate()$value and integrate()$abs.error
  int last, ier;  // integrate()$subdivision and integrate()$message
  int neval; // number of evaluation of integrand
  // working arrays
  int lenw = 4 * limit, *iwork = Calloc(limit, int);
  double *work = Calloc(lenw,  double);
  /* done setting for Rdqags */

  ex[1] = sigma; ex[2] = lambda0; ex[3] = lambda1; ex[4] = lambda2;
  ex[5] = p; ex[6] = (double) dim;
  NumericVector value(n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < dim; j++) ex[7 + j] = x(i, j);
    // do integrate
    b = t[i]; ex[0] = t[i];
    Rdqags(ths_f10, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	   &limit, &lenw, &last, iwork, work);
    value[i] = result;
  }
  // free memory
  Free(ex); Free(iwork); Free(work);
  return(value);
}

// [[Rcpp::export]]
NumericVector ths_h11(NumericMatrix x, NumericVector t, NumericVector theta,
		  NumericVector integrControl) {
  int dim = x.ncol(), n = x.nrow();
  double lambda0 = theta[0], lambda1 = theta[1], lambda2 = theta[2];
  double sigma = theta[3], p = theta[4];

  /* set up for Rdqags */
  double *ex = Calloc(7 + dim, double);
  double a = 0., b; // = t;
  // input
  double epsabs = integrControl[0], epsrel = integrControl[1];
  int limit = (int) integrControl[2]; // subdivision in argument of R integrate
  // output
  double result, abserr; // integrate()$value and integrate()$abs.error
  int last, ier;  // integrate()$subdivision and integrate()$message
  int neval; // number of evaluation of integrand
  // working arrays
  int lenw = 4 * limit, *iwork = Calloc(limit, int);
  double *work = Calloc(lenw,  double);
  /* done setting for Rdqags */

  ex[1] = sigma; ex[2] = lambda0; ex[3] = lambda1; ex[4] = lambda2;
  ex[5] = p; ex[6] = (double) dim;
  NumericVector value(n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < dim; j++) ex[7 + j] = x(i, j);
    // do integrate
    b = t[i]; ex[0] = t[i];
    Rdqags(ths_f11, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	   &limit, &lenw, &last, iwork, work);
    value[i] = result;
  }
  // free memory
  Free(ex); Free(iwork); Free(work);
  return(value);
}

// [[Rcpp::export]]
NumericVector ths_h12(NumericMatrix x, NumericVector t, NumericVector theta,
		  NumericVector integrControl) {
  int dim = x.ncol(), n = x.nrow();
  double lambda0 = theta[0], lambda1 = theta[1], lambda2 = theta[2];
  double sigma = theta[3], p = theta[4];

  /* set up for Rdqags */
  double *ex = Calloc(7 + dim, double);
  double a = 0., b; // = t;
  // input
  double epsabs = integrControl[0], epsrel = integrControl[1];
  int limit = (int) integrControl[2]; // subdivision in argument of R integrate
  // output
  double result, abserr; // integrate()$value and integrate()$abs.error
  int last, ier;  // integrate()$subdivision and integrate()$message
  int neval; // number of evaluation of integrand
  // working arrays
  int lenw = 4 * limit, *iwork = Calloc(limit, int);
  double *work = Calloc(lenw,  double);
  /* done setting for Rdqags */

  ex[1] = sigma; ex[2] = lambda0; ex[3] = lambda1; ex[4] = lambda2;
  ex[5] = p; ex[6] = (double) dim;
  NumericVector value(n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < dim; j++) ex[7 + j] = x(i, j);
    // do integrate
    b = t[i]; ex[0] = t[i];
    Rdqags(ths_f12, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	   &limit, &lenw, &last, iwork, work);
    value[i] = result;
  }
  // free memory
  Free(ex); Free(iwork); Free(work);
  return(value);
}

// [[Rcpp::export]]
NumericVector ths_h20(NumericMatrix x, NumericVector t, NumericVector theta,
		  NumericVector integrControl) {
  int dim = x.ncol(), n = x.nrow();
  double lambda0 = theta[0], lambda1 = theta[1], lambda2 = theta[2];
  double sigma = theta[3], p = theta[4];

  /* set up for Rdqags */
  double *ex = Calloc(7 + dim, double);
  double a = 0., b; // = t;
  // input
  double epsabs = integrControl[0], epsrel = integrControl[1];
  int limit = (int) integrControl[2]; // subdivision in argument of R integrate
  // output
  double result, abserr; // integrate()$value and integrate()$abs.error
  int last, ier;  // integrate()$subdivision and integrate()$message
  int neval; // number of evaluation of integrand
  // working arrays
  int lenw = 4 * limit, *iwork = Calloc(limit, int);
  double *work = Calloc(lenw,  double);
  /* done setting for Rdqags */

  ex[1] = sigma; ex[2] = lambda0; ex[3] = lambda1; ex[4] = lambda2;
  ex[5] = p; ex[6] = (double) dim;
  NumericVector value(n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < dim; j++) ex[7 + j] = x(i, j);
    // do integrate
    b = t[i]; ex[0] = t[i];
    Rdqags(ths_f20, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	   &limit, &lenw, &last, iwork, work);
    value[i] = result;
  }
  // free memory
  Free(ex); Free(iwork); Free(work);
  return(value);
}

// [[Rcpp::export]]
NumericVector ths_h21(NumericMatrix x, NumericVector t, NumericVector theta,
		  NumericVector integrControl) {
  int dim = x.ncol(), n = x.nrow();
  double lambda0 = theta[0], lambda1 = theta[1], lambda2 = theta[2];
  double sigma = theta[3], p = theta[4];

  /* set up for Rdqags */
  double *ex = Calloc(7 + dim, double);
  double a = 0., b; // = t;
  // input
  double epsabs = integrControl[0], epsrel = integrControl[1];
  int limit = (int) integrControl[2]; // subdivision in argument of R integrate
  // output
  double result, abserr; // integrate()$value and integrate()$abs.error
  int last, ier;  // integrate()$subdivision and integrate()$message
  int neval; // number of evaluation of integrand
  // working arrays
  int lenw = 4 * limit, *iwork = Calloc(limit, int);
  double *work = Calloc(lenw,  double);
  /* done setting for Rdqags */

  ex[1] = sigma; ex[2] = lambda0; ex[3] = lambda1; ex[4] = lambda2;
  ex[5] = p; ex[6] = (double) dim;
  NumericVector value(n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < dim; j++) ex[7 + j] = x(i, j);
    // do integrate
    b = t[i]; ex[0] = t[i];
    Rdqags(ths_f21, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	   &limit, &lenw, &last, iwork, work);
    value[i] = result;
  }
  // free memory
  Free(ex); Free(iwork); Free(work);
  return(value);
}


// [[Rcpp::export]]
NumericVector ths_h22(NumericMatrix x, NumericVector t, NumericVector theta,
		  NumericVector integrControl) {
  int dim = x.ncol(), n = x.nrow();
  double lambda0 = theta[0], lambda1 = theta[1], lambda2 = theta[2];
  double sigma = theta[3], p = theta[4];

  /* set up for Rdqags */
  double *ex = Calloc(7 + dim, double);
  double a = 0., b; // = t;
  // input
  double epsabs = integrControl[0], epsrel = integrControl[1];
  int limit = (int) integrControl[2]; // subdivision in argument of R integrate
  // output
  double result, abserr; // integrate()$value and integrate()$abs.error
  int last, ier;  // integrate()$subdivision and integrate()$message
  int neval; // number of evaluation of integrand
  // working arrays
  int lenw = 4 * limit, *iwork = Calloc(limit, int);
  double *work = Calloc(lenw,  double);
  /* done setting for Rdqags */

  ex[1] = sigma; ex[2] = lambda0; ex[3] = lambda1; ex[4] = lambda2;
  ex[5] = p; ex[6] = (double) dim;
  NumericVector value(n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < dim; j++) ex[7 + j] = x(i, j);
    // do integrate
    b = t[i]; ex[0] = t[i];
    Rdqags(ths_f22, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	   &limit, &lenw, &last, iwork, work);
    value[i] = result;
  }
  // free memory
  Free(ex); Free(iwork); Free(work);
  return(value);
}
