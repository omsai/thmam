#include <Rcpp.h>
#include <coga.h>
#include <R_ext/Applic.h>

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


/****************************** calculate p02 *****************************/

double new_p02(double s, double t,
	       double lambda0, double lambda1, double lambda2,
	       double p) {
  p = 1 - p;
  double cart = lambda2;
  lambda2 = lambda1;
  lambda1 = cart;
  return new_p01(s, t, lambda0, lambda1, lambda2, p);
}

// [[Rcpp::export]]
NumericVector new_vp02(NumericVector vs, double t,
		       double lambda0, double lambda1, double lambda2,
		       double p) {
  int n = vs.size();
  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    result[i] = new_p02(vs[i], t, lambda0, lambda1, lambda2, p);
  }
  return result;
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

// [[Rcpp::export]]
NumericVector new_vp10(NumericVector vs, double t,
		       double lambda0, double lambda1, double lambda2,
		       double p) {
  int n = vs.size();
  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    result[i] = new_p10(vs[i], t, lambda0, lambda1, lambda2, p);
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

// [[Rcpp::export]]
NumericVector new_vp11(NumericVector vs, double t,
		       double lambda0, double lambda1, double lambda2,
		       double p) {
  int n = vs.size();
  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    result[i] = new_p11(vs[i], t, lambda0, lambda1, lambda2, p);
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

// [[Rcpp::export]]
NumericVector new_vp12(NumericVector vs, double t,
		       double lambda0, double lambda1, double lambda2,
		       double p) {
  int n = vs.size();
  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    result[i] = new_p12(vs[i], t, lambda0, lambda1, lambda2, p);
  }
  return result;
}

/****************************** calculate p20 *****************************/

double new_p20(double s, double t,
	       double lambda0, double lambda1, double lambda2,
	       double p) {
  double cart = lambda1;
  lambda1 = lambda2;
  lambda2 = cart;
  p = 1 - p;

  return new_p10(s, t, lambda0, lambda1, lambda2, p);
}

// [[Rcpp::export]]
NumericVector new_vp20(NumericVector vs, double t,
		       double lambda0, double lambda1, double lambda2,
		       double p) {
  int n = vs.size();
  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    result[i] = new_p20(vs[i], t, lambda0, lambda1, lambda2, p);
  }
  return result;
}

/****************************** calculate p21 *****************************/

double new_p21(double s, double t,
	       double lambda0, double lambda1, double lambda2,
	       double p) {
  double cart = lambda1;
  lambda1 = lambda2;
  lambda2 = cart;
  p = 1 - p;

  return new_p11(s, t, lambda0, lambda1, lambda2, p);
}

// [[Rcpp::export]]
NumericVector new_vp21(NumericVector vs, double t,
		       double lambda0, double lambda1, double lambda2,
		       double p) {
  int n = vs.size();
  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    result[i] = new_p21(vs[i], t, lambda0, lambda1, lambda2, p);
  }
  return result;
}

/****************************** calculate p22 *****************************/

double new_p22(double s, double t,
	       double lambda0, double lambda1, double lambda2,
	       double p) {
  double cart = lambda1;
  lambda1 = lambda2;
  lambda2 = cart;
  p = 1 - p;

  return new_p12(s, t, lambda0, lambda1, lambda2, p);
}

// [[Rcpp::export]]
NumericVector new_vp22(NumericVector vs, double t,
		       double lambda0, double lambda1, double lambda2,
		       double p) {
  int n = vs.size();
  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    result[i] = new_p22(vs[i], t, lambda0, lambda1, lambda2, p);
  }
  return result;
}

/************************************************************************************
 ****************************** prepare for integrate  ******************************
 ************************************************************************************/

void new_f00(double *s, int n, void *ex) {
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
    double temp = new_p00(s[i], t, lambda0, lambda1, lambda2, p);
    double sd = sigma * sqrt(s[i]);
    for (int j = 0; j < dim; j++) temp *= R::dnorm(x[j], 0.0, sd, 0);
    s[i] = temp;
  }
}

void new_f01(double *s, int n, void *ex) {
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
    double temp = new_p01(s[i], t, lambda0, lambda1, lambda2, p);
    double sd = sigma * sqrt(s[i]);
    for (int j = 0; j < dim; j++) temp *= R::dnorm(x[j], 0.0, sd, 0);
    s[i] = temp;
  }
}

void new_f02(double *s, int n, void *ex) {
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
    double temp = new_p02(s[i], t, lambda0, lambda1, lambda2, p);
    double sd = sigma * sqrt(s[i]);
    for (int j = 0; j < dim; j++) temp *= R::dnorm(x[j], 0.0, sd, 0);
    s[i] = temp;
  }
}

void new_f10(double *s, int n, void *ex) {
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
    double temp = new_p10(s[i], t, lambda0, lambda1, lambda2, p);
    double sd = sigma * sqrt(s[i]);
    for (int j = 0; j < dim; j++) temp *= R::dnorm(x[j], 0.0, sd, 0);
    s[i] = temp;
  }
}

void new_f11(double *s, int n, void *ex) {
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
    double temp = new_p11(s[i], t, lambda0, lambda1, lambda2, p);
    double sd = sigma * sqrt(s[i]);
    for (int j = 0; j < dim; j++) temp *= R::dnorm(x[j], 0.0, sd, 0);
    s[i] = temp;
  }
}

void new_f12(double *s, int n, void *ex) {
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
    double temp = new_p12(s[i], t, lambda0, lambda1, lambda2, p);
    double sd = sigma * sqrt(s[i]);
    for (int j = 0; j < dim; j++) temp *= R::dnorm(x[j], 0.0, sd, 0);
    s[i] = temp;
  }
}

void new_f20(double *s, int n, void *ex) {
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
    double temp = new_p20(s[i], t, lambda0, lambda1, lambda2, p);
    double sd = sigma * sqrt(s[i]);
    for (int j = 0; j < dim; j++) temp *= R::dnorm(x[j], 0.0, sd, 0);
    s[i] = temp;
  }
}

void new_f21(double *s, int n, void *ex) {
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
    double temp = new_p21(s[i], t, lambda0, lambda1, lambda2, p);
    double sd = sigma * sqrt(s[i]);
    for (int j = 0; j < dim; j++) temp *= R::dnorm(x[j], 0.0, sd, 0);
    s[i] = temp;
  }
}

void new_f22(double *s, int n, void *ex) {
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
    double temp = new_p22(s[i], t, lambda0, lambda1, lambda2, p);
    double sd = sigma * sqrt(s[i]);
    for (int j = 0; j < dim; j++) temp *= R::dnorm(x[j], 0.0, sd, 0);
    s[i] = temp;
  }
}


/*****************************************************************************
 ****************************** integrate for h ******************************
 *****************************************************************************/

// [[Rcpp::export]]
NumericVector new_h00(NumericMatrix x, NumericVector t, NumericVector theta,
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
    Rdqags(new_f00, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	   &limit, &lenw, &last, iwork, work);
    value[i] = result + prod;
  }
  // free memory
  Free(ex); Free(iwork); Free(work);
  return(value);
}


// [[Rcpp::export]]
NumericVector new_h01(NumericMatrix x, NumericVector t, NumericVector theta,
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
    Rdqags(new_f01, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	   &limit, &lenw, &last, iwork, work);
    value[i] = result;
  }
  // free memory
  Free(ex); Free(iwork); Free(work);
  return(value);
}

// [[Rcpp::export]]
NumericVector new_h02(NumericMatrix x, NumericVector t, NumericVector theta,
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
    Rdqags(new_f02, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	   &limit, &lenw, &last, iwork, work);
    value[i] = result;
  }
  // free memory
  Free(ex); Free(iwork); Free(work);
  return(value);
}


// [[Rcpp::export]]
NumericVector new_h10(NumericMatrix x, NumericVector t, NumericVector theta,
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
    Rdqags(new_f10, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	   &limit, &lenw, &last, iwork, work);
    value[i] = result;
  }
  // free memory
  Free(ex); Free(iwork); Free(work);
  return(value);
}

// [[Rcpp::export]]
NumericVector new_h11(NumericMatrix x, NumericVector t, NumericVector theta,
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
    Rdqags(new_f11, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	   &limit, &lenw, &last, iwork, work);
    value[i] = result;
  }
  // free memory
  Free(ex); Free(iwork); Free(work);
  return(value);
}

// [[Rcpp::export]]
NumericVector new_h12(NumericMatrix x, NumericVector t, NumericVector theta,
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
    Rdqags(new_f12, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	   &limit, &lenw, &last, iwork, work);
    value[i] = result;
  }
  // free memory
  Free(ex); Free(iwork); Free(work);
  return(value);
}

// [[Rcpp::export]]
NumericVector new_h20(NumericMatrix x, NumericVector t, NumericVector theta,
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
    Rdqags(new_f20, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	   &limit, &lenw, &last, iwork, work);
    value[i] = result;
  }
  // free memory
  Free(ex); Free(iwork); Free(work);
  return(value);
}

// [[Rcpp::export]]
NumericVector new_h21(NumericMatrix x, NumericVector t, NumericVector theta,
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
    Rdqags(new_f21, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	   &limit, &lenw, &last, iwork, work);
    value[i] = result;
  }
  // free memory
  Free(ex); Free(iwork); Free(work);
  return(value);
}


// [[Rcpp::export]]
NumericVector new_h22(NumericMatrix x, NumericVector t, NumericVector theta,
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
    Rdqags(new_f22, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
	   &limit, &lenw, &last, iwork, work);
    value[i] = result;
  }
  // free memory
  Free(ex); Free(iwork); Free(work);
  return(value);
}


/************************************************************************************
 ************* negative log likelihood function via forward algorithm ***************
 ************************************************************************************/

// [[Rcpp::export]]
double nllk_fwd_ths(NumericVector &theta, NumericMatrix &data,
	       NumericVector &integrControl) {
  // theta lambda0, lambda1, lambda2, sigma, p
  // data diff of t and x
  int n = data.nrow(); int dim = data.ncol() - 1;
  double lambda0 = theta[0], lambda1 = theta[1], lambda2 = theta[2];
  double ps0 = 1. / lambda0 / (1. / lambda0 + 1. / lambda1 + 1. / lambda2);
  double ps1 = 1. / lambda1 / (1. / lambda0 + 1. / lambda1 + 1. / lambda2);
  double ps2 = 1. / lambda2 / (1. / lambda0 + 1. / lambda1 + 1. / lambda2);
  NumericVector tt = data.column(0);
  NumericMatrix x = data(Range(0, n - 1), Range(1, dim));
  NumericVector
    hresult00 = new_h00(x, tt, theta, integrControl),
    hresult01 = new_h01(x, tt, theta, integrControl),
    hresult02 = new_h02(x, tt, theta, integrControl),
    hresult10 = new_h10(x, tt, theta, integrControl),
    hresult11 = new_h11(x, tt, theta, integrControl),
    hresult12 = new_h12(x, tt, theta, integrControl),
    hresult20 = new_h20(x, tt, theta, integrControl),
    hresult21 = new_h21(x, tt, theta, integrControl),
    hresult22 = new_h22(x, tt, theta, integrControl);
  double alpha0 = ps0, alpha1 = ps1, alpha2 = ps2;

  // do forward algorithm
  double llk = 0.;
  for (int i = 0; i < n; i++) {
    NumericVector crow = x.row(i);
    if (is_true(all(crow == 0.))) {
      hresult00[i] = 0.;
      hresult01[i] = 0.;
      hresult02[i] = 0.;
      hresult10[i] = 0.;
      hresult11[i] = exp(-lambda1 * tt[i]);
      hresult12[i] = 0.;
      hresult20[i] = 0.;
      hresult21[i] = 0.;
      hresult22[i] = exp(-lambda2 * tt[i]);
    }
    double sumf0 = alpha0 * hresult00[i] + alpha1 * hresult10[i] + alpha2 * hresult20[i];
    double sumf1 = alpha0 * hresult01[i] + alpha1 * hresult11[i] + alpha2 * hresult21[i];
    double sumf2 = alpha0 * hresult02[i] + alpha1 * hresult12[i] + alpha2 * hresult22[i];

    double dx = sumf0 + sumf1 + sumf2;
    alpha0 = sumf0 / dx;
    alpha1 = sumf1 / dx;
    alpha2 = sumf2 / dx;
    llk += log(dx);
  }
  return(-llk);
}
