// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::export]]
double ssebCPP(arma::vec B1, arma::mat B2, arma::vec Y, int exclude,
               double alpha, double lambda) {

  double sse = 0.0;  int ns = B2.n_rows; int nknots = B2.n_cols;
  double alpha_inv = 1 / alpha;
  double echat;

  B2 = pow(B2, alpha_inv);

  uword i; uword k;
  for (i = 0; i < ns; i++) {
    echat = 0;
    for (k = 0; k < nknots; k++) {
      // B2 is already raised to alpha_inv
      echat += pow(pow(B1(k), alpha_inv) + B2(i, k), alpha);
    }

    if (i != (exclude - 1)) {
      sse += pow(Y(i) - echat, 2);
    }
  }

  if (accu(B1) != 1) {
    sse += lambda * pow((accu(B1) - 1), 2);
  }

  return sse;
}

// [[Rcpp::export]]
Rcpp::NumericVector rowSumsC(arma::mat X) {
  int nrow = X.n_rows;

  Rcpp::NumericVector out(nrow);
  for (uword i = 0; i < nrow; i++) {
    out[i] = accu(X.row(i));
  }

  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector colSumsC(arma::mat X) {
  int ncol = X.n_cols;

  Rcpp::NumericVector out(ncol);
  for (uword i = 0; i < ncol; i++) {
    out[i] = accu(X.col(i));
  }

  return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix sweepC1plus(Rcpp::NumericMatrix X,
                                Rcpp::NumericVector y) {
  int nrow = X.nrow(); int ncol = X.ncol();

  Rcpp::NumericMatrix out(nrow, ncol);
  for (uword i = 0; i < ncol; i++) {
    out(_, i) = X(_, i) + y;
  }

  return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix sweepC2plus(Rcpp::NumericMatrix X,
                                Rcpp::NumericVector y) {
  int nrow = X.nrow(); int ncol = X.ncol();

  Rcpp::NumericMatrix out(nrow, ncol);
  for (uword i = 0; i < nrow; i++) {
    out(i, _) = X(i, _) + y;
  }

  return out;
}

// [[Rcpp::export]]
arma::mat sweepC1times(arma::mat X,
                       arma::colvec y) {
  int nrow = X.n_rows; int ncol = X.n_cols;

  arma::mat out(nrow, ncol);
  for (uword i = 0; i < ncol; i++) {
    out.col(i) = X.col(i) % y;
  }

  return out;
}

// [[Rcpp::export]]
arma::mat powC(arma::mat X, double alpha) {
  return exp(alpha * log(X));
}

// [[Rcpp::export]]
arma::mat sweepC2times(arma::mat X,
                       arma::rowvec y) {
  int nrow = X.n_rows; int ncol = X.n_cols;

  arma::mat out(nrow, ncol);
  for (uword i = 0; i < nrow; i++) {
    out.row(i) = X.row(i) % y;
  }

  return out;
}

// [[Rcpp::export]]
arma::mat getECRhoAlphaC(arma::mat w, double alpha) {
  int n = w.n_rows;

  arma::mat ec(n, n); double ec_ij;
  for (uword i = 0; i < (n - 1); i++) {
    for (uword j = (i + 1); j < n; j++) {
      ec_ij = accu(pow(w.row(i) + w.row(j), alpha));
      ec(i, j) = ec_ij;
      ec(j, i) = ec_ij;
    }
  }

  return(ec);
}