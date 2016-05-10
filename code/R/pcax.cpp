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
arma::vec rowSumsC(arma::mat X) {
  int nrow = X.n_rows; int ncol = X.n_cols;

  arma::vec out(nrow);
  for (uword i = 0; i < nrow; i++) {
    double total = 0;
    for (uword j = 0; j < ncol; j++) {
      total += X(i, j);
    }
    out[i] = total;
  }

  return out;
}

// [[Rcpp::export]]
arma::mat sweepC2plus(arma::mat X, arma::rowvec y) {
  int nrow = X.n_rows; int ncol = X.n_cols;

  arma::mat out(nrow, ncol);
  for (uword i = 0; i < nrow; i++) {
    out.row(i) = X.row(i) + y;
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