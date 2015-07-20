// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::export]]
arma::mat ifelsematCPP(arma::mat x, double tol) {
  uword ncols = x.n_cols;
  uword nrows = x.n_rows;
  
  for (uword i = 0; i < nrows; i++) {
    for (uword j = 0; j < ncols; j++) {
      x(i, j) = fabs(x(i, j)) < tol ? 0 : x(i, j);
    }
  }
  
  return x;
}

// [[Rcpp::export]]
arma::vec ifelsevecCPP(arma::vec x, double tol) {
  uword n = x.n_elem;
  
  for (uword i = 0; i < n; i++) {
    x[i] = fabs(x[i]) < tol ? 0 : x[i];
  }
  
  return x;
}

// [[Rcpp::export]]
arma::mat getLevelCPP(arma::mat a, arma::vec cuts) {
  uword ncols = a.n_cols;
  uword nrows = a.n_rows;
  uword ncuts = cuts.n_elem;
  
  arma::mat lev(nrows, ncols);
  for (uword i = 0; i < nrows; i++) {
    for (uword j = 0; j < ncols; j++) {
      int lev_ij = 1;
      double a_ij = a(i, j);
      for (uword k = 0; k < ncuts; k++) {
        lev_ij = a_ij > cuts(k) ? lev_ij + 1 : lev_ij;
      }
      lev(i, j) = lev_ij;
    }
  }
  
  return lev;
}