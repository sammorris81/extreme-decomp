// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::export]]
NumericMatrix madogramCPP(NumericMatrix data) {
  // Computing the madogram with thresholded data
  
  uword ns = data.nrow(); uword nt = data.ncol();
  uword i; uword j; uword k;  // site 1, site 2, years
  NumericMatrix mado(ns, ns);

  for (i = 0; i < (ns - 1); i++) {
    for (j = i + 1; j < ns; j++) {
      int nobs = 0;
      for (k= 0; k < nt; k++) {
        if (!Rcpp::NumericMatrix::is_na(data(i, k)) && 
            !Rcpp::NumericMatrix::is_na(data(j, k))) {  // obs at both sites
            if (data(i, k) > 0 || data(j, k) > 0) {  // at least one exceedance
              mado(i, j) += fabs(data(i, k) - data(j, k));
              nobs ++;
            }
        } 
      }
      mado(i, j) *= 0.5 / (double) nobs;
      mado(j, i) = mado(i, j);
    }
  }
    
  return mado;
}

