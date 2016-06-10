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
arma::mat getLevelCPP(arma::mat A, arma::vec cuts) {
  uword ncols = A.n_cols;
  uword nrows = A.n_rows;
  uword ncuts = cuts.n_elem;

  arma::mat lev(nrows, ncols);
  for (uword i = 0; i < nrows; i++) {
    for (uword j = 0; j < ncols; j++) {
      int lev_ij = 1;
      double a_ij = A(i, j);
      for (uword k = 0; k < ncuts; k++) {
        lev_ij = a_ij > cuts(k) ? lev_ij + 1 : lev_ij;
      }
      lev(i, j) = lev_ij;
    }
  }

  return lev;
}

// [[Rcpp::export]]
arma::mat dPSCPP(arma::mat A, double alpha, arma::vec mid_points,
                 arma::vec bin_width, int threads = 1) {

  uword ns = A.n_rows; uword nt = A.n_cols;
  uword nbins = mid_points.n_elem;
  uword s; uword t; uword i;
  double integral; double llst; double psi; double logc; double logint;
  double ast;
  arma::mat ll(ns, nt);

  for (s = 0; s < ns; s++) {
    for (t = 0; t < nt; t++) {
      ast = A(s, t);
      llst = log(alpha) - log(1 - alpha) - log(ast) / (1 - alpha);
      integral = 0;
#pragma omp parallel for private(psi, logc, logint) reduction(+:integral) schedule(dynamic)
      for (i = 0; i < nbins; i++) {
        psi = PI * mid_points[i];
        logc = (log(sin(alpha * psi)) - log(sin(psi))) / (1 - alpha) +
          log(sin((1 - alpha) * psi)) - log(sin(alpha * psi));
        logint = logc - exp(logc) * pow(ast, (- alpha / (1 - alpha)));
        integral += exp(logint) * bin_width[i];
      }
      ll(s, t) = llst + log(integral);
    }
  }

  return ll;
}

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