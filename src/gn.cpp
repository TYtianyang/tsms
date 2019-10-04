#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include <math.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec grad(arma::vec X,int U,arma::vec beta
             ,int p
             ,int q
             , IntegerVector S1,int r1
             , IntegerVector S2,int r2
             , arma::mat W){
  
  int N = X.size();
  int d = W.n_rows;
  arma::vec epsi = arma::zeros(N);
  arma::mat Z = arma::zeros(p+q+r1+r2+d,N);
  arma::vec gradient = arma::zeros(p+q+r1+r2+d);
  
  for (int t=U;t<N;t++){
    arma::vec X_G(p+q+r1+r2+d);
    for (int i=0;i<p;i++){
      X_G[i] = X[t-i-1];
    }
    for (int j=0;j<q;j++){
      X_G[p+j] = epsi[t-j-1];
    }
    for (int k=0;k<r1;k++){
      X_G[p+q+k] = X[t-S1[k]];
    }
    for (int l=0;l<r2;l++){
      X_G[p+q+r1+l] = epsi[t-S2[l]];
    }
    for (int m=0;m<d;m++){
      X_G[p+q+r1+r2+m] = W(m,t);
    }
    double ac=0;
    for (int i=0;i<(p+q+r1+r2);i++){
      ac += beta[i]*X_G[i];
    }
    epsi[t] = X[t] - ac;
    
    arma::vec acc = arma::zeros(p+q+r1+r2+d);
    for (int j=0;j<q;j++){
      acc += beta[p+j]*Z.col(t-j-1);
    }
    for (int k=0;k<r2;k++){
      acc += beta[p+q+r1+k]*Z.col(t-S2[k]);
    }
    Z.col(t) = - X_G - acc;
    
    gradient += epsi[t]*Z.col(t);
  }

  return gradient;
}





