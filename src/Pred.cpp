#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include <math.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec Pred(arma::vec X,int U,int predict_t,
                  arma::vec phi,int p,
                  arma::vec psi,int q,
                  arma::vec tau1, IntegerVector S1,int r1,
                  arma::vec tau2, IntegerVector S2,int r2,
                  arma::vec gamma, arma::mat W){
  
  int N = X.size();
  int d = gamma.size();
  arma::vec mid_X = arma::zeros(N+predict_t);
  arma::vec predict_X = arma::zeros(N+predict_t);
  arma::vec epsi = arma::zeros(N+predict_t);
  
  for (int i=U;i<N;i++){
    mid_X[i] = X[i];
  }
  
  for (int i=U;i<N;i++){
    double X_lag = 0;
    for (int j=0;j<p;j++){
      X_lag = X_lag + phi[j]*X[i-j-1];
    }
    
    double epsi_lag = 0;
    for (int k=0;k<q;k++){
      epsi_lag = epsi_lag + psi[k]*epsi[i-k-1];
    }
    
    double X_sea = 0;
    for (int l=0;l<r1;l++){
      X_sea = X_sea + tau1[l]*X[i-S1[l]];
    }
    
    double epsi_sea = 0;
    for (int m=0;m<r2;m++){
      epsi_sea = epsi_sea + tau2[m]*epsi[i-S2[m]];
    }
    
    double exo = 0;
    for (int o=0;o<d;o++){
      exo = exo + gamma[o]*W(o,i);
    }
    
    epsi[i] = X[i] - X_lag - epsi_lag - X_sea - epsi_sea - exo;
    predict_X[i] = X_lag + epsi_lag + X_sea + epsi_sea + exo;
  }
  for (int i=N;i<(N+predict_t);i++){
    double X_lag = 0;
    for (int j=0;j<p;j++){
      X_lag = X_lag + phi[j]*mid_X[i-j-1];
    }
    
    double epsi_lag = 0;
    for (int k=0;k<q;k++){
      epsi_lag = epsi_lag + psi[k]*epsi[i-k-1];
    }
    
    double X_sea = 0;
    for (int l=0;l<r1;l++){
      X_sea = X_sea + tau1[l]*mid_X[i-S1[l]];
    }
    
    double epsi_sea = 0;
    for (int m=0;m<r2;m++){
      epsi_sea = epsi_sea + tau2[m]*epsi[i-S2[m]];
    }
    
    double exo = 0;
    for (int o=0;o<d;o++){
      exo = exo + gamma[o]*W(o,i);
    }
    
    predict_X[i] = X_lag + epsi_lag + X_sea + epsi_sea + exo;
    mid_X[i] = predict_X[i];
  }
  
  return predict_X;
}





