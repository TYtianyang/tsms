garma.fit = function(X,U,p,q,S1,S2,W=NULL,k=2){
  return(FSARMA.fit(X,U,p,q,S1,S2,W,k))
}

garma.pred = function(obj,X,W=NULL,pred_t){
  return(FSARMA.pred(obj,X,W,pred_t))
}

garma.gen = function(n, U, p, q, S1 = NULL, S2 = NULL, phi = NULL, psi = NULL, tau1 = NULL, tau2 = NULL, noise = 1){
  return(genTS(n, U, p, q, S1, S2, phi, psi, tau1, tau2, noise))
}

MS = function(X,W=NULL,U,p_range,q_range,r=1,S=NULL,
                       blur.out=c(2,2),sar=T,sma=F,sfourier=T,order=1,pred_t=0,
                       k=2,level=0.05,multicore=T){
  return(HSARMA.auto(X,W,U,p_range,q_range,r,S,
                       blur.out,sar,sma,sfourier,order,pred_t,
                       k,level,multicore))
}
