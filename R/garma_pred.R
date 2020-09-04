# freq predict
garma.pred = function(obj,X,W=NULL,pred_t){
  
  if (dim(X)==1){
    X = as.matrix(X) 
  }
  
  # correct pred_t & W
  if (is.null(W)){
    W = t(as.matrix(rep(0,nrow(X)+pred_t)))
  }
  
  pred_t = min(pred_t,ncol(W)-ncol(X))
  
  # init arguments
  n = nrow(X)
  U = obj[[1]]
  p = obj[[2]]
  phi = obj[[3]]
  q = obj[[4]]
  psi = obj[[5]]
  r1 = obj[[6]]
  S1 = obj[[7]]
  tau1 = obj[[8]]
  r2 = obj[[9]]
  S2 = obj[[10]]
  tau2 = obj[[11]]
  gamma = obj[[12]]
  sigma = obj[[13]]
  
  pred_X = Pred(X,U,pred_t,
                   as.double(phi),p,
                   as.double(psi),q,
                   as.double(tau1), as.integer(S1),r1,
                   as.double(tau2), as.integer(S2),r2,
                   as.double(gamma), W)
  
  return(pred_X)
}
