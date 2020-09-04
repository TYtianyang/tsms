garma.gen = function(n,U,p,q,S1=NULL,S2=NULL,phi=NULL,psi=NULL,tau1=NULL,tau2=NULL,
                 noise=1){
  if (length(phi)==0){
    rand_phi = T
  }else{rand_phi=F}
  if (length(psi)==0){
    rand_psi = T
  }else{rand_psi=F}
  if (length(tau1)==0){
    rand_tau1 = T
  }else{rand_tau1=F}
  if (length(tau2)==0){
    rand_tau2 = T
  }else{rand_tau2=F}
  
  stop = 100
  while (stop>=99){
    X = matrix(0,n,1)
    epsi = matrix(rnorm(n,sd=noise),n,1)
    X[1:(U-1),1] = rnorm(U-1)
    
    if (rand_phi){
      for (i in 1:p){
        if (i==1){
          phi[i] = 0.5
        }else{
          phi[i] = runif(1,-1,1)*phi[i-1]
        }
      }
    }
    if (rand_psi){
      psi = rnorm(q,sd=0.3)
    }
    if (rand_tau1){
      tau1 = rep(0.2,length(S1))
    }
    if (rand_tau2){
      tau2 = rep(0.2,length(S2))
    }
    
    beta = matrix(0,p+q+length(S1)+length(S2),1)
    beta[1:p,1] = phi
    beta[(p+1):(p+q),1] = psi
    if (length(S1)!=0){
      beta[(p+q+1):(p+q+length(S1)),1] = tau1
    }
    if (length(S2)!=0){
      beta[(p+q+length(S1)+1):(p+q+length(S1)+length(S2))] = tau2
    }
    for (t in U:n){
      X_G = matrix(0,p+q+length(S1)+length(S2),1)
      X_G[1:p,1] = X[(t-1):(t-p),1]
      X_G[(p+1):(p+q),1] = epsi[(t-1):(t-q),1]
      if (length(S1)!=0){
        for (i in 1:length(S1)){
          X_G[p+q+i,1] = X[t-S1[i],1]
        }
      }
      if (length(S2)!=0){
        for (i in 1:length(S2)){
          X_G[p+q+length(S1)+i,1] = epsi[t-S2[i],1]
        }
      }
      X[t,1] = t(X_G)%*%beta + epsi[t,1]
    }
    stop = max(abs(X[(n-10):n,1]))
  }
  obj = list(X=X,p=p,q=q,S1=S1,S2=S2,beta=beta)
  return(obj)
}
