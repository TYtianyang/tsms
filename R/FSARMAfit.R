# freq fitting
FSARMA.fit =function(X,U,p,q,S1,S2,W=NULL,crit='AIC'){

  if (is.null(W)){
    W = t(as.matrix(rep(0,nrow(X))))
  }
  r1 = length(S1)
  r2 = length(S2)
  S1.int = as.integer(S1)
  S2.int = as.integer(S2)
  d = nrow(W)

  Rfit = function(beta){
    SS = fit(X,U,beta,p,q,S1.int,r1,S2.int,r2,W)
    return(SS)
  }
  Rgrad = function(beta){
    gradient = grad(X,U,beta,p,q,S1.int,r1,S2.int,r2,W)
    return(gradient)
  }

  ult = tryCatch({
    beta = rep(0,p+q+r1+r2+d)
    beta = optim(beta,Rfit,Rgrad,method='BFGS')$par
    if (p==0){phi=NULL
    }else{phi=beta[1:p]}
    if (q==0){psi=NULL
    }else{psi=beta[(p+1):(p+q)]}
    if (r1==0){tau1=NULL
    }else{tau1=beta[(p+q+1):(p+q+r1)]}
    if (r2==0){tau2=NULL
    }else{tau2=beta[(p+q+r1+1):(p+q+r1+r2)]}
    if (is.null(d)){gamma=NULL
    }else{gamma = beta[(p+q+r1+r2+1):(p+q+r1+r2+d)]}
    sigma = sqrt(Rfit(beta)/(nrow(X)-U))

    mdl = likeli(X,U,
                 as.double(phi),p,
                 as.double(psi),q,
                 as.double(tau1),S1.int,r1,
                 as.double(tau2),S2.int,r2,
                 as.double(gamma), W, as.double(sigma))
    dim = (p+q+length(S1.int)+length(S2.int)+ifelse(is.null(W),0,d))

    if (crit=='AIC'){
      ic = - 2*mdl + 2*dim
    }

    if (crit=='BIC'){
      ic = - 2*mdl + log(nrow(X))*dim
    }

    if (crit=='BC'){
      ic = -2*mdl + (nrow(X)^(2/3))*(sum(1/c(1:dim)))
    }

    ult = list(U,p,phi,q,psi,r1,S1,tau1,r2,S2,tau2,gamma,sigma,dim,ic)
  },error=function(cond){
    if (p>0){phi=rep(0,p)
    }else{phi=NULL}
    if (q>0){psi=rep(0,q)
    }else{psi=NULL}
    if (r1>0){tau1=rep(0,r1)
    }else{tau1=NULL}
    if (r2>0){tau2=rep(0,r2)
    }else{tau2=NULL}
    if (!is.null(d)){gamma=rep(0,d)
    }else{gamma=NULL}
    dim = (p+q+length(S1.int)+length(S2.int)+ifelse(is.null(W),0,d))
    sigma = 1
    ic = Inf
    ult = list(U,p,phi,q,psi,r1,S1,tau1,r2,S2,tau2,gamma,sigma,dim,ic)
    # message(paste(c("Failure of fitting:    ",
    #                 'p=',p,' ','q=',q,' ',
    #                 'S1=c(',paste(S1,collapse=','),') ',
    #                 'S2=c(',paste(S2,collapse=','),')')))
    return(ult)
  })

  names(ult) = c('U','p','phi','q','psi','r1','S1','tau1','r2','S2','tau2','gamma',
                 'sigma','dim','ic')
  return(ult)

}

