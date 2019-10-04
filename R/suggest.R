suggest = function(X,level=0.05){
  
  spec = spectrum(X,log='no',plot=FALSE)
  spx = spec$freq
  spy = 2*spec$spec
  z.score = scale(spy,scale=TRUE)
  p.value = 2*(1-pnorm(abs(z.score)))
  spx = spx[p.value<=level]
  spy = spy[p.value<=level]
  fmin = min(diff(spec$freq))
  S = round(1/spx)
  S = unique(S[(S>(max(c(p_range,q_range))+max(blur.out)))&(S<(U-max(blur.out)))])
  
  ind = diff(1/S)<(2*fmin)
  if (sum(ind)>0){
    ind_l = S[c(ind,FALSE)]
    ind_r = S[c(FALSE,ind)]
    S = S[!c(ind,FALSE)]
    ind_mat = matrix(nrow=length(ind_l),ncol=2)
    ind_mat[,2] = ind_l
    ind_mat[,1] = ind_r
    
    cut = 2000
    cl = c()
    cr = c()
    for (i in 1:nrow(ind_mat)){
      row.temp = seq(from=ind_mat[i,1],to=ind_mat[i,2],by=cut)
      if (length(row.temp)==1){
        cl = c(ind_mat[i,1],cl)
        cr = c(ind_mat[i,2],cr)
      }else{
        cl = c(ind_mat[i,1],row.temp[2:length(row.temp)],cl)
        cr = c(row.temp[2:length(row.temp)],ind_mat[i,2],cr)
      }
    }
    
    ind_mat = matrix(nrow=length(cl),ncol=2)
    ind_mat[,1] = sort(cl,decreasing=T)
    ind_mat[,2] = sort(cr,decreasing=T)
    
    vp = function(P,Y){
      m = lm(Y~matrix(c(sin(2*pi*c(1:nrow(Y))/P),cos(2*pi*c(1:nrow(Y))/P)),ncol=2))
      cf = as.vector(m$coefficients[-1])
      p = sum(cf^2)
      return(p)
    }
    
    S.add = rep(0,nrow(ind_mat))
    for (i in 1:nrow(ind_mat)){
      range = c(ind_mat[i,1]:ind_mat[i,2])
      p.list = apply(as.matrix(range),1,FUN=vp,Y=X)
      S.add[i] = range[which.max(p.list)]
    }
    
    S = unique(c(S,S.add))
  }
  return(sort(S))
}