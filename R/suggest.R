suggest = function(X,order=NULL,level=0.05){

  spec = spectrum(X,log='no',plot=FALSE)
  spx = spec$freq
  spy = 2*spec$spec
  z.score = scale(spy,scale=TRUE)
  p.value = 2*(1-pnorm(abs(z.score)))
  spx = spx[p.value<=level]
  spy = spy[p.value<=level]
  S = round(1/spx)[order(spy,decreasing=T)[1:ifelse(is.null(order),length(spy),min(length(spy),order))]]

  return(sort(S))
}
