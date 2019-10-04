moving = function(X,width){
  if (width %% 2 == 0){message('Your width is even number!')}
  else{
    trend = filter(X,rep(1/width,width),sides=2)
    trend[1:((width-1)/2)] = trend[((width-1)/2)+1]
    trend[(nrow(X)-((width-1)/2)+1):nrow(X)] = trend[nrow(X)-((width-1)/2)]
    return(trend)
  }
}