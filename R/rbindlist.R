rbindlist = function(W.list,W.id){
  ids = strsplit(W.id,split=' ')[[1]]
  W = W.list[[ids[1]]]
  if (length(ids)>1){
    for (i in 2:length(ids)){
      W = rbind(W,W.list[[ids[i]]])
    }
  }
  return(W)
}