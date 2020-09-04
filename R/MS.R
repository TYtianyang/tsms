MS = function(X,W=NULL,U,p_range,q_range,r=1,S=NULL,
                       blur.out=c(2,2),sar=T,sma=T,
                       sfourier=F,order=1,
                       sfactor=F,pred_t=0,
                       crit='BC',level=0.05,multicore=T){

  # step 1: prepare
  X = as.matrix(X)

  # step 2: suggest seasonality
  if (is.null(S)){
    S = suggest(X,level)
  }
  S = unique(S[(S>(max(p_range,q_range)+blur.out[2]))&(S<(U-blur.out[1]))])

  # step 3: prune seasonality
  S.nat = sort(S)
  S.refer = list()
  j=1
  for (k in 1:length(S.nat)){
    if (k==1){
      S.refer[[j]] = S.nat[k]
      j=j+1
    }else{
      S.temp = c()
      for (i in 1:length(S.refer[[j-1]])){
        S.temp = c(S.temp,c((S.refer[[j-1]][i]-sum(blur.out)):
                              S.refer[[j-1]][i]+sum(blur.out)))
      }
      if (S.nat[k]%in%S.temp){
        S.refer[[j-1]] = c(S.refer[[j-1]],S.nat[k])
      }else{
        S.refer[[j]] = S.nat[k]
        j = j+1
      }
    }
  }
  S.pick = list()
  for (k in 1:length(S.refer)){
    S.temp = c()
    for (j in 1:length(S.refer[[k]])){
      S.temp = sort(unique(c(S.temp,c((S.refer[[k]][j]-blur.out[1]):
                                        (S.refer[[k]][j]+blur.out[2])))))
    }
    S.pick[[k]] = S.temp
  }
  r.nomial = min(r,length(S.pick))

  # step 4: merge W.list
  W.list = list(exo = W)
  if (sfourier|sfactor){
    for (i in 1:length(S.pick)){
      S = mean(S.pick[[i]])
      sea = matrix(nrow=as.numeric(sfourier)*order*2+as.numeric(sfactor)*S,ncol=nrow(X)+pred_t)
      if (sfourier){
        for (j in 1:order){
          sea[j,] = sin(2*pi*c(1:(nrow(X)+pred_t))/(S/j))
          sea[j+order,] = cos(2*pi*c(1:(nrow(X)+pred_t))/(S/j))
        }
      }
      if (sfactor){
        start = as.numeric(sfourier)*order*2
        base.vec = rep(0,nrow(X)+pred_t)
        base.vec[1+S*(c(1:floor((nrow(X)+pred_t)/S))-1)] = 1
        for (j in 1:S){
          sea[start+j,] = c(rep(0,j-1),base.vec[1:(nrow(X)+pred_t-j+1)])
        }
      }
      W.list[[paste(c('sea',i),collapse='')]] = sea
    }
  }

  # step 5: construct para_panel
  S_n = ncol(combn(c(1:length(S.pick)),r.nomial))
  para_panel = data.frame(matrix(nrow=length(p_range)*length(q_range)*S_n,ncol=7))
  colnames(para_panel) = c('p','q','S1','S2','W','ic','dim')

  comb = combn(c(1:length(S.pick)),r.nomial)
  S_vec = c()
  W_vec = c()
  for (i in 1:S_n){
    S_vec[i] = paste(unlist(S.pick[as.vector(comb[,i])]),collapse=' ')
    W_vec[i] = paste(c(paste('sea',as.vector(comb[,i]),sep=''),'exo'),collapse=' ')
  }
  para_panel[,5] = rep(W_vec,length(p_range)*length(q_range))
  if (sar){
    para_panel[,3] = rep(S_vec,length(p_range)*length(q_range))
  }else{
    para_panel[,3] = NA
  }
  if (sma){
    para_panel[,4] = rep(S_vec,length(p_range)*length(q_range))
  }else{
    para_panel[,4] = NA
  }
  para_panel[,1] = rep(p_range,rep(S_n*length(q_range),length(p_range)))
  para_panel[,2] = rep(rep(q_range,rep(S_n,length(q_range))),length(p_range))

  # step 6: fit and select optimal
  para_panel = parallel.fit(X,W.list,para_panel,U,'AIC',multicore)

  if (crit=='BC'){
    dim_up = para_panel[which.min(para_panel[,6]),'dim']
    for (i in 1:nrow(para_panel)){
      dim = para_panel[i,'dim']
      para_panel[i,'ic'] = para_panel[i,'ic'] - 2*dim + (nrow(X)^(2/3))*(sum(1/c(1:dim)))
    }
    temp_para_panel = para_panel[para_panel$dim<=dim_up,]
    best.row = temp_para_panel[which.min(temp_para_panel[,6]),1:5]
  }else{
    if (crit=='BIC'){
      para_panel$ic = para_panel$ic - 2*para_panel$dim + log(nrow(X))*para_panel$dim
    }
    best.row = para_panel[which.min(para_panel[,6]),1:5]
  }

  p = best.row[1,1]
  q = best.row[1,2]
  if (is.na(best.row[1,3])){S1=NULL
  }else{S1 = as.integer(strsplit(best.row[1,3],' ')[[1]])}
  if (is.na(best.row[1,4])){S2=NULL
  }else{S2 = as.integer(strsplit(best.row[1,4],' ')[[1]])}
  W.temp = rbindlist(W.list,best.row[1,5])
  best.m = garma.fit(X,U,p,q,S1,S2,W.temp,crit)
  prediction = garma.pred(best.m,X,W.temp,pred_t)

  output = list(obj=best.m,ic_panel=para_panel,prediction=prediction)
}
