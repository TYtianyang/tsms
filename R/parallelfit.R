parallel.fit = function(X,W.list,para_panel,U,crit,multicore=F){

  if (!multicore){
    # single core
    for (i in 1:nrow(para_panel)){
      temp_row = para_panel[i,]
      p = temp_row[1,1]
      q = temp_row[1,2]
      if (is.na(temp_row[1,3])){S1=NULL
      }else{S1 = as.integer(strsplit(temp_row[1,3],' ')[[1]])}
      if (is.na(temp_row[1,4])){S2=NULL
      }else{S2 = as.integer(strsplit(temp_row[1,4],' ')[[1]])}
      W.temp = rbindlist(W.list,temp_row[1,5])

      m = FSARMA.fit(X,U,p,q,S1,S2,W.temp,crit)
      para_panel[i,6] = m$ic
      para_panel[i,7] = m$dim
    }

  }else{
    # multiple cores
    ncores = min(detectCores(),nrow(para_panel))
    cl = makePSOCKcluster(ncores)

    l = floor(nrow(para_panel)/ncores)
    assign_loc = rep(l,ncores)
    if (nrow(para_panel)%%ncores!=0){
      for (i in 1:(nrow(para_panel)-l*ncores)){
        assign_loc[i] = assign_loc[i] + 1
      }
    }
    assign_list = list()
    assign_loc_cums = cumsum(assign_loc)
    a = 1
    for (i in 1:ncores){
      assign_list[[i]] = c(a:assign_loc_cums[i])
      a = 1 + assign_loc_cums[i]
    }

    singleFit = function(core_ind,assign_list,para_panel){
      index_set = assign_list[[core_ind]]
      ic_set = rep(0,length(index_set))
      dim_set = rep(0,length(index_set))
      for (i in 1:length(index_set)){
        index = index_set[i]
        temp_row = para_panel[index,]
        p = temp_row[1,1]
        q = temp_row[1,2]
        if (is.na(temp_row[1,3])){S1=NULL
        }else{S1 = as.integer(strsplit(temp_row[1,3],' ')[[1]])}
        if (is.na(temp_row[1,4])){S2=NULL
        }else{S2 = as.integer(strsplit(temp_row[1,4],' ')[[1]])}
        W.temp = rbindlist(W.list,temp_row[1,5])

        m = FSARMA.fit(X,U,p,q,S1,S2,W.temp,crit)
        ic = m$ic
        dim = m$dim
        ic_set[i] = ic
        dim_set[i] = dim
      }
      output = list(ic_set=ic_set,dim_set=dim_set)
      return(output)
    }

    clusterExport(cl, list("FSARMA.fit","FSARMA.pred", "rbindlist","singleFit", "assign_list"
                           ,"para_panel","W.list","X","U","crit")
                  , envir=environment())
    pout = parLapply(cl, c(1:ncores), singleFit, assign_list=assign_list
                     ,para_panel=para_panel)
    ic.list = list()
    dim.list = list()
    for (i in 1:length(pout)){
      ic.list[[i]] = pout[[i]][[1]]
      dim.list[[i]] = pout[[i]][[2]]
    }
    ic_col = unlist(ic.list)
    dim_col = unlist(dim.list)
    para_panel[,6] = ic_col
    para_panel[,7] = dim_col
    stopCluster(cl)
  }

  return(para_panel)
}
