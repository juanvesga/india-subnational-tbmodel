return_output<-function(fx,x,ref){
  
  i <- ref$i; s <- ref$s
  # Pre allocate list of results
  
  
  
  model <- fx(x)
  
  t<-model$t
  sln<-model$soln
  
  tmp <- rowSums(sln[,1:i$nstates])
  pop <- (tmp[1:length(tmp)-1]+tmp[2:length(tmp)])/2
  
  tmp    <- rowSums(sln[,s$hi])
  pop_hi <- (tmp[1:length(tmp)-1]+tmp[2:length(tmp)])/2
  
  tmp_n <- apply(sln[,i$aux$inc],2,diff)
  inc   <- 1e5*tmp_n/cbind(pop,pop_hi,pop)
  
  inc_all<-inc[,1]
  inc_hi <-inc[,2]
  inc_mdr<-inc[,3]
  
  tmp <- rowSums(sln[,s$prevalent])
  pre <- (tmp[1:length(tmp)-1]+tmp[2:length(tmp)])/2;
  prev= 1e5*(pre/pop)
  
  tmp <- rowSums(sln[,intersect(s$hi,s$prevalent)])
  pre <- (tmp[1:length(tmp)-1]+tmp[2:length(tmp)])/2;
  prev_hi= 1e5*(pre/pop_hi)
  
  tmp_n <- diff(sln[,i$aux$mort])
  mort  <- 1e5*tmp_n/pop
  
  tmp_n <- diff(sln[,i$aux$noti[1]])
  not  <- 1e5*tmp_n/pop
  
  
  TBnoti<-model$TBnoti
  
  
  
  # Preventive cascade
  sfin <- sln[nrow(sln),]
  
  
  #Arrays
  res<-list()
  res$popu<-pop
  res$inc_all<-inc_all
  res$inc_hi<-inc_hi
  res$inc_mdr<-inc_mdr
  res$prev_all<-prev
  res$prev_hi<-prev_hi
  res$mort_tb<-mort
  res$notif_all<-not
  res$sfin<-sfin
  
  res$llk<-model$llk
  res$prpe<-model$TBepi[6]
  res$pr_mdr<-model$TBepi[7]
  
  

  # filename<-paste("res/","output","_",location,"_",runtype,sep="")
  # saveRDS(res,filename)
  # 
  return(res)
}


