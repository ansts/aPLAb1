ngbtbl <- function(AL, fl) {
  require(parallel)
  t0=rep(0,7)
  names(t0)=1:7
  
  cl=makeCluster(16)
  clusterExport(cl, c("fl", "t0"), envir = environment())
  Ngb=t(parSapply(cl,AL, function(l){
    lf=fl[l]
    t1=table(lf)
    t0[names(t1)]=t1
    return(t0)
  }))
  stopCluster(cl)
  return(Ngb)
}