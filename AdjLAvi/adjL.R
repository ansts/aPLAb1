# Builds an adjacency list from a list of sequences with defined metrics and t as threshold
# 

adjL=function(L, t=5, d="lcs",subs=NULL){
  require(parallel)
  require(stringdist)
  if (length(subs)>0) s=subs else s=seq_along(L)
  cl = makeCluster(16)
  clusterExport(cl,varlist=c("L","t", "s"), envir = environment())
  clusterEvalQ(cl, library("stringdist"))
  aL=parSapply(cl, s, function(i){
    p=L[i]
    l=stringdist(p,L[L!=p],method=d, nthread=1)
    l=L[L!=p][l<t]
    return(l)
  })
  stopCluster(cl)
  names(aL)=L[s]
  return(aL)
}