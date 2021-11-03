# Builds an adjacency list from a list of sequences with defined metrics and t as threshold
# 

adjL=function(L, t=5, d="lcs"){
  require(parallel)
  require(stringdist)
  
  cl = makeCluster(getOption("cl.cores", 4))
  clusterExport(cl,varlist=c("L","t"), envir = environment())
  clusterEvalQ(cl, library("stringdist"))
  aL=parSapply(cl,L, function(p){
    l=stringdist(p,L[L!=p],method=d)
    l=L[L!=p][l<t]
    return(l)
  })
  stopCluster(cl)
  names(aL)=L
  return(aL)
}