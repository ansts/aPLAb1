# Function generating 20 ROC curves testing a predictor of 
# single mAb mimotope pairs based on a given string metric
# by comparing pairs of known mimotopes to the same mabsim
# and pairs of any one of those with irrelevant mimotope 
# from our glbal IgM library.

rocBootstrap=function(m=mims,L7=L7c,n=20, mthd="hamming", q=1){
  require(pROC)
  require(stringdist)
  require(parallel)
  proct=proc.time()
  
  cl = makeCluster(getOption("cl.cores", 4))
  clusterExport(cl,varlist=c("L7", "m", "mthd"), envir = environment())
  clusterEvalQ(cl, library("stringdist"))
  clusterEvalQ(cl, library("pROC"))
  
  BS=parLapply(cl, 1:n, function(i){
    #print(i)
    ho=c()
    he=c()
    for (i in seq_along(L7[,1])) {
      s=seq_along(L7[,1])
      n=0
      for (j in s[s!=i]) {
        if (L7[i,2]==L7[j,2]) {
          d=stringdist(L7[i,1],L7[j,1], method = mthd,q=q)
          ho=c(ho,d)
          n=n+1
        }
      }
      for (j in 1:n){
        p=sample(m,1)
        d=stringdist(L7[i,1],p, method = mthd,q=q)
        he=c(he,d)
      }
    }
    
    hohe=cbind(c(ho,he),c(rep(1,length(ho)),rep(0,length(he))))
    
    r=roc(hohe[,2], hohe[,1])
    return(r)
  })
  stopCluster(cl)
  print(proc.time()-proct)
  return(BS)
}