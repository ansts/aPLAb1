# Generates ROC curves to compare different string metrics (editing distances) 
# applied to the original AA sequences as well as to recoded AA sequences 
# to emphasize AA pairs of common properties, as well as pairwise alignment data

require(stringdist)
require(pROC)
require(Biostrings)
require(parallel)
require(matrixStats)

mims=read.table(file="mims.txt", quote = NULL, header = F,row.names = NULL)
mims=mims[,1]
mthds=c("jaccard","cosine","hamming","qgram","lv","dl","osa","lcs")

rjac=rocBootstrap(mthd="jaccard", q=2)
plot(rjac[[1]])
for (i in 2:20){
    plot(rjac[[i]], add=T)
}
x=sapply(rjac,function(cu){
  y=cu$auc
  return(y)
})
aucBxp=x

rcos=rocBootstrap(mthd="cosine", q=2)
#plot(rcos[[1]])
for (i in 2:20){
  plot(rcos[[i]], add=T, col=2)
}
x=sapply(rcos,function(cu){
  y=cu$auc
  return(y)
})
aucBxp=cbind(aucBxp,cos=x)
colnames(aucBxp)[1]="jac"

rham=rocBootstrap(mthd="hamming")             
#plot(rham[[1]])
for (i in 2:20){
  plot(rham[[i]], add=T, col=3)
}
x=sapply(rham,function(cu){
  y=cu$auc
  return(y)
})
aucBxp=cbind(aucBxp,ham=x)

rqgr=rocBootstrap(mthd="qgram",q=2)
#plot(rqgr[[1]])
for (i in 2:20){
  plot(rqgr[[i]], add=T, col=4)
}
x=sapply(rqgr,function(cu){               
  y=cu$auc
  return(y)
})
aucBxp=cbind(aucBxp,qgr=x)

rlv=rocBootstrap(mthd="lv")              
#plot(rlv[[1]])
for (i in 2:20){
  plot(rlv[[i]], add=T, col=5)
}
x=sapply(rlv,function(cu){
  y=cu$auc
  return(y)
})
aucBxp=cbind(aucBxp,lv=x)

rdl=rocBootstrap(mthd="dl")                 
#plot(rdl[[1]])
for (i in 2:20){
  plot(rdl[[i]], add=T, col=6)
}
x=sapply(rdl,function(cu){
  y=cu$auc
  return(y)
})
aucBxp=cbind(aucBxp,dl=x)


rosa=rocBootstrap(mthd="osa")            
for (i in 1:20){
  plot(rosa[[i]], add=T, col=7)
}
x=sapply(rosa,function(cu){
  y=cu$auc
  return(y)
})
aucBxp=cbind(aucBxp,osa=x)

rlcs=rocBootstrap(mthd="lcs")             
for (i in 1:20){
  plot(rlcs[[i]], add=T,col=8)
}
x=sapply(rlcs,function(cu){
  y=cu$auc
  print(y)
  return(y)        
})
aucBxp=cbind(aucBxp,lcs=x)
legend("bottomright",legend=mthds, fill = 1:8)

rocs=list(rjac, rcos,rham,rqgr,rlv,rdl,rosa,rlcs)
roctests=sapply(rocs, function(ri){
  sapply(rocs, function(rj){
    roc.test(ri[[1]],rj[[1]])$p.value
  })
})
x=p.adjust(roctests)
roctests=array(x, dim = dim(roctests))
rownames(roctests)=c("rjac","rcos","rham","rqgr","rlv","rdl","rosa","rlcs")
colnames(roctests)=c("rjac","rcos","rham","rqgr","rlv","rdl","rosa","rlcs")
boxplot(aucBxp, notch=T, main="AUC Comparison")

SpSe=sapply(rocs,function(r){
  x=sapply(r, function(y){
    x=coords(y,x="best")
    if (!is.null(dim(x))) x=x[,1]
    return(x)
  })
  m=rowMeans(x)
  s=rowSds(x[2:3,])
  return(c(m,s))
})
colnames(SpSe)=c("rjac","rcos","rham","rqgr","rlv","rdl","rosa","rlcs")


aafmx=matrix("X",20,20)
rownames(aafmx)=aa
colnames(aafmx)=aa
aafm=list(M=c("M","I","L","V"), "F"=c("F","W","Y"), Y=c("Y","H"), D=c("D","E"), N=c("D","N"), H=c("H","N"), S=c("S","N"), E=c("E","Q"), K=c("E","K"), R=c("R", "K","Q"), A=c("S","A"), T=c("T","S"), G="G", C="C",P="P")
for (i in names(aafm)){
    aafmx[aafm[i][[1]],aafm[i][[1]]]=i
}

cl = makeCluster(getOption("cl.cores", 4))
clusterExport(cl,varlist=c("L7c", "mims", "daafm"), envir = environment())
clusterEvalQ(cl, library("pROC"))

rdaa=lapply(mthds, function(mthd){
  proct=proc.time()
  print(mthd)
  roci=parLapply(cl, 1:20, function(i){
  ho=c()
  he=c()
  for (i in seq_along(L7c[,1])) {
  s=seq_along(L7c[,1])
  n=0
  for (j in s[s!=i]) {
    if (L7c[i,2]==L7c[j,2]) {
      d=daafm(as.character(L7c[i,1]),as.character(L7c[j,1]),mth = mthd)
      ho=c(ho,d)
      n=n+1
    }
  }
  for (j in 1:n){
    p=sample(mims,1)
    d=daafm(as.character(L7c[i,1]),p,mth = mthd)
    he=c(he,d)
    }
  }

  hohe=cbind(c(ho,he),c(rep(1,length(ho)),rep(0,length(he))))
  x=roc(hohe[,2], hohe[,1])
  return(x)
 })
  print(proc.time()-proct)
  return(roci)
})

stopCluster(cl)
names(rdaa)=mthds
aucrdaa=sapply(seq_along(rdaa), function(j){
  ri=rdaa[[j]]
  if (i==1) plot(ri[[1]])
  for (i in 2:20){
    plot(ri[[i]], add=T, col=j)
  }
  sapply(ri,function(cu){
    return(cu$auc)
  })
})
legend("bottomright",legend=mthds, fill = 1:8)
rdaamthds=paste(mthds,"_daa", sep="")
colnames(aucrdaa)=rdaamthds
bestauc=apply(cbind(aucBxp,aucrdaa),2, which.max)
bauc=apply(cbind(aucBxp,aucrdaa),2,function(co){
  max(co)
})

for (j in mthds){
for (i in 1:20){
  if (j=="jaccard"&i==1) plot(rdaa[j][[1]][[i]]) else plot(rdaa[j][[1]][[i]], add=T)
}
}

x=sapply(rdaa,function(m){
  sapply(m, function(cu){
  y=cu$auc
  return(y)
})
})
aucBxp2=cbind(aucBxp,x)
boxplot(aucBxp2, notch=T)
################################################################################
data(BLOSUM62)
cl = makeCluster(getOption("cl.cores", 4))
clusterExport(cl,varlist=c("L7c", "mims", "daafm","BLOSUM62"), envir = environment())
clusterEvalQ(cl, library("Biostrings"))
clusterEvalQ(cl, library("pROC"))

rpwa=parLapply(cl, 1:20, function(i){
  print(i)
  ho=c()
  he=c()
  for (i in seq_along(L7c[,1])) {
  s=seq_along(L7c[,1])
  n=0
  for (j in s[s!=i]) {
    if (L7c[i,2]==L7c[j,2]) {
      d=pairwiseAlignment(as.character(L7c[i,1]),as.character(L7c[j,1]), substitutionMatrix=BLOSUM62)@score
      ho=c(ho,d)
      n=n+1
    }
  }
  for (j in 1:n){
    p=sample(mims,1)
    d=pairwiseAlignment(as.character(L7c[i,1]),p, substitutionMatrix=BLOSUM62)@score
    he=c(he,d)
  }
  }

  hohepwa=cbind(c(ho,he),c(rep(1,length(ho)),rep(0,length(he))))
  x=roc(hohepwa[,2], hohepwa[,1])
  return(x)
})
stopCluster(cl)

plot(rpwa[[1]])
for (i in 2:20){
  plot(rpwa[[i]], add=T)
}
x=sapply(rpwa,function(cu){
  y=cu$auc
  return(y)
})
aucBxp3=cbind(aucBxp2,pwa=x)

aucBxp=aucBxp[,c(1,2,6,7,8,5,4,9,10,3)]
boxplot(aucBxp3, notch=T)

plot(rjac[[1]])
for (i in 2:20){
  plot(rjac[[i]], add=T)
}
for (i in 1:20){
  plot(rcos[[i]], add=T, col=2)
}
plot(rdl[[1]])
for (i in 1:20){
  plot(rdl[[i]], add=T, col=3)
}
for (i in 1:20){
  plot(rosa[[i]], add=T, col=4)
}
for (i in 1:20){
  plot(rlcs[[i]], add=T, col=5)
}
for (i in 1:20){
  plot(rlv[[i]], add=T, col=6)
}
for (i in 1:20){
  plot(rqgr[[i]], add=T, col=8)
}
for (i in 1:20){
  plot(rpwa[[i]], add=T, col=9)
}
for (i in 1:20){
  plot(rdaa[[i]], add=T, col=10)
}


legend("bottomright",legend=c("dl","lcs","aa alignment"), fill=c(3,5,2))
