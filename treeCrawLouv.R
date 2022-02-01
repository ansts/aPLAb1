treeCrawLouv=function(L, G, s=1000, r=NULL, nn="C"){
  require(igraph)
  require(stringi)
  
  g=induced.subgraph(G,L)
  
  cl=cluster_louvain(g)$memberships
  if (is.null(r)) nom=nrow(cl) else nom=r     
  cl=cl[nom,]
  if (max(cl)==1) {
    return("atomic")
  }
  cl=aggregate(L,by=list(cl), "list")$x
  l=lengths(cl)
  ni=paste(nn,seq_along(cl), sep="_")
  print(ni)
  x=lapply(seq_along(cl),function(i){
    co=cl[i]
    if (length(co[[1]])>s) {              
        nii=ni[i]
        y=treeCrawLouv(co[[1]],g, s=s, nn=nii)
        return(y)
      }
      else {
        names(co)=ni[i]
        co=lapply(co[[1]],function(y) list(Seq=y))
        names(co)=paste(ni[i],seq_along(co), sep="_")
        # }
        #  print(names(co))
        return(co)
      }
  })
  names(x)=ni
  return(x)
}
